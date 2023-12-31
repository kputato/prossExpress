#   Design a protein using MPNN and other intermediate steps:
#   Produce a PSSM matrix to ensure that conserved areas of the protein aren't altered
#   Remove signal peptide (if present)
#   Run MPNN to obtain individual mutations
#   Obtain top five mutations 
#   Iteratively add five additional mutations 

# ! NOTE: this whole file is zero indexed when talking about sequence positions!


from colabfold.colabfold import run_mmseqs2
from Bio import SeqIO, AlignIO
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.mpnn.model import aa_order

import os
from os import listdir
import json
import subprocess
import pandas as pd
import time
import shutil
import requests
from io import StringIO
from collections import defaultdict
from multiprocessing.pool import ThreadPool

import make_pssm_dict as prep_mpnn



class runMPNN():
    
    def __init__(self, pdb_id: str):
        self.name = pdb_id.lower() 
        self.native = ""
        self.seq_len = 0
        
        
    def run(self):
            
        # create a temp file for this tool that will be deleted later
        os.makedirs(f"{self.name}", exist_ok=True)
        os.makedirs(f"{self.name}/results", exist_ok=True)
        os.makedirs(".temp", exist_ok=True) 
        
        
        # get fastas from rcsb database
        # self.get_fasta(self.name, f"example_run/{self.name}.fasta")
        # run fastas on colabfold alphafold2_batch googlecolab and place folded pdbs in example_run
        
        # update native sequence and length
        with open(f"benchmarked/{self.name}.fasta", "r") as fasta:
            self.native = fasta.readlines()[-1].replace("\n","")
            self.seq_len = len(self.native)

        # run sequence alignment using colabdesign api
        run_mmseqs2(self.native, self.name)

        # find and parse a3m file to write to string
        seq_str = ''
        records = SeqIO.parse(f"{self.name}_env/uniref.a3m", "fasta")

        for record in records:
            seq_str += (">" + str(record.id) + '\n')
            seq_str += (str(record.seq).replace("-","") + '\n')

        #run ClustalO on mmseq2 sequences
        child = subprocess.Popen(['mafft', '--quiet', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        child.stdin.write(seq_str.encode())
        child_out = child.communicate()[0].decode('utf8')
        alignment = AlignIO.read(StringIO(child_out), 'fasta')
        child.stdin.close()

        # take aligned target sequence, and get all positions with gaps
        # remove these gap positions from other sequences
        reference = alignment[0].seq
        to_del = []
        for i in range(len(reference)):
            if reference[i] == "-":
                to_del.append(i)
                
        for align in alignment:
            pos = {}
            for i in range(len(align.seq)):
                pos[i] = align.seq[i]
            for key in to_del:
                pos.pop(key)
            align.seq = "".join(list(pos.values()))
        
        # run MPNN and update json database
        status_one = self.get_pssm(alignment=alignment)
        print(status_one)
        status_two = self.run_mpnn()
        print(status_two)

        # obtain top 10 mutations
        data = self.get_muts()
        data.to_csv(f"{self.name}/results/{self.name}.csv", index=False)
        
        with open(f"{self.name}/results/{self.name}_sequences.json", "w") as f:
            json.dump(self.get_seqs(data), f)
        
        # delete temp folder and MSA folder
        shutil.rmtree("./.temp")
        shutil.rmtree(f"./{self.name}_env")

        return "MPNN has successfully finished making designed sequences."

    # TODO: get fasta from user given pdb file
    def get_fasta(self, pdb_id, output_file):
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for any HTTP error
        except requests.exceptions.RequestException as e:
            print(f"Failed to retrieve FASTA for PDB ID {pdb_id}: {e}")
            return
        
        fasta_text = response.text.strip()
        if fasta_text:
            with open(output_file, "w") as f:
                f.write(fasta_text)
        else:
            print(f"No FASTA data found for PDB ID {pdb_id}")

    def get_pssm(self, alignment):
        """Get a PSSM matrix by running psi blast against MSA object; process it for MPNN"""
        
        # Write MSA to FASTA file and turn it into a database
        db = f".temp/database.fasta"
        seq_records = [SeqIO.SeqRecord(seq=record.seq, id=record.id, description="") for record in alignment]
        tmp = open(db, 'w')
        for records in seq_records:
            tmp.write(">" + records.id + "\n")
            tmp.write(records.seq + "\n")
        tmp.close()
        subprocess.run("makeblastdb -in .temp/database.fasta -dbtype prot -out .temp/db", shell=True)

        # prepare for psi blast pssm
        query = f"benchmarked/{self.name}.fasta"

        # obtain a PSSM matrix using psi blast
        subprocess.run(f'psiblast -query {query} -db .temp/db -num_iterations 3 -out_ascii_pssm \
                       .temp/pssm.txt -outfmt 0', shell=True)
        
        return "PSSM matrix has been successfully produced and saved in .temp/pssm.txt"

    def _task(self, i):
        # multithreaded process of runing MPNN
        
        pos_data = defaultdict(list)
        
        # get PSSM matrix
        with open('.temp/pssm_dict.jsonl','r') as pssm:
            temp = json.load(pssm)
            bias = temp[self.name]["A"]["pssm_bias"]
            
        bias_matrix = bias[i]
        pos = i + 1
        
        if i < 1:
            fixed = f"{pos+1}-{len(bias)}"
        elif pos >= len(bias):
            fixed = f"1-{i}"
        elif i == 1:
            fixed = f"1, {pos+1}-{len(bias)}"
        elif ((pos+1) == len(bias)):
            fixed = f"1-{i},{len(bias)}"
        else:
            fixed = f"1-{i},{pos+1}-{len(bias)}"
            
        # create mpnn model
        mpnn_model = mk_mpnn_model()
        mpnn_model.prep_inputs(pdb_filename=f"{self.name}/{self.name}.pdb", 
                            fix_pos=fixed)
        
        # adjust PSSM probabilities
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        for aa in range(len(alphabet)):
            mpnn_model._inputs["bias"][i,aa_order[alphabet[aa]]] = bias_matrix[aa]
        
        samples = mpnn_model.sample_parallel(batch=5)
        
        # record the data
        ctr = 0
        for sequence in samples["seq"]:
            for j in range(self.seq_len):
                if self.native[j] != sequence[j]:
                    pos_data["sequence"].append(sequence)
                    mut = f"{self.native[j]}    {sequence[j]}" # record what native position has been mutated to
                    pos_data["mutation"].append(mut) 
                    pos_data["score"].append(samples["score"][ctr])
            ctr += 1   
        
        if pos_data:
            df = pd.DataFrame(pos_data)
            df['score'] = pd.to_numeric(df['score'])
            df2 = df.sort_values(by='score', ascending=False).reset_index(drop=True)
            df3 = df2.head(1)
            df3.to_csv(f".temp/csv/{self.name}_{i}.csv", index=False, header=False)
            
        return "One position has been designed, and the amino acid mutation that led \
            to the greatest improvement in thermostability has been recorded."

    def run_mpnn(self):
        """
        Run helper scripts and protein MPNN
        """
        
        os.makedirs(f".temp/csv", exist_ok=True)

        # Run helper scripts to set up for MPNN
        subprocess.run(f"python ProteinMPNN/helper_scripts/parse_multiple_chains.py \
                        --input_path {self.name} \
                        --output_path .temp/parsed_pdbs.jsonl", shell=True)
        
        prep_mpnn.make_dict(self.seq_len)
        
        start = time.time()

        # Run MPNN on sequence and mutating +_scoring one by one
        print("Mutating positions one by one...")
        pool = ThreadPool()
        pool.map(self._task, range(self.seq_len))
        pool.close()
        pool.join()
        
        end = time.time()
        
        return f"MPNN has finished running in {end-start} seconds, and all data has been recorded in .temp/csv."
   
    def get_muts(self):
        
        all_data = {
            "Position": [],
            "Sequence": [],
            "Wildtype": [],
            "Mutant": [],
            "Score": []
        }
        
        # iterate through the csv directory
        filenames = listdir(f".temp/csv")
        files = [filename for filename in filenames if filename.endswith(".csv")]
        
        for file in files:
            all_data["Position"].append(int(file.replace(".csv","").split("_")[-1]) + 1)
            with open(os.path.join(f".temp/csv", file), "r") as f:
                data = f.readlines()[0].split(",")
                all_data["Sequence"].append(data[0])
                all_data["Wildtype"].append(data[1].split(" ")[0])
                all_data["Mutant"].append(data[1].split(" ")[-1])
                all_data["Score"].append(float(data[2].replace("\n","")))
        
        df = pd.DataFrame.from_dict(all_data)
        df.sort_values("Score", ascending=False, inplace=True,ignore_index=True)
        
        return df.head(10)
    
    def get_seqs(self, data):
        print(data)
        # create sequence with the individual mutations
        to_mutate = dict(zip(data['Position'], data['Mutant']))
        mut_seq = list(self.native)
        for i in range(len(self.native)):
            if (i+1) in to_mutate:
                mut_seq[i] = f"<{to_mutate[i+1]}>"
        return "".join(mut_seq)




test_runs = ["4qum"]
for tests in test_runs:
    x = runMPNN(tests)
    x.run()

#TODO: install streamlit interface after benchmarking
