#   Design a protein using MPNN and other intermediate steps:
#   Produce a PSSM matrix to ensure that conserved areas of the protein aren't altered
#   Remove signal peptide (if present)
#   Run MPNN to obtain individual mutations
#   Obtain top five mutations 
#   Iteratively add five additional mutations 


from colabfold.colabfold import run_mmseqs2
from Bio import SeqIO, AlignIO
from colabdesign.mpnn import mk_mpnn_model
from colabdesign.mpnn.model import aa_order
from io import StringIO
import os
from os import listdir
import shutil
import json
from Bio import SeqIO
import subprocess
from collections import defaultdict
import pandas as pd
import time
from multiprocessing.pool import ThreadPool
import make_pssm_dict as prep_mpnn
from tqdm import tqdm




class runMPNN():
    
    def __init__(self, pdb_path: str):
        self.target = ""
        self.name = pdb_path.split("/")[-1].replace(".pdb","")
        self.native = ""
        self.seq_len = 0
        
        
    def run(self):
            
        # create a temp file for this tool that will be deleted later
        os.makedirs(".temp", exist_ok=True)
        
        # clean pdb and get native sequence and length
        print(f"python rosetta_helper/clean_pdb.py example_run/{self.name}.pdb A")
        subprocess.run(f"python rosetta_helper/clean_pdb.py example_run/{self.name}.pdb A", shell=True)
        
        # update native sequence and length
        with open(f"{self.name}_A.fasta", "r") as fasta:
            self.native = fasta.readlines()[-1].replace("\n","")
            self.seq_len = len(self.native)
            self.target = f"{self.name}_A.pdb"
            self.name = self.name + "_A"

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
        self._getPSSM(alignment=alignment)
        status = self._runMPNN()
        print(status)

        # obtain mutation sequences for top 10 mutations
        data = self._getSeqs()
        
        os.makedirs("results", exist_ok=True)
        data.to_csv(f"results/{self.name}.csv", index=False)
        
        # # delete temp folder and MSA folder
        # shutil.rmtree("./.temp")
        shutil.rmtree(f"./{self.name}_env")
        # TODO: clean up .fasta 'n .pdb

        return "MPNN has successfully finished making designed sequences."

    
    def _getPSSM(self, alignment):
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
        query = f"{self.name}.fasta"

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
        mpnn_model.prep_inputs(pdb_filename=self.target, 
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



    def _runMPNN(self):
        """
        Run helper scripts and protein MPNN
        """
        
        os.makedirs(".temp/csv", exist_ok=True)

        # Run helper scripts to set up for MPNN
        subprocess.run(f"python ProteinMPNN/helper_scripts/parse_multiple_chains.py \
                        --input_path ./ \
                        --output_path .temp/parsed_pdbs.jsonl", shell=True)
        
        prep_mpnn.make_dict(self.seq_len)
        
        start = time.time()

        # Run MPNN on sequence and mutating +_scoring one by one
        print("Mutating positions one by one...")
        # Set the number of threads to the number of cores on your machine -1  
        pool = ThreadPool()
        tasks = range(self.seq_len)
        
        # Wrap the task range with tqdm for progress bar
        for _ in tqdm(pool.imap_unordered(self._task, tasks), total=len(tasks)):
            pass

        pool.close()
        pool.join()
        
        end = time.time()
        
        return f"MPNN has finished running in {end-start} seconds, \
            and all data has been recorded in .temp/csv"

   
   
    def _getSeqs(self):
        
        all_data = {
            "Position": [],
            "Sequence": [],
            "Wildtype": [],
            "Mutant": [],
            "Score": []
        }
        
        # iterate through the csv directory
        filenames = listdir(".temp/csv")
        files = [filename for filename in filenames if filename.endswith(".csv")]
        
        for file in files:
            all_data["Position"].append(int(file.replace(".csv","").split("_")[-1]) + 1)
            with open(".temp/csv/"+file, "r") as f:
                data = f.readlines()[0].split(",")
                all_data["Sequence"].append(data[0])
                all_data["Wildtype"].append(data[1].split(" ")[0])
                all_data["Mutant"].append(data[1].split(" ")[-1])
                all_data["Score"].append(float(data[2].replace("\n","")))
        
        df = pd.DataFrame.from_dict(all_data)
        df.sort_values("Score", ascending=False, inplace=True,ignore_index=True)
        
        return df.head(10)
    #     # create sequence with top 5 scoring mutations
    #     muts = mut.head(5)
    #     to_mutate = dict(zip(muts['position'], muts['mutated_resi']))
    #     mpnn_seqs = []
    #     for i in range(len(native)):
    #         if (i+1) in to_mutate:
    #             mut_seq[i] = to_mutate[i+1]
    #     mpnn_seqs.append(''.join(mut_seq))
        
    #     # iteratively add next five mutations
    #     iterative_muts = mut.tail(5)
    #     to_mutate = dict(zip(iterative_muts['position'], iterative_muts['mutated_resi']))
    #     for key in to_mutate:
    #         mut_seq[key-1] = to_mutate[key]
    #         mpnn_seqs.append(''.join(mut_seq))
            
    #     return mpnn_seqs


x = runMPNN("example_run/7nei.pdb")
x.run()


# benchmark
#TODO: install streamlit interface after benchmarking
