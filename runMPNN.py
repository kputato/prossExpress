# Design a protein using MPNN and other intermediate steps:
#   Produce a PSSM matrix to ensure that conserved areas of the protein aren't altered
#   Remove signal peptide (if present)
#   Run MPNN to obtain individual mutations
#   Obtain top five mutations 
#   Iteratively add five additional mutations 


from colabfold.colabfold import run_mmseqs2
from Bio import SeqIO, AlignIO
import subprocess
from io import StringIO
import os
import shutil
import numpy as np
import json
from Bio import SeqIO
import subprocess
from collections import defaultdict
import pandas as pd
import xml.etree.ElementTree as ET
import Bio.PDB as bpdb
import requests
from langchain.tools import BaseTool

class mpnnDesignTool():
    name = "mpnn_design"
    description = "This tool will be used after a PDB file is successfully extracted. We will first run MSA against the protein \
sequence to find conserved areas, resulting in a PSSM matrix representation. If the protein contains a signal peptide, the signal peptide \
will be removed. This tool will then iterate through the protein sequence and run MPNN, along with the PSSM matrix, on the protein, \
where each iteration will only allow MPNN to design a single specific position. The individual mutations are scored with MPNN and ranked \
and the top five mutations will be placed in the first designed sequence. The next five mutations will then be iteratively added, in addition \
to the first five mutatios, creating five more designed sequences. This tool is expected to place the designed sequences in a json file."

    def __str__(self):
        return f"{self.name}: {self.description}"
    
    def __init__(self, pdb_path: str):
        self.target = pdb_path
        self.name = pdb_path.split("/")[-1].replace(".pdb","")
        
    def run(self):
            
        # create a temp file for this tool that will be deleted later
        os.makedirs(".temp", exist_ok=True)
        
        # get sequence from pdb file and make fasta
        chain = {record.id: str(record.seq) for record in SeqIO.parse(self.target, 'pdb-seqres')}
        query = "".join(list(chain.values()))

        # run sequence alignment using colabdesign api
        run_mmseqs2(query, self.name)

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

        # find signal peptide as well
        start, end = self._get_peptide(self.name)

        for align in alignment:
            pos = {}
            for i in range(len(align.seq)):
                pos[i] = align.seq[i]
            for key in to_del:
                pos.pop(key)
            align.seq = "".join(list(pos.values()))
            align.seq = align.seq [:start-1] + align.seq[end:] # remove signal peptide

        # remove signal peptide from pdb
        self._process_pdb(pdb=self.target, start_resi=start, end_resi=end)

        # create a fasta file for the pdb
        f = open(f"./.temp/{self.name}.fasta", "w")
        with open(self.target, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                f.write('>wildtype\n')
                f.write(str(record.seq))
        f.close()

        # run MPNN and update json database
        self._getPSSM(alignment=alignment, uniprot=self.name)
        top_ten = self._runMPNN(len(str(record.seq)))

        # obtain mutation sequences for top 10 mutations
        print(self._getSeqs(top_ten))

        # delete temp folder and MSA folder
        shutil.rmtree("./.temp")
        shutil.rmtree(f"./{self.name}_env")

        return "MPNN has successfully finished making designed sequences. The designed sequences are \
available in the json file."

    def _get_peptide(self, uniprot_id):
    # Retrieve positions of the signal peptide inthe protein

        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"

        # Send a GET request to the UniProt API
        response = requests.get(url)

        # Check if the request was successful
        if response.ok:
            # Parse the XML response
            xml_data = response.content.decode()

            root = ET.fromstring(xml_data)

            # Find the feature elements
            feature_elements = root.findall(".//{http://uniprot.org/uniprot}feature")

            # Find the signal peptide feature
            for elements in feature_elements:
                if elements.get('type') == "signal peptide":
                    start = int(elements[0][0].attrib['position'])
                    end = int(elements[0][1].attrib['position'])
            return start, end  
        
    def _process_pdb(self, pdb, start_resi, end_resi):
        # modify pdb file in place to remove signal peptide

        s = bpdb.PDBParser().get_structure("tmp", pdb)

        # subclass that will tell pdb output class which residues to save
        class ResSelect(bpdb.Select):
            def accept_residue(self, res):
                if res.id[1] >= start_resi and res.id[1] <= end_resi:
                    return False
                else:
                    return True
                
        io = bpdb.PDBIO()
        io.set_structure(s)
        io.save(pdb, ResSelect())
        
        # renumber the pdb file to account for signal peptide removal
        io = bpdb.PDBIO()
        parser = bpdb.PDBParser()

        my_pdb_structure = parser.get_structure('test', 'example_run/A0A0K8P6T7.pdb')

        # renumber residue in my_pdb_structure
        residue_N = 1
        for model in my_pdb_structure:
            for chain in model:
                    for residue in chain:
                        print(residue.id)
                        if 'A' in residue.id[2]:
                            residue.id = (residue.id[0], residue_N-1, residue.id[2])       
                        else:
                            residue.id = (residue.id[0], residue_N, residue.id[2])
                            residue_N += 1
                
        io.set_structure(my_pdb_structure) 
        io.save('example_run/A0A0K8P6T7.pdb',  preserve_atom_numbering=True) 

        return None
    
    def _getPSSM(self, alignment, uniprot):
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
        query = f".temp/{uniprot}.fasta"

        # obtain a PSSM matrix using psi blast
        subprocess.run(f'psiblast -query {query} -db .temp/db -num_iterations 3 -out_ascii_pssm \
                       .temp/pssm.txt -outfmt 0', shell=True)

    def _runMPNN(self, seq_len):
        """
        Run helper scripts and protein MPNN
        """
        
        pdb_path = "/".join((self.target).split('/')[:-1])

        #Run helper scripts to set up for MPNN
        subprocess.run(f"python ProteinMPNN/helper_scripts/parse_multiple_chains.py \
                       --input_path {pdb_path} \
                        --output_path .temp/parsed_pdbs.jsonl", shell=True)
        
        subprocess.run("python make_pssm_dict.py", shell=True)
        
        # dictionaries to keep track of scores/sequences of mutated positions
        sequences = defaultdict(list)
        mutations = defaultdict(list)
        scores = {}
        
        # get native sequence for future comparison
        with open(f".temp/{self.name}.fasta", "r") as f:
                native = f.readlines()[1]
           
        # iterating through sequence and mutating +_scoring one by one
        for pos in range(1, seq_len+1):
            
            # fixed positions
            fixed = ""
            for i in range(1, seq_len+1):
                if i != pos:
                    fixed += f" {str(i)}"
                    
            # prepare MPNN 
            subprocess.run(f"python ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path .temp/parsed_pdbs.jsonl \
                --output_path .temp/fixed_pdbs.jsonl --chain_list A --position_list '{fixed}'", shell=True)
            
            # make sure you have the necessary files
            assert os.path.exists('.temp/parsed_pdbs.jsonl') 
            assert os.path.exists('.temp/pssm_dict.jsonl') 
            assert os.path.exists('.temp/fixed_pdbs.jsonl') 
                    
            # Run MPNN
            subprocess.run(f"python ProteinMPNN/protein_mpnn_run.py \
                                --jsonl_path .temp/parsed_pdbs.jsonl\
                                --out_folder ./.temp \
                                --sampling_temp '0.1' \
                                --pdb_path {self.target} \
                                --pssm_jsonl .temp/pssm_dict.jsonl\
                                --pssm_multi .8 \
                                --num_seq_per_target 32 \
                                --pssm_bias_flag 1 \
                                --pssm_log_odds_flag 1 \
                                --save_score 1 \
                                --fixed_positions_jsonl .temp/fixed_pdbs.jsonl", shell=True)
            
            os.remove(".temp/fixed_pdbs.jsonl")

            ctr = False # don't include native sequence 
            for record in SeqIO.parse(f"./.temp/seqs/{self.name}.fa", "fasta"):
                if ctr:
                    sequence = str(record.seq)
                    # ensure you are only recording sequences that have been mutated or are not empty
                    for i in range(len(native)):
                        if native[i] != sequence[i]:
                            sequences[pos].append(sequence)
                            
                            # record what native position has been mutated to
                            mut = f"{native[i]}    {sequence[i]}" 
                            mutations[pos].append(mut) 
                            
                            # record scores
                            s = np.load(f"./.temp/scores/{self.name}.npz")
                            scores[pos] = (s['score'].astype(str)).tolist() 
                ctr = True

        # write these to a pandas dataframe    
        # Create a DataFrame
        data = []
        for key in sequences.keys():
            sequence_list = sequences[key]
            score_list = scores[key]
            amino_acid_list = mutations[key]

            # Determine the length of the longest list
            max_length = max(len(sequence_list), len(score_list), len(amino_acid_list))

            for i in range(max_length):
                sequence = sequence_list[i] if i < len(sequence_list) else ''
                score = score_list[i] if i < len(score_list) else ''
                amino_acid = amino_acid_list[i] if i < len(amino_acid_list) else ''

                # Append the data as a dictionary to the list
                data.append({
                    'position': key,
                    'sequence': sequence,
                    'scores': score,
                    'wt_resi': amino_acid[0] if amino_acid else '',
                    'mutated_resi': amino_acid[-1] if amino_acid else ''
                })

        
        os.makedirs("results", exist_ok=True)
        
        df = pd.DataFrame(data)
        df.to_csv(f'results/{self.name}_mpnn.csv',index=False)
        df2 = pd.read_csv(f"results/{self.name}_mpnn.csv")
        df2 = df2.dropna()
        
        # Rank mutations based on the score in descending order
        df2['scores'] = pd.to_numeric(df2['scores'])  # Convert 'scores' column to numeric type
        df3 = df2.sort_values(by='scores', ascending=False).reset_index(drop=True)
        df3 = df2.drop_duplicates('position', keep='first')
        df3.to_csv(f'results/{self.name}_mpnn.csv',index=False)

        return df3.head(10)
   
    def _getSeqs(self, mut: pd):
        
        # get native sequence for future comparison
        with open(f".temp/A0A0K8P6T7.fasta", "r") as f:
            native = f.readlines()[1]
            mut_seq = list(native)
        
        # create sequence with top 5 scoring mutations
        muts = mut.head(5)
        to_mutate = dict(zip(muts['position'], muts['mutated_resi']))
        mpnn_seqs = []
        for i in range(len(native)):
            if (i+1) in to_mutate:
                mut_seq[i] = to_mutate[i+1]
        mpnn_seqs.append(''.join(mut_seq))
        
        # iteratively add next five mutations
        iterative_muts = mut.tail(5)
        to_mutate = dict(zip(iterative_muts['position'], iterative_muts['mutated_resi']))
        for key in to_mutate:
            mut_seq[key-1] = to_mutate[key]
            mpnn_seqs.append(''.join(mut_seq))
            
        return mpnn_seqs


x = mpnnDesignTool("./example_run/A0A0K8P6T7.pdb")
x.run()
    
