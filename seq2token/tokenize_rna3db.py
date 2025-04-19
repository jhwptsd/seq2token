# Import libraries
import os
import json
import urllib.request
import re
from Bio.PDB import PDBParser, PDBIO, MMCIFParser

# Parse RNA3db Sequences file tree
def parse_json(n, max_len=150):

    if not os.path.exists("./rna3db-jsons"):
        print("Downloading RNA3Db...")
        os.system('wget https://github.com/marcellszi/rna3db/releases/download/incremental-update/rna3db-jsons.tar.gz')
        print("Extracting sequence files...")
        os.system('tar -xzf rna3db-jsons.tar.gz')
        os.system('rm rna3db-jsons.tar.gz')
    path = "rna3db-jsons/cluster.json"

    num = -1
    seqs = {}
    f = open(path)
    data = json.load(f)
    for _, j_dict in data.items():
        for _, k_dict in j_dict.items():
            for k, details in k_dict.items():
                num += 1
                if num > n:
                    break
                if details["length"] > max_len:
                    n -= 1
                    continue
                seqs[k] = details["sequence"]
    f.close()
    return seqs

# Retrieve PDB files from RCSB
def download_file(pdb_id):
    urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdb_id}.cif', f'files/{pdb_id}.pdb')

# Download PDB file for each sequence
def download_rna3db(num):
    seqs = parse_json(num)
    os.makedirs("files", exist_ok=True)
    print("Downloading PDB files...")
    count = 0
    for _, (pdb_id, _) in enumerate(seqs.items()):
        _location = pdb_id.find("_")
        id = pdb_id.upper()[:_location]
        if not os._exists(f"./files/{id}.pdb"):
            try:
                download_file(id)
                convert_cif_to_pdb(f'files/{id}.pdb')
                count += 1
            except:
                continue
    print(f"Downloaded {count} PDB files.")
    return seqs

# Convert downloaded CIF files to PDB format
def convert_cif_to_pdb(filename):
    try:
        parser = MMCIFParser()
        structure = parser.get_structure("structure_id", filename)
    except Exception as e:
         parser = PDBParser()
         structure = parser.get_structure("structure_id", filename)

    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)

if __name__=='__main__':
    seqs = download_rna3db(2)

    # Edit test_pdb.yaml file for correct filenames and directories
    rep1 = re.compile('results_dir: results/pdbs/.*')
    rep2 = re.compile('pdb_path: files/.*.pdb')
    for i, (pdb_id, seq) in enumerate(seqs.items()):
        _location = pdb_id.find("_")
        id = pdb_id.upper()[:_location]
        with open('./configs/test_pdb.yaml', 'r+') as f:
            data = rep1.sub(f'results_dir: results/pdbs/{id}', f.read())
            data = rep2.sub(f'pdb_path: files/{id}.pdb', data)
            f.seek(0)
            f.write(data)
            f.truncate()
        
        # Run bio2token for each pdb file
        os.system('uv run scripts/test_pdb.py --config test_pdb.yaml')

        # Uncharted territory!
        # Next step is to figure out the structure of bio2token output files so I can condense it into a csv/dict with PDB ID and sequence