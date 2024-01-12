import pandas as pd
import numpy as np
import requests
import cx_Oracle
import re
import ast

def fetch_fasta_sequence(uniprot_accession,begin,end):
    fasta_url = "https://rest.uniprot.org/uniprotkb/"+uniprot_accession+".fasta"
    # Send a GET request to the URL and fetch the content
    response = requests.get(fasta_url)
    if response.status_code == 200:
        fasta_content = response.text
        # Split the content into lines
        lines = fasta_content.split("\n")
        # Initialize variables to store the header and sequence
        header = ""
        sequence = ""
        for line in lines:
            if line.startswith(">"):
                # This is the header line
                header = line[0:]+"\n"
            else:
                # This is a sequence line, concatenate it to the sequence variable
                sequence += line
        # Extract the sequence from residue 9 to 485
        print(begin,end)
        
        sequence = sequence[begin:(end+1)]
        return header+sequence
    
def fetch_aligned_regions(uniprot_accession, pdb_code):
    # Define the GraphQL query
    query = """
    {
        polymer_entities(entity_ids: ["%s_1"]) {
            rcsb_id
            rcsb_polymer_entity_container_identifiers {
                auth_asym_ids
            }
            rcsb_polymer_entity_align {
                reference_database_name
                reference_database_accession
                aligned_regions {
                    entity_beg_seq_id
                    length
                    ref_beg_seq_id
                }
            }
            polymer_entity_instances{
                rcsb_polymer_entity_instance_container_identifiers {
                    asym_id
                    entity_id
                    auth_asym_id
                }
                rcsb_polymer_instance_annotation {
                    type
                    name
                    provenance_source
                    annotation_lineage {
                        id
                        name
                        depth
                    }
                }
            }
        }
    }
    """ % pdb_code

    # Define the API endpoint
    api_endpoint = "https://data.rcsb.org/graphql"

    # Prepare the request
    data = {"query": query}
    response = requests.post(api_endpoint, json=data)

    if response.status_code == 200:
        data = response.json()
        # Find the polymer entity for the specified UniProt accession
        for polymer_entity in data["data"]["polymer_entities"]:
            chain=polymer_entity["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"][0]
            min_begin=[]
            max_end=[]
            ref_begin_seq_id=[]
            for align_entity in polymer_entity["rcsb_polymer_entity_align"]:
                if uniprot_accession == align_entity["reference_database_accession"]:                    
                    aligned_regions = align_entity["aligned_regions"]
                    print(aligned_regions)
                    for region in aligned_regions:
                        min_begin.append(region["entity_beg_seq_id"])
                        ref_begin_seq_id.append(region["ref_beg_seq_id"])
                        max_end.append(region["ref_beg_seq_id"]+region["length"])
            scop2list=[]
            for polymer_entity_instance in polymer_entity["polymer_entity_instances"]:
                if polymer_entity_instance["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"]==chain:
                    if(polymer_entity_instance["rcsb_polymer_instance_annotation"]!=None):
                        for annotation in polymer_entity_instance["rcsb_polymer_instance_annotation"]:
                            if (annotation["type"]=="SCOP2") and (len(annotation["annotation_lineage"])==4) and scop2list==[]:
                                scop2list=[scopclass["id"] for scopclass in annotation["annotation_lineage"]]
                            elif(annotation["type"]=="SCOP2") and (len(annotation["annotation_lineage"])==5):
                                scop2list=[scopclass["id"] for scopclass in annotation["annotation_lineage"]]
            return (chain, np.min(ref_begin_seq_id), np.max(max_end),scop2list)
        
    else:
        print(f"Failed to retrieve data from RCSB GraphQL API. Status code: {response.status_code}")
    return(None)

#Replace 'your_file_path' with the actual path to your file
file_path = 'data_generic.csv'
# Read the data into a Pandas DataFrame
df = pd.read_csv(file_path, delimiter=';')

sequencefile=open("sequences.fa","w")
residuelistfile=open("residueList.tsv","w")
residuelistfile.write("Accession;Residues;SCOP2\n")

for index, row in df.iterrows():
    accession = row['Accession']
    pdb_code = row['PDBCode']
    pdb_str_res_num = ast.literal_eval(row['ResidueList'])
    #1S3B,
    # if accession=='P07900':
    print(accession)
    seq_range=fetch_aligned_regions(accession, pdb_code)
    scop2_list=seq_range[3]
    if(seq_range!=None):
        fasta=fetch_fasta_sequence(accession,seq_range[1],seq_range[2])
        if(len(pdb_str_res_num)>5):
            positions=[int(res[0]) for res in pdb_str_res_num]
            diffs=np.diff(positions)
            regexp=r""
            for ridx,res in enumerate(pdb_str_res_num):
                if(ridx>0):
                    curdiff=diffs[ridx-1]
                    if(curdiff>1):
                        regexp+='.'*(curdiff-1)+res[1]
                    else:
                        regexp+=res[1]
                else:
                    regexp=res[1]            
                
            sequence=fasta.split("\n")[1:][0]
            matches = [match.start() for match in re.finditer(regexp, sequence)]
            # Print the positions of the pattern in the input string
            if matches:
                print(regexp)
                print("Pattern found at positions:", matches)
                residuelist=[matches[0]+1]
                for diff in diffs:
                    residuelist.append(residuelist[-1]+diff)
                residuelistfile.write(accession+";"+str(residuelist)+";"+str(scop2_list)+"\n")
                sequencefile.write(fasta+"\n")
        else:
            print("check pdb code for cavity residue mappings: "+pdb_code)
