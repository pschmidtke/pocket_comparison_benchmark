import pandas  as pd
from chembl_webresource_client.new_client import new_client
import numpy as np

file_path = 'data.csv'
outpath = 'chembl_compounds.csv'
# Read the data into a Pandas DataFrame
df = pd.read_csv(file_path, delimiter=';')

targets_api = new_client.target
bioactivities_api = new_client.activity
compounds_api = new_client.molecule

with open(outpath,'w') as out:
    out.write('accession;compound_id;smiles\n')
    for index, row in df.iterrows():
        accession=row["Accession"]
        print(accession)
        # if accession=='P07900':
        targets=targets_api.get(target_components__accession=accession, type__in=["SINGLE PROTEIN","PROTEIN FAMILY","SELECTIVITY GROUP"]).only("target_chembl_id", "pref_name", "target_type","comment")
        if targets:
            target_ids = [target['target_chembl_id'] for target in targets]
        # print(target_ids)
        # for target_id in target_ids:
        bioactivities = bioactivities_api.filter(target_chembl_id__in=target_ids, type__in=["IC50","Kd","Ki","Potency","XC50","EC50","AC50"],confidence_score__in=[6,7,8,9], relation="=", assay_type="B").only(
            "activity_id",
            "activity_comment"
            "molecule_chembl_id",
            "type",
            "standard_units",
            "relation",
            "standard_value",
            "pchembl_value"
        )
        # print(bioactivities)
        compound_ids=list(np.unique([activity["molecule_chembl_id"] for activity in bioactivities if(activity["pchembl_value"] is not None and float(activity["pchembl_value"])>=5)]))
        # print(compound_ids)
        compounds_provider = compounds_api.filter(molecule_chembl_id__in=compound_ids).only("molecule_chembl_id", "molecule_structures")
        for compound in compounds_provider:
                if compound["molecule_structures"] is not None and compound["molecule_structures"]["canonical_smiles"] is not None:
                    out.write(f'{accession};{compound["molecule_chembl_id"]};{compound["molecule_structures"]["canonical_smiles"]}\n')