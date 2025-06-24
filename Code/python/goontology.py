import requests, sys, json

file_path = 'Data/mart_Cre_Uniprot.txt'

with open(file_path, 'r') as file:
    lines = file.readlines()

uniprot_ids = [line.split('\t')[1] for line in lines if line.strip()]

for uniID in uniprot_ids[1:]:
    uniID = uniID.strip()
    requestURL = f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?includeFields=goName&selectedFields=goId&selectedFields=goAspect&geneProductId={uniID}"

    print(f"Searching for: {uniID}")
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    responseBody = r.text
    data = json.loads(responseBody)
    results = data.get('results', [])

    if results:
        print(f"GO terms found for {uniID}")
        # extract the information from the results
        go_ids = []
        go_names = []
        go_aspects = []
    
        for result in results:
            go_id = result.get('goId', 'N/A')
            go_ids.append(go_id)

            go_name = result.get('goName', 'N/A')
            go_names.append(go_name)

            go_aspect = result.get('goAspect', 'N/A')
            go_aspects.append(go_aspect)
    
        # Save the results into a table
        with open('Data/GO_table.txt', 'a') as f:
            if f.tell() == 0:  # Check if the file is empty
                f.write("UniProtID\tGO_ID\tGO_Name\tGO_Aspect\n")
            for go_id, go_name, go_aspect in zip(go_ids, go_names, go_aspects):
                f.write(f"{uniID}\t{go_id}\t{go_name}\t{go_aspect}\n")


# Read the GO_table.txt file
go_table_path = 'Data/GO_table.txt'

with open(go_table_path, 'r') as go_file:
    table = go_file.readlines()

# Find UniProtIDs with GO_Name containing 'activity'
matching_uniprot_ids = []

for line in table[1:]:  # Skip header
    columns = line.strip().split('\t')
    uniprot_id = columns[0]
    go_id = columns[1]
    go_name = columns[2]
    go_aspect = columns[3]
        
    if 'activity' in go_name.lower():
        matching_uniprot_ids.append(uniprot_id)

        # Save the results into a table
        with open('Data/GO_table_filtered.txt', 'a') as f:
            if f.tell() == 0:  # Check if the file is empty
                f.write("UniProtID\tGO_ID\tGO_Name\tGO_Aspect\n")
            f.write(f"{uniprot_id}\t{go_id}\t{go_name}\t{go_aspect}\n")

len(matching_uniprot_ids)