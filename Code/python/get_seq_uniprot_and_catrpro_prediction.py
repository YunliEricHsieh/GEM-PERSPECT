import requests
import io
import pandas as pd

# Read the CSV file
csv_file = 'Results/CataPro/enzymes_and_substrates_for_reference.csv'
df = pd.read_csv(csv_file)

# Extract the uniprot IDs from the appropriate column
uniprot_IDs = df['Enzyme'].dropna().tolist()

for uniID in uniprot_IDs:
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniID}')
    seq = response.json()['sequence']['value']
    df.loc[df['Enzyme'] == uniID, 'Sequence'] = seq

# Save the updated DataFrame back to the CSV file
df.to_csv(csv_file, index=False)

## Predict kcat values for a given dataset using the CataPro model.
python CataPro/inference/predict.py \
            -inp_fpath Results/CataPro/enzymes_and_substrates_for_reference.csv \
            -model_dpath CataPro/models \
            -batch_size 128 \
            -device cpu \
            -out_fpath Results/CataPro/enzymes_and_substrates_for_reference_kcats.csv    