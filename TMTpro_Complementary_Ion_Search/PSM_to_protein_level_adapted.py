# =============================================================================
# Adapted from:
#   Title:    TMTproC-DIA Complementary-Ion Search Python Scripts
#   Authors:  Zicong Wang, Peng-Kai Liu, Haiyan Lu, Lingjun Li
#             (University of Wisconsin–Madison)
#   Source:   https://zenodo.org/records/14567079
#             DOI: 10.5281/zenodo.14567079
#   Paper:    Wang Z, Liu P-K, Lu H, Li L.
#             "Multiplexed Data-Independent Acquisition-Based Proteomics
#             Enabled by TMTpro Complementary Ions."
#             Analytical Chemistry. 2025;97(48):26419.
#             DOI: 10.1021/acs.analchem.5c03563
#   License:  MIT (see LICENSE file)
#
# Modifications by Tim Thierer, IMP Vienna, 2026:
#   - Simplified load_data() to Excel-only (pd.read_excel) to match pipeline
#   - Added ranking_columns parameter to aggregate_psm_to_peptide() and
#     aggregate_peptide_to_protein(): PSM/peptide ranking excludes comp_ion3
#     (TMT 126 channel, which is empty in this experimental setup)
#   - Hyperscore is no longer a mandatory required_column; normalization
#     factor is determined flexibly via determine_normalization_factor()
#   - Changed output format from CSV to XLSX
#   - Updated file paths to match output filenames from script 1
# =============================================================================
import pandas as pd
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_data(file_path):
    try:
        data = pd.read_excel(file_path)
        logging.info(f"Successfully loaded data from {file_path}")
        return data
    except FileNotFoundError:
        logging.error(f"Error: File not found {file_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading file: {e}")
        sys.exit(1)


def check_columns(data, required_columns):
    # Remove extra spaces from column names
    data.columns = [col.strip() for col in data.columns]
    logging.info(f"Column names in data: {list(data.columns)}")

    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        logging.error(f"Error: Missing columns in input file: {missing_columns}")
        sys.exit(1)
    else:
        logging.info("All required columns are present.")


def determine_normalization_factor(data):
    if 'Hyperscore' in data.columns:
        normalization_factor = 'Hyperscore'
    elif 'Probability' in data.columns:
        normalization_factor = 'Probability'
    else:
        logging.warning("Neither 'Hyperscore' nor 'Probability' column found. Creating default 'Hyperscore' column and setting it to 1.0")
        data['Hyperscore'] = 1.0
        normalization_factor = 'Hyperscore'
    logging.info(f"Using {normalization_factor} as the normalization factor")
    return normalization_factor


def preprocess_data(data, reporter_columns):
    # Ensure reporter intensity columns are numeric and handle missing values
    for col in reporter_columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    data[reporter_columns] = data[reporter_columns].fillna(0)
    logging.info("Reporter intensity columns converted to numeric and missing values handled.")

    # Print the first 10 values of the 'Is Unique' column to check original data
    logging.info("First 10 values of the original 'Is Unique' column:")
    logging.info(data['Is Unique'].head(10))

    # Convert 'Is Unique' column to boolean type
    data['Is Unique'] = data['Is Unique'].astype(str).str.strip().str.upper().replace({
        '1': True,
        'YES': True,
        'Y': True,
        'TRUE': True,
        '0': False,
        'NO': False,
        'N': False,
        'FALSE': False
    })
    data['Is Unique'] = data['Is Unique'].astype(bool)
    logging.info("'Is Unique' column converted to boolean type.")

    # Verify the converted 'Is Unique' column
    unique_values = data['Is Unique'].unique()
    value_counts = data['Is Unique'].value_counts()
    logging.info(f"Unique values in the converted 'Is Unique' column: {unique_values}")
    logging.info(f"Value counts in the converted 'Is Unique' column:\n{value_counts}")

    # Print the first 10 rows of data to verify conversion
    logging.info("First 10 rows of data after conversion:")
    logging.info(data[['Is Unique', 'Protein', 'Peptide']].head(10))

    # Handle missing 'Protein' column entries
    protein_missing = data['Protein'].isnull().sum()
    if protein_missing > 0:
        logging.warning(f"There are {protein_missing} rows with missing 'Protein' column. These rows will be removed.")
        data = data.dropna(subset=['Protein'])
    else:
        logging.info("No missing values in the 'Protein' column.")

    # Check for missing values in 'Peptide' and 'Protein' columns
    peptide_missing = data['Peptide'].isnull().sum()
    protein_missing = data['Protein'].isnull().sum()
    logging.info(f"Number of missing values in 'Peptide' column: {peptide_missing}")
    logging.info(f"Number of missing values in 'Protein' column: {protein_missing}")

    return data


def aggregate_psm_to_peptide(data, reporter_columns, ranking_columns, normalization_factor, top_n=3):
    logging.info("Starting aggregation from PSM to Peptide level")

    # -----------------------------
    # IMPORTANT CHANGE:
    # 126 channel is comp_ion3_intensity and is empty -> exclude it from ranking only.
    # We still compute/report comp_ion3, but ranking is based on comp_ion1+comp_ion2.
    # -----------------------------
    data['ranking_intensity'] = data[ranking_columns].sum(axis=1)

    # Sort PSMs by ranking intensity (NOT by the sum including comp_ion3)
    data_sorted = data.sort_values(by=['Peptide', 'ranking_intensity'], ascending=[True, False])

    # Select the top N PSMs for each peptide
    top_psms = data_sorted.groupby('Peptide').head(top_n)
    logging.info(f"Selected the top {top_n} PSMs for each peptide (ranking excludes 126/comp_ion3)")

    # Aggregate reporter intensities and normalization factor
    aggregation_columns = reporter_columns + [normalization_factor]
    peptide_level = top_psms.groupby('Peptide')[aggregation_columns].mean().reset_index()

    logging.info(f"Number of aggregated peptides: {len(peptide_level)}")

    if normalization_factor not in peptide_level.columns:
        logging.error(f"Aggregation resulted in missing '{normalization_factor}' column.")
        sys.exit(1)

    # Add Protein information
    peptide_level = pd.merge(peptide_level, data[['Peptide', 'Protein']].drop_duplicates(), on='Peptide')
    logging.info("Completed aggregation from PSM to Peptide level")
    logging.info(f"Preview of aggregated peptide data:\n{peptide_level.head()}")

    return peptide_level


def aggregate_peptide_to_protein(data, reporter_columns, ranking_columns, normalization_factor, top_n_peptides=3):
    data = data.copy()
    logging.info("Starting aggregation from Peptide to Protein level")

    # -----------------------------
    # IMPORTANT CHANGE:
    # Rank peptides per protein using ONLY comp_ion1+comp_ion2 (exclude comp_ion3 / 126).
    # -----------------------------
    data['sum_intensity'] = data[ranking_columns].sum(axis=1)

    # Sort peptides based on the ranking sum (excluding 126)
    data_sorted = data.sort_values(by=['Protein', 'sum_intensity'], ascending=[True, False])

    # Select the top N peptides for each protein
    top_peptides = data_sorted.groupby('Protein').head(top_n_peptides)
    logging.info(f"Selected the top {top_n_peptides} peptides for each protein (ranking excludes 126/comp_ion3)")

    # Calculate weighted intensities (keep ALL reporter columns in the final quant table)
    for col in reporter_columns:
        top_peptides[f'{col}_weighted'] = top_peptides[col] * top_peptides[normalization_factor]
        logging.debug(f"Calculated weighted intensity: {col}_weighted")

    # Aggregate weighted intensities and normalization factor
    agg_dict = {f'{col}_weighted': 'sum' for col in reporter_columns}
    agg_dict[normalization_factor] = 'sum'

    protein_level = top_peptides.groupby('Protein').agg(agg_dict).reset_index()
    logging.info(f"Number of aggregated proteins: {len(protein_level)}")

    # Check if normalization factor is zero
    zero_norm = protein_level[normalization_factor] == 0
    if zero_norm.any():
        logging.warning(f"Proteins with a total normalization factor of zero found in '{normalization_factor}'. Replacing with 1 to avoid division by zero.")
        protein_level.loc[zero_norm, normalization_factor] = 1

    # Calculate normalized intensities
    for col in reporter_columns:
        protein_level[col] = protein_level[f'{col}_weighted'] / protein_level[normalization_factor]
        logging.debug(f"Calculated normalized intensity: {col}")

    # Select the final protein-level data
    protein_level = protein_level[['Protein'] + reporter_columns]
    logging.info("Completed aggregation from Peptide to Protein level")
    logging.info(f"Preview of aggregated protein data:\n{protein_level.head()}")

    return protein_level


def save_data(data, filename):
    # -----------------------------
    # IMPORTANT CHANGE (your request):
    # Save as XLSX instead of CSV.
    # -----------------------------
    try:
        data.to_excel(filename, index=False)
        logging.info(f"Data saved to: {filename}")
    except Exception as e:
        logging.error(f"Error saving file {filename}: {e}")


def main():
    logging.info("Script execution started")

    # Parameter settings
    file_path = 'DIA1_psm_int.xlsx'  # Change to your file path

    # Reporter columns from the first script (comp_ion3 == 126)
    reporter_columns = ['comp_ion1_intensity', 'comp_ion2_intensity', 'comp_ion3_intensity']

    # Ranking columns (exclude comp_ion3 / 126)
    ranking_columns = ['comp_ion1_intensity', 'comp_ion2_intensity']

    # -----------------------------
    # IMPORTANT CHANGE:
    # Do NOT require Hyperscore as a mandatory column here,
    # because the script can use Probability or create a default Hyperscore.
    # -----------------------------
    required_columns = ['Is Unique', 'Protein', 'Peptide',
                        'comp_ion1_intensity', 'comp_ion2_intensity', 'comp_ion3_intensity']

    top_n_psm = 3
    top_n_peptides = 3

    # Load data
    data = load_data(file_path)

    # Check for required columns
    check_columns(data, required_columns)

    # Determine normalization factor
    normalization_factor = determine_normalization_factor(data)

    # Preprocess data
    data = preprocess_data(data, reporter_columns)

    # Filter unique peptides
    data_unique = data[data['Is Unique'] == True].copy()
    logging.info(f"Number of unique peptides after filtering: {len(data_unique)}")

    if data_unique.empty:
        logging.warning("No data with unique peptides after filtering. Please check the values in the 'Is Unique' column.")
        sys.exit(1)

    # Aggregate from PSM to Peptide level (ranking excludes 126)
    peptide_level_data = aggregate_psm_to_peptide(
        data_unique,
        reporter_columns=reporter_columns,
        ranking_columns=ranking_columns,
        normalization_factor=normalization_factor,
        top_n=top_n_psm
    )

    # Aggregate from Peptide to Protein level (ranking excludes 126)
    protein_level_data = aggregate_peptide_to_protein(
        peptide_level_data,
        reporter_columns=reporter_columns,
        ranking_columns=ranking_columns,
        normalization_factor=normalization_factor,
        top_n_peptides=top_n_peptides
    )

    peptide_level_file = 'DIA1_pept_level.xlsx'
    protein_level_file = f'DIA1_prot_level_{normalization_factor.lower()}.xlsx'

    save_data(peptide_level_data, peptide_level_file)
    save_data(protein_level_data, protein_level_file)

    logging.info("Peptide-level data preview:")
    print(peptide_level_data.head())
    logging.info("Protein-level data preview:")
    print(protein_level_data.head())

    logging.info("Script execution completed")


if __name__ == "__main__":
    main()