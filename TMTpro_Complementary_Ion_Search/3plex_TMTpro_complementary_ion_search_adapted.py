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
#   - Raised ION_INTENSITY_THRESHOLD from 0 to 2500
#   - Added M1_INTENSSITY_THRESHOLD = 250 for isotopic peak filtering
#   - Added EMPTY_CHANNEL_COMP_ION_INDEX = 2: the TMT 126 channel (comp_ion3)
#     is empty in this experimental setup; it is computed and reported but not
#     required to exist and not used as a hard filter or ranking criterion
#   - Added total_intensity_no126 (comp_ion1 + comp_ion2) for PSM ranking,
#     so the empty 126 channel does not distort selection of top PSMs
#   - Changed output format from CSV to XLSX
# =============================================================================
import pandas as pd
import numpy as np
import math
from math import comb, log, exp
from pyteomics import mzml


columns_to_keep = ['Spectrum', 'Spectrum File', 'Peptide', 'Modified Peptide', 'Extended Peptide', 'Prev AA',
                   'Next AA', 'Peptide Length', 'Charge', 'Retention', 'Observed Mass', 'Calibrated Observed Mass',
                   'Observed M/Z', 'Calibrated Observed M/Z', 'Calculated Peptide Mass', 'Calculated M/Z',
                   'Delta Mass', 'SpectralSim', 'RTScore', 'Expectation', 'Hyperscore', 'Nextscore', 'Probability',
                   'Number of Enzymatic Termini', 'Number of Missed Cleavages', 'Protein Start', 'Protein End',
                   'Intensity', 'Assigned Modifications', 'Observed Modifications', 'Purity', 'Is Unique',
                   'Protein', 'Protein ID', 'Entry Name', 'Gene', 'Protein Description', 'Mapped Genes',
                   'Mapped Proteins']

# Set the intensity threshold for individual ions
ION_INTENSITY_THRESHOLD = 2500
M1_INTENSSITY_THRESHOLD = 250

# -----------------------------
# IMPORTANT CHANGE:
# 126 channel is EMPTY and maps to comp_ion3 in this script.
# Python indices are 0-based:
#   comp_ion1 -> index 0
#   comp_ion2 -> index 1
#   comp_ion3 -> index 2  (THIS is 126 and is empty)
#
# We STILL compute and report comp_ion3, but:
#   - we do NOT require it to exist
#   - we do NOT require it to pass intensity thresholds
#   - we do NOT use it for ranking/top3 selection later
# -----------------------------
EMPTY_CHANNEL_COMP_ION_INDEX = 2  # 0-based index for comp_ion3 == 126

# Define the mass difference and natural abundance of isotopes
isotope_abundances = {
    'C13': {'element': 'C', 'mass_diff': 1, 'abundance': 0.0107},
    'N15': {'element': 'N', 'mass_diff': 1, 'abundance': 0.00364},
    'H2':  {'element': 'H', 'mass_diff': 1, 'abundance': 0.00015},
    'O17': {'element': 'O', 'mass_diff': 1, 'abundance': 0.000373},
    'O18': {'element': 'O', 'mass_diff': 2, 'abundance': 0.00205}
}

def calculate_elemental_composition(peptide_sequence):
    """Calculate the elemental composition of the peptide"""
    amino_acid_composition = {
        'A': {'C': 3, 'H': 5, 'N': 1, 'O': 1},
        'R': {'C': 6, 'H': 12, 'N': 4, 'O': 1},
        'N': {'C': 4, 'H': 6, 'N': 2, 'O': 2},
        'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3},
        'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},
        'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3},
        'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2},
        'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1},
        'H': {'C': 6, 'H': 7, 'N': 3, 'O': 1},
        'I': {'C': 6, 'H': 11, 'N': 1, 'O': 1},
        'L': {'C': 6, 'H': 11, 'N': 1, 'O': 1},
        'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1},
        'M': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},
        'F': {'C': 9, 'H': 9, 'N': 1, 'O': 1},
        'P': {'C': 5, 'H': 7, 'N': 1, 'O': 1},
        'S': {'C': 3, 'H': 5, 'N': 1, 'O': 2},
        'T': {'C': 4, 'H': 7, 'N': 1, 'O': 2},
        'W': {'C': 11, 'H': 10, 'N': 2, 'O': 1},
        'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2},
        'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1}
    }

    element_counts = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
    for aa in peptide_sequence:
        if aa in amino_acid_composition:
            for element, count in amino_acid_composition[aa].items():
                element_counts[element] += count
        else:
            print(f"Unknown amino acid: {aa}")
            return None

    element_counts['H'] += 2
    element_counts['O'] += 1
    return element_counts

def log_comb(n, k):
    """Calculate the logarithm of the combination number to avoid numerical issues caused by large factorials"""
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def binomial_probability(n, k, p):
    """Calculate the probability of the binomial distribution using logarithmic operations to improve numerical stability"""
    if k > n or k < 0:
        return 0.0
    try:
        log_coeff = log_comb(n, k)
        log_prob = log_coeff + k * log(p) + (n - k) * log(1 - p)
        prob = exp(log_prob)
        return prob
    except (ValueError, OverflowError):
        return 0.0

def calculate_m4_ratio(element_counts):
    """Use the binomial distribution to calculate the ratio of M+4 isotope peaks relative to the monoisotopic peak,
       considering combinations of 13C, 15N, 2H, 17O, and 18O"""
    # Get the number of atoms and isotope abundances
    nC = element_counts.get('C', 0)
    nN = element_counts.get('N', 0)
    nH = element_counts.get('H', 0)
    nO = element_counts.get('O', 0)

    p_C13 = isotope_abundances['C13']['abundance']
    p_N15 = isotope_abundances['N15']['abundance']
    p_H2 = isotope_abundances['H2']['abundance']
    p_O17 = isotope_abundances['O17']['abundance']
    p_O18 = isotope_abundances['O18']['abundance']

    # Calculate the probability P0 of the monoisotopic peak
    P0 = (1 - p_C13) ** nC * (1 - p_N15) ** nN * (1 - p_H2) ** nH * (1 - p_O17 - p_O18) ** nO

    if P0 == 0:
        return 0.0  # If P0 is 0, return 0 directly to avoid division by zero

    # Initialize total probability
    total_probability = 0.0

    # Iterate through all possible isotope substitution combinations with a total mass difference of 4
    max_total_substitutions = 4
    for c13 in range(0, min(nC, max_total_substitutions) + 1):
        for n15 in range(0, min(nN, max_total_substitutions - c13) + 1):
            for h2 in range(0, min(nH, max_total_substitutions - c13 - n15) + 1):
                for o17 in range(0, min(nO, max_total_substitutions - c13 - n15 - h2) + 1):
                    # Calculate the remaining mass difference
                    remaining_mass_diff = max_total_substitutions - (c13 + n15 + h2 + o17)
                    if remaining_mass_diff >= 0 and remaining_mass_diff % 2 == 0:
                        o18 = remaining_mass_diff // 2  # The mass difference for O18 is 2
                        if o18 <= nO - o17:  # Ensure the combination of O18 and O17 is valid
                            # Calculate the probability of each isotope substitution
                            prob_C13 = binomial_probability(nC, c13, p_C13)
                            prob_N15 = binomial_probability(nN, n15, p_N15)
                            prob_H2 = binomial_probability(nH, h2, p_H2)
                            prob_O17 = binomial_probability(nO, o17, p_O17)
                            prob_O18 = binomial_probability(nO - o17, o18, p_O18 / (1 - p_O17))

                            # Accumulate to the total probability
                            total_probability += prob_C13 * prob_N15 * prob_H2 * prob_O17 * prob_O18

    # Calculate the correction ratio
    ratio = total_probability / P0 if P0 > 0 else 0.0
    return ratio

def calculate_m1_ratio(element_counts):
    # (Kept same style as your original code)

    nC = element_counts.get('C', 0)
    nN = element_counts.get('N', 0)
    nH = element_counts.get('H', 0)
    nO = element_counts.get('O', 0)

    p_C13 = isotope_abundances['C13']['abundance']
    p_N15 = isotope_abundances['N15']['abundance']
    p_H2 = isotope_abundances['H2']['abundance']
    p_O17 = isotope_abundances['O17']['abundance']

    # Calculate the probability P0 of the monoisotopic peak
    P0 = (1 - p_C13) ** nC * (1 - p_N15) ** nN * (1 - p_H2) ** nH * (1 - p_O17) ** nO

    if P0 == 0:
        return 0.0  # If P0 is 0, return 0 directly to avoid division by zero

    # Initialize total probability
    total_probability = 0.0

    # Iterate through all possible isotope substitution combinations with a total mass difference of 1
    max_total_substitutions = 1
    for c13 in range(0, min(nC, max_total_substitutions) + 1):
        for n15 in range(0, min(nN, max_total_substitutions - c13) + 1):
            for h2 in range(0, min(nH, max_total_substitutions - c13 - n15) + 1):
                for o17 in range(0, min(nO, max_total_substitutions - c13 - n15 - h2) + 1):
                    # Calculate the remaining mass difference
                    remaining_mass_diff = max_total_substitutions - (c13 + n15 + h2 + o17)
                    if remaining_mass_diff == 0:
                        # Only single substitutions that exactly match the mass difference
                        # Calculate the probability of each isotope substitution
                        prob_C13 = binomial_probability(nC, c13, p_C13)
                        prob_N15 = binomial_probability(nN, n15, p_N15)
                        prob_H2 = binomial_probability(nH, h2, p_H2)
                        prob_O17 = binomial_probability(nO, o17, p_O17)

                        # Since M+1 cannot include O18, we only consider O17
                        total_probability += prob_C13 * prob_N15 * prob_H2 * prob_O17

    # Calculate the correction ratio
    ratio = total_probability / P0 if P0 > 0 else 0.0
    return ratio

def load_fragpipe_psm(file_path):
    """Read FragPipe's PSM.xlsx file, extract peptide information, and retain specified columns"""
    df = pd.read_excel(file_path)

    peptides = []
    for _, row in df.iterrows():
        peptide_info = {col: row[col] for col in columns_to_keep if col in df.columns}
        peptide_info['Peptide'] = row['Peptide']
        peptide_info['charge'] = int(row['Charge'])
        peptide_info['spectrum_id'] = row['Spectrum']

        if 'Observed M/Z' in df.columns:
            precursor_mz = float(row['Observed M/Z'])
            peptide_info['precursor_mz'] = precursor_mz
        else:
            print("The 'Observed M/Z' column was not found in the PSM.xlsx file. The program will exit.")
            return []

        # Keep original behavior: if charge==1, try charge 2/3 versions
        if peptide_info['charge'] == 1:
            new_charges = [2, 3]
            for new_charge in new_charges:
                new_info = peptide_info.copy()
                new_info['precursor_mz'] = (precursor_mz + (new_charge - 1) * 1.00783) / new_charge
                new_info['charge'] = new_charge
                new_info['conversion_flag'] = True
                peptides.append(new_info)
        else:
            peptide_info['conversion_flag'] = False
            peptides.append(peptide_info)
    return peptides

def extract_complementary_ion_intensity(peptides, dia_file_path, ppm_tolerance=20):
    """Extract the intensity of complementary ions and their M+1 isotope peaks for each PSM in the DIA data
       and perform isotope peak corrections"""
    complementary_ion_intensities = []

    # Read the DIA mzML file
    spectra = mzml.read(dia_file_path)
    spectrum_dict = {}
    for spectrum in spectra:
        if spectrum['ms level'] == 2:
            spectrum_id_str = spectrum['id']
            if 'scan=' in spectrum_id_str:
                scan_num = int(spectrum_id_str.split('scan=')[-1])
                spectrum_dict[scan_num] = spectrum

    for peptide in peptides:
        precursor_mz = peptide['precursor_mz']
        charge = peptide['charge']
        spectrum_id = peptide['spectrum_id']
        peptide_sequence = peptide['Peptide']

        try:
            spectrum_id_int = int(spectrum_id.split('.')[-3])
        except ValueError:
            continue

        # Original complementary peak calculation (same as original)
        complementary_ions = [
            (precursor_mz * charge - 163.14648) / (charge - 1),  # Ion 1, lowest comp m/z (135N)
            (precursor_mz * charge - 159.13308) / (charge - 1),  # Ion 2 (131N)
            (precursor_mz * charge - 155.12599) / (charge - 1)   # Ion 3, highest comp m/z (126)
        ]

        # Isotopic complementary peak calculation (same structure)
        iso_mass_diff = 1.00335  # Mass difference of 13C
        isotopic_complementary_ions = [
            complementary_ions[0] + iso_mass_diff / (charge - 1),
            complementary_ions[1] + iso_mass_diff / (charge - 1),
            complementary_ions[2] + iso_mass_diff / (charge - 1)
        ]

        if spectrum_id_int in spectrum_dict:
            spectrum = spectrum_dict[spectrum_id_int]
            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']

            all_peaks_exist = True
            ion_data = []

            # Calculate elemental composition
            element_counts = calculate_elemental_composition(peptide_sequence)
            if element_counts is None:
                continue  # Skip peptides with unknown amino acids

            # Calculate isotope peak ratios (same as original)
            m_plus_4_ratio = calculate_m4_ratio(element_counts)
            m_plus_1_ratio = calculate_m1_ratio(element_counts)

            # Initialize adjusted intensities
            adjusted_intensities = [0.0] * 3

            # Initialize M+1 intensities
            m1_intensities = [0.0] * 3

            # Process ions in order from low mass to high mass (same structure)
            for i in range(3):
                ion_mz = complementary_ions[i]
                tolerance_da = ion_mz * ppm_tolerance / 1e6
                idx = np.where(np.abs(mz_array - ion_mz) <= tolerance_da)[0]

                if len(idx) > 0:
                    max_intensity_idx = idx[np.argmax(intensity_array[idx])]
                    observed_intensity = intensity_array[max_intensity_idx]

                    if i == 0:
                        # The lowest mass ion does not require correction
                        adjusted_intensity = observed_intensity
                    else:
                        # Subtract the M+4 peak contribution of lower mass ions from the current ion
                        M_plus_4_contribution = 0
                        for j in range(i):
                            power = i - j
                            contribution = adjusted_intensities[j] * (m_plus_4_ratio ** power)
                            M_plus_4_contribution += contribution
                        adjusted_intensity = observed_intensity - M_plus_4_contribution
                        adjusted_intensity = max(0.0, adjusted_intensity)  # Ensure non-negative

                    # -----------------------------
                    # IMPORTANT CHANGE :
                    # If this is comp_ion3 (126), do NOT enforce the intensity threshold as a hard criterion.
                    # -----------------------------
                    if i != EMPTY_CHANNEL_COMP_ION_INDEX:
                        if adjusted_intensity < ION_INTENSITY_THRESHOLD:
                            all_peaks_exist = False
                            break

                    adjusted_intensities[i] = adjusted_intensity

                else:
                    # -----------------------------
                    # IMPORTANT CHANGE:
                    # If this is comp_ion3 (126), allow it to be missing (set to 0) and keep going.
                    # Otherwise, keep original behavior (fail).
                    # -----------------------------
                    if i == EMPTY_CHANNEL_COMP_ION_INDEX:
                        adjusted_intensities[i] = 0.0
                    else:
                        all_peaks_exist = False
                        break

            # Extract M+1 isotope intensities for complementary ions (same structure)
            for i in range(3):
                m1_ion_mz = isotopic_complementary_ions[i]
                tolerance_da = m1_ion_mz * ppm_tolerance / 1e6
                idx_m1 = np.where(np.abs(mz_array - m1_ion_mz) <= tolerance_da)[0]
                if len(idx_m1) > 0:
                    m1_max_intensity = np.max(intensity_array[idx_m1])
                    m1_intensities[i] = m1_max_intensity
                else:
                    m1_intensities[i] = 0.0  # If M+1 peak not found, set intensity to 0

            # Check isotopic complementary peaks exist and have sufficient intensity
            # -----------------------------
            # IMPORTANT CHANGE (your request):
            # Ignore comp_ion3 (126) here as a hard criterion.
            # -----------------------------
            for i, ion_mz in enumerate(isotopic_complementary_ions):
                if i == EMPTY_CHANNEL_COMP_ION_INDEX:
                    continue  # do NOT require the 126 isotopic complementary peak

                tolerance_da = ion_mz * ppm_tolerance / 1e6
                idx = np.where(np.abs(mz_array - ion_mz) <= tolerance_da)[0]
                if len(idx) == 0 or np.max(intensity_array[idx]) <= M1_INTENSSITY_THRESHOLD:
                    all_peaks_exist = False
                    break

            # Add outputs if criteria are met
            if all_peaks_exist:
                # -----------------------------
                # IMPORTANT CHANGE :
                # Add a "no126" total intensity to use for ranking/top3 PSMs later
                # (ranking = comp_ion1 + comp_ion2 only).
                # Keep the original total_intensity = sum of all 3 for reference.
                # -----------------------------
                total_intensity_no126 = adjusted_intensities[0] + adjusted_intensities[1]

                peptide.update({
                    'comp_ion1_intensity': adjusted_intensities[0],
                    'comp_ion2_intensity': adjusted_intensities[1],
                    'comp_ion3_intensity': adjusted_intensities[2],  # 126 is still computed/reported
                    'comp_ion1_m1_intensity': m1_intensities[0],
                    'comp_ion2_m1_intensity': m1_intensities[1],
                    'comp_ion3_m1_intensity': m1_intensities[2],
                    'total_intensity': sum(adjusted_intensities),       # all 3 ions
                    'total_intensity_no126': total_intensity_no126,     # ions 1+2 only (for ranking)
                    'm_plus_4_ratio': m_plus_4_ratio,
                    'm_plus_1_ratio': m_plus_1_ratio
                })
                complementary_ion_intensities.append(peptide)

    df = pd.DataFrame(complementary_ion_intensities)
    if not df.empty:
        print(f"Total number of PSMs found: {len(df)}")
    else:
        print("No PSMs met the criteria.")
    return df

def peptide_level_summary(psm_data):
    """Calculate the average intensities of the three complementary ions and average m+1 and m+4 ratios for each peptide"""
    peptide_summary = []

    # Group by Peptide
    grouped = psm_data.groupby('Peptide')

    for peptide, group in grouped:
        # -----------------------------
        # IMPORTANT CHANGE:
        # Sort by total_intensity_no126 (comp_ion1 + comp_ion2) instead of total_intensity.
        # This ensures the empty 126 channel (comp_ion3) is NOT used for ranking.
        # -----------------------------
        if 'total_intensity_no126' in group.columns:
            group = group.sort_values(by='total_intensity_no126', ascending=False)
        else:
            # Fallback (should not happen if extract_complementary_ion_intensity ran)
            group = group.sort_values(by='total_intensity', ascending=False)

        # Select the top three PSMs with the highest intensity (by the chosen ranking metric)
        top_psms = group.head(3)

        # Calculate the averages
        avg_comp_ion1 = top_psms['comp_ion1_intensity'].mean()
        avg_comp_ion2 = top_psms['comp_ion2_intensity'].mean()
        avg_comp_ion3 = top_psms['comp_ion3_intensity'].mean()
        avg_m_plus_4_ratio = top_psms['m_plus_4_ratio'].mean()
        avg_m_plus_1_ratio = top_psms['m_plus_1_ratio'].mean()

        # Calculate average M+1 intensities
        avg_comp_ion1_m1 = top_psms['comp_ion1_m1_intensity'].mean()
        avg_comp_ion2_m1 = top_psms['comp_ion2_m1_intensity'].mean()
        avg_comp_ion3_m1 = top_psms['comp_ion3_m1_intensity'].mean()

        peptide_summary.append({
            'Peptide': peptide,
            'Average Comp Ion 1 Intensity': avg_comp_ion1,
            'Average Comp Ion 2 Intensity': avg_comp_ion2,
            'Average Comp Ion 3 Intensity': avg_comp_ion3,  # 126 is still reported
            'Average m+4 Ratio': avg_m_plus_4_ratio,
            'Average m+1 Ratio': avg_m_plus_1_ratio,
            'Average Comp Ion 1 M+1 Intensity': avg_comp_ion1_m1,
            'Average Comp Ion 2 M+1 Intensity': avg_comp_ion2_m1,
            'Average Comp Ion 3 M+1 Intensity': avg_comp_ion3_m1
        })

    return pd.DataFrame(peptide_summary)

# Main program
if __name__ == "__main__":
    fragpipe_psm_path = 'HelaYeast_DIA1.xlsx'  # Please replace with your actual file path
    dia_mzml_path = 'HelaYeast_DIA1.mzML'      # Please replace with your actual file path

    # Read PSM information
    peptides = load_fragpipe_psm(fragpipe_psm_path)

    # Extract complementary ion intensities
    psm_intensity_data = extract_complementary_ion_intensity(peptides, dia_mzml_path, ppm_tolerance=20)

    # Export PSM-level results as XLSX (your request)
    psm_intensity_data.to_excel('DIA1_psm_int.xlsx', index=False)

    # Generate Peptide-level statistics
    peptide_summary = peptide_level_summary(psm_intensity_data)

    # Export Peptide-level results as XLSX (your request)
    peptide_summary.to_excel('DIA1_pept_summary.xlsx', index=False)