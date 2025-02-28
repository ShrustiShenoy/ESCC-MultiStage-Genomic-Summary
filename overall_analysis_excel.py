import os
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from collections import defaultdict

# ---- CONFIG ----
METADATA_FILE = 'Filtered_clinical_path.xlsx'
BASE_FOLDER = './grade_generalised'
STAGEII_SUBTYPES = ['II', 'IIA', 'IIB']
OUTPUT_EXCEL_FILE = 'ESCC_Overall_Genomic_Summary_2025-02-28.xlsx'
FAILED_FILES_LOG = 'failed_files.log'

# ---- File Readers ----
def read_cnv_file(file_path):
    try:
        cnv_df = pd.read_csv(file_path, sep='\t')
        return cnv_df
    except Exception as e:
        print(f"❌ Failed to read CNV file {file_path}: {e}")
        return None

def read_maf_file(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]
        return pd.read_csv(StringIO("\n".join(lines)), sep='\t')
    except Exception as e:
        print(f"❌ Failed to read MAF file {file_path}: {e}")
        return None

# ---- Helper Function to Find Segment Mean Column ----
def get_segment_mean_column(cnv_df):
    for col in cnv_df.columns:
        if col.lower() in ['segment_mean', 'segmentmean']:
            return col
    return None

# ---- Process One Stage ----
def process_stage(stage, stage_folder):
    stage_data = {
        'cnv_segment_means': [],
        'snp_data': [],
        'cnv_segment_means_by_sample': defaultdict(list),
        'snp_data_by_sample': defaultdict(list),
    }

    for sample_folder in os.listdir(stage_folder):
        sample_path = os.path.join(stage_folder, sample_folder)
        if not os.path.isdir(sample_path):
            continue

        print(f"📂 Processing sample: {sample_folder}")

        for file in os.listdir(sample_path):
            file_path = os.path.join(sample_path, file)

            if 'Copy_Number_Variation' in file:
                cnv_df = read_cnv_file(file_path)
                if cnv_df is not None:
                    segment_col = get_segment_mean_column(cnv_df)
                    if segment_col:
                        segment_means = cnv_df[segment_col].tolist()
                        stage_data['cnv_segment_means'].extend(segment_means)
                        stage_data['cnv_segment_means_by_sample'][sample_folder].extend(segment_means)
                    else:
                        print(f"⚠️ Warning: No Segment Mean column found in {file_path}.")
                        log_failed_file(file_path, "Missing Segment Mean column")

            elif 'Simple_Nucleotide_Variation' in file or file.endswith('.maf'):
                snp_df = read_maf_file(file_path)
                if snp_df is not None:
                    stage_data['snp_data'].append(snp_df)
                    stage_data['snp_data_by_sample'][sample_folder].append(snp_df)

    return stage_data

# ---- Logging ----
def log_failed_file(file_path, reason):
    with open(FAILED_FILES_LOG, 'a') as log:
        log.write(f"{file_path}\t{reason}\n")

# ---- Create Excel Summary ----
def save_summary_to_excel(cnv_segment_means, snp_data_by_sample):
    with pd.ExcelWriter(OUTPUT_EXCEL_FILE) as writer:

        # --- CNV Segment Summary ---
        if cnv_segment_means:
            cnv_stats = pd.Series(cnv_segment_means).describe()
            cnv_stats.to_frame(name='CNV Segment Mean Statistics').to_excel(writer, sheet_name='CNV_Segment_Stats')

        # --- Combined SNP Data ---
        if snp_data_by_sample:
            combined_snp_df = pd.concat([df for dfs in snp_data_by_sample.values() for df in dfs], ignore_index=True)

            # Top 10 mutated genes
            top_genes = combined_snp_df['Hugo_Symbol'].value_counts().head(10)
            top_genes.to_frame(name='Mutation_Count').to_excel(writer, sheet_name='Top_Mutated_Genes')

            # Mutation classification counts
            mutation_class_counts = combined_snp_df['Variant_Classification'].value_counts()
            mutation_class_counts.to_frame(name='Count').to_excel(writer, sheet_name='Mutation_Classification')

    print(f"✅ Summary saved to '{OUTPUT_EXCEL_FILE}'")

# ---- Main ----
def main():
    print("🔎 Starting analysis across all stages in:", BASE_FOLDER)

    all_cnv_segment_means = []
    all_snp_data_by_sample = defaultdict(list)

    # Clear any old failed file log
    if os.path.exists(FAILED_FILES_LOG):
        os.remove(FAILED_FILES_LOG)

    for stage in os.listdir(BASE_FOLDER):
        stage_folder = os.path.join(BASE_FOLDER, stage)

        if not os.path.isdir(stage_folder):
            continue

        print(f"\n🧬 Processing Stage: {stage}")
        stage_data = process_stage(stage, stage_folder)

        all_cnv_segment_means.extend(stage_data['cnv_segment_means'])
        for sample, snp_dfs in stage_data['snp_data_by_sample'].items():
            all_snp_data_by_sample[sample].extend(snp_dfs)

    # Save all summaries to Excel
    save_summary_to_excel(all_cnv_segment_means, all_snp_data_by_sample)

    print(f"\n✅ Analysis complete. Check '{OUTPUT_EXCEL_FILE}' for results.")
    print(f"⚠️ Check '{FAILED_FILES_LOG}' for any problematic files.")

if __name__ == "__main__":
    main()
