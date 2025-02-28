import os
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from collections import defaultdict
from datetime import datetime

# ---- CONFIG ----
BASE_FOLDER = './grade_generalised'

# Create unique filename with timestamp
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
OUTPUT_EXCEL = f'ESCC_Genomic_Summary_{timestamp}.xlsx'

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
        'snp_data': []
    }

    for sample_folder in os.listdir(stage_folder):
        sample_path = os.path.join(stage_folder, sample_folder)
        if not os.path.isdir(sample_path):
            continue

        for file in os.listdir(sample_path):
            file_path = os.path.join(sample_path, file)

            if 'Copy_Number_Variation' in file:
                cnv_df = read_cnv_file(file_path)
                if cnv_df is not None:
                    segment_col = get_segment_mean_column(cnv_df)
                    if segment_col:
                        stage_data['cnv_segment_means'].extend(cnv_df[segment_col].tolist())

            elif 'Simple_Nucleotide_Variation' in file or file.endswith('.maf'):
                snp_df = read_maf_file(file_path)
                if snp_df is not None:
                    stage_data['snp_data'].append(snp_df)

    return stage_data

# ---- Create Summary Tables for Each Stage ----
def create_stage_summary(stage, stage_data):
    summary = {}

    # CNV Segment Means Summary
    if stage_data['cnv_segment_means']:
        cnv_stats = pd.Series(stage_data['cnv_segment_means']).describe()
        summary['CNV_Segment_Stats'] = cnv_stats.to_frame(name='Value')

    # SNP Mutation Summary
    if stage_data['snp_data']:
        combined_snp_data = pd.concat(stage_data['snp_data'], ignore_index=True)

        # Top 10 mutated genes
        top_genes = combined_snp_data['Hugo_Symbol'].value_counts().head(10)
        summary['Top_Mutated_Genes'] = top_genes.to_frame(name='Mutation_Count')

        # Mutation Classification Distribution
        mutation_classification = combined_snp_data['Variant_Classification'].value_counts()
        summary['Mutation_Classification'] = mutation_classification.to_frame(name='Count')

    return summary

# ---- Write Summary to Excel ----
def save_summaries_to_excel(stage_summaries):
    with pd.ExcelWriter(OUTPUT_EXCEL) as writer:
        for stage, summaries in stage_summaries.items():
            for sheet_name, df in summaries.items():
                sheet_title = f'{stage}_{sheet_name}'
                
                # Excel has a 31 character limit on sheet names
                if len(sheet_title) > 31:
                    sheet_title = sheet_title[:28] + '...'

                df.to_excel(writer, sheet_name=sheet_title)

    print(f"✅ All stage summaries saved to '{OUTPUT_EXCEL}'")

# ---- Main ----
def main():
    print("🔎 Starting analysis across all stages in:", BASE_FOLDER)

    stage_summaries = {}

    for stage in os.listdir(BASE_FOLDER):
        stage_folder = os.path.join(BASE_FOLDER, stage)
        if not os.path.isdir(stage_folder):
            continue

        print(f"📊 Processing Stage: {stage}")
        stage_data = process_stage(stage, stage_folder)

        stage_summary = create_stage_summary(stage, stage_data)
        stage_summaries[stage] = stage_summary

    save_summaries_to_excel(stage_summaries)

if __name__ == "__main__":
    main()
