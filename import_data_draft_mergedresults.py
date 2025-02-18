import psycopg2
import pandas as pd
import numpy as np  # Required for handling infinity values

# PostgreSQL connection parameters
db_params = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '115.146.85.41',
    'port': 5432
}

# File containing draft genome data
draft_genomes_path = "merged_results.tsv"  # Update with actual file path

# Import draft genome data
print(f"Importing data from {draft_genomes_path}")

# Load data
draft_genomes = pd.read_csv(draft_genomes_path, sep="\t")

# Normalize column names (remove spaces & make lowercase)
draft_genomes.columns = draft_genomes.columns.str.strip().str.lower()

# Print column names to verify
print("✅ Columns in file:", draft_genomes.columns.tolist())

# Function to remove % and convert to float
def clean_percentage(column):
    return column.astype(str).str.replace("%", "", regex=True).astype(float)

# Clean percentage columns Apply to all percentage-containing columns
percentage_columns = ["homozygosity", "heterozygosity", "modelfit", "errorrate", "percent_gaps"]  # Add more if needed

for col in percentage_columns:
    if col in draft_genomes.columns:
        print(f"🔹 Cleaning percentage column: {col}")
        draft_genomes[col] = clean_percentage(draft_genomes[col])

# Replace 'NA', 'N/A', and other non-standard missing values with NaN
draft_genomes.replace(["NA", "N/A", "NULL", "missing"], np.nan, inplace=True)

# Convert remaining object (string) columns to numeric if possible
for col in draft_genomes.columns:
    if draft_genomes[col].dtype == "object":
        draft_genomes[col] = pd.to_numeric(draft_genomes[col], errors="coerce")

# Replace remaining NaN with None for PostgreSQL
draft_genomes = draft_genomes.astype(object).where(pd.notna(draft_genomes), None)

# Define integer-like columns
int_columns = [
    "passed_filter_reads", "low_quality_reads", "too_many_n_reads", "too_short_reads", "too_long_reads",
    "raw_total_reads", "raw_total_bases", "raw_q20_bases", "raw_q30_bases",
    "raw_read1_mean_length", "raw_read2_mean_length", "total_reads", "total_bases", "q20_bases", "q30_bases",
    "read1_mean_length", "read2_mean_length", "genomesize", "repeatsize", "uniquesize", "num_contigs",
    "num_contigs_mitochondrion", "num_contigs_plastid", "num_contigs_prokarya",
    "bp_mitochondrion", "bp_plastid", "bp_prokarya", "n_markers",
    "number_of_scaffolds", "number_of_contigs", "total_length", "scaffold_n50", "contigs_n50",
    "unique_k_mers_assembly", "k_mers_total", "solid_k_mers", "total_k_mers", "readbp", "estgenomesize"
]

# Print missing values BEFORE conversion
print("\n🔍 Missing values before conversion:")
print(draft_genomes.replace([np.inf, -np.inf], np.nan).isna().sum())


# Define decimal-like columns
decimal_columns = [
    "raw_q20_rate", "raw_q30_rate", "raw_gc_content", "q20_rate", "q30_rate", "gc_content",
    "homozygosity", "heterozygosity", "modelfit", "errorrate",
    "complete", "single_copy", "multi_copy", "fragmented", "missing",
    "percent_gaps", "qv", "error", "completeness", "mapadjust", "scdepth"
]

# Convert to numeric before checking for infinite values
for col in int_columns + decimal_columns:
    if col in draft_genomes.columns:
        print(f"🔹 Processing column: {col}")

        # Convert column to numeric first to avoid type errors
        draft_genomes[col] = pd.to_numeric(draft_genomes[col], errors="coerce")

        # Detect and print inf values before conversion
        if draft_genomes[col].dtype in ["float64", "int64"]:  # Only apply np.isinf() to numeric columns
            if np.isinf(draft_genomes[col]).sum() > 0:
                print(f"⚠️ WARNING: Column '{col}' contains infinite values. Replacing with NaN.")
                draft_genomes[col] = draft_genomes[col].replace([np.inf, -np.inf], np.nan)

        # Convert to Int64 safely if it's an integer column
        if col in int_columns:
            draft_genomes[col] = draft_genomes[col].astype("Int64")



# Convert decimal columns safely
for col in decimal_columns:
    if col in draft_genomes.columns:
        # Detect and print inf values before conversion
        if np.isinf(draft_genomes[col]).sum() > 0:
            print(f"⚠️ WARNING: Column '{col}' contains infinite values. Replacing with NaN.")

        draft_genomes[col] = draft_genomes[col].replace([np.inf, -np.inf], np.nan)
        draft_genomes[col] = pd.to_numeric(draft_genomes[col], errors="coerce")

# Replace NaN with None for database compatibility
draft_genomes = draft_genomes.astype(object).where(pd.notna(draft_genomes), None)

# Print summary of changes
print("\n🔍 Final dataset summary:")
print(draft_genomes.describe())

try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    row_count = 0  # Track number of processed rows

    for index, row in draft_genomes.iterrows():
        values = tuple(row.values)

        # Extract primary key values
        sample, seq_date = row["sample"], row["seq_date"]

        # UPSERT: Insert if not exists, otherwise update
        upsert_query = """
        INSERT INTO draft_genomes (
            sample, seq_date, mach, initial, passed_filter_reads, low_quality_reads, too_many_n_reads, too_short_reads, too_long_reads,
            raw_total_reads, raw_total_bases, raw_q20_bases, raw_q30_bases, raw_q20_rate, raw_q30_rate, raw_read1_mean_length, raw_read2_mean_length,
            raw_gc_content, total_reads, total_bases, q20_bases, q30_bases, q20_rate, q30_rate, read1_mean_length, read2_mean_length, gc_content,
            homozygosity, heterozygosity, genomesize, repeatsize, uniquesize, modelfit, errorrate, num_contigs, num_contigs_mitochondrion,
            num_contigs_plastid, num_contigs_prokarya, bp_mitochondrion, bp_plastid, bp_prokarya, complete, single_copy, multi_copy, fragmented, missing,
            n_markers, domain, number_of_scaffolds, number_of_contigs, total_length, percent_gaps, scaffold_n50, contigs_n50, unique_k_mers_assembly,
            k_mers_total, qv, error, k_mer_set, solid_k_mers, total_k_mers, completeness, depmethod, adjust, readbp, mapadjust, scdepth, estgenomesize
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (sample, seq_date) DO UPDATE SET
            completeness = EXCLUDED.completeness;
        """

        cursor.execute(upsert_query, values)
        row_count += 1  

    conn.commit()
    print(f"✅ Successfully processed {row_count} rows!")

except Exception as e:
    conn.rollback()
    print(f"❌ Error: {e}")

finally:
    cursor.close()
    conn.close()
