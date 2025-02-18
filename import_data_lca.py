import psycopg2
import pandas as pd

# PostgreSQL connection parameters
db_params = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '115.146.85.41',  # Update if necessary
    'port': 5432  # Default PostgreSQL port
}

# File containing LCA data
lca_path = "lca.tsv"  # Update with actual file path

# Import LCA Data
print(f"Importing data from {lca_path}")

# Load and preprocess LCA data
lca = pd.read_csv(lca_path, sep="\t", usecols=[
    "og_id", "tech", "seq_date", "code", "annotation", "taxonomy", "lca", "percent_match",
    "length", "lca_run_date", "region"
])

# Convert integer-like columns to actual integers (handling missing values)
int_columns = ["length"]
for col in int_columns:
    lca[col] = pd.to_numeric(lca[col], errors='coerce').astype('Int64')  # Convert to nullable Int64

# Replace NaN values with None for database compatibility
lca = lca.replace({pd.NA: None})

try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    row_count = 0  # Track number of processed rows

    for index, row in lca.iterrows():
        values = tuple(row.values)
        # Extract primary key values
        og_id, tech, seq_date, code = row['og_id'], row['tech'], row['seq_date'], row['code']

        # UPSERT: Insert if not exists, otherwise update
        upsert_query = """
        INSERT INTO lca (
            og_id, tech, seq_date, code, annotation, taxonomy, lca, percent_match, length, lca_run_date, region
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (og_id, tech, seq_date, code, annotation,region) 
        DO UPDATE SET
            annotation = EXCLUDED.annotation,
            taxonomy = EXCLUDED.taxonomy,
            lca = EXCLUDED.lca,
            percent_match = EXCLUDED.percent_match,
            length = EXCLUDED.length,
            region = EXCLUDED.region;
        """

        # Debugging Check
        print("\n🔹 DEBUG: Checking SQL Query and Values")
        print("Number of %s placeholders:", upsert_query.count("%s"))
        print("Number of values in row:", len(values))
        print("Values:", values)
        
        # Ensure placeholders match values before executing
        if len(values) != upsert_query.count("%s"):
            print("❌ ERROR: Mismatch in number of placeholders and values!")
            continue  # Skip this row to avoid SQL error

        cursor.execute(upsert_query, values)
        row_count += 1  # Increment successful row counter
    
    print(f"✅ Successfully processed {row_count} rows. Committing changes...")

    conn.commit()
    print("Database update complete.")

except Exception as e:
    print(f"❌ Error: {e}")

finally:
    if cursor:
        cursor.close()
    if conn:
        conn.close()
    print("🔹 Connection closed.")
