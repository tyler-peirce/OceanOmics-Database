import psycopg2

# PostgreSQL connection parameters
db_params = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '115.146.85.41',  # Change to your database host if necessary
    'port': 5432  # Default PostgreSQL port
}

# Function to update the database
def update_database_from_file(file_path):
    try:
        # Connect to PostgreSQL
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()

        # Read the input file
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split("\t")  # Tab-separated
                if len(parts) != 2:
                    print(f"Skipping invalid line: {line.strip()}")
                    continue

                # Extract relevant values
                bankit_id, composite_key, pp_number = parts[0], parts[0].split()[1], parts[1]

                # Split composite_key into individual fields
                key_parts = composite_key.split(".")
                if len(key_parts) < 4:
                    print(f"Skipping invalid composite key: {composite_key}")
                    continue

                og_id, tech, date, code = key_parts[0], key_parts[1], key_parts[2], key_parts[3]

                # Perform the database update
                update_query = """
                    UPDATE your_table_name
                    SET pp_number_column = %s
                    WHERE OG_id = %s AND tech = %s AND date = %s AND code = %s;
                """
                cursor.execute(update_query, (pp_number, og_id, tech, date, code))

        # Commit the transaction
        conn.commit()
        print("Database update complete.")

    except Exception as e:
        print(f"Error: {e}")

    finally:
        # Close database connection
        if cursor:
            cursor.close()
        if conn:
            conn.close()

# Call the function with your input file
update_database_from_file("your_input_file.txt")
