try:
    import name_convert
    import schema
    import queries
    import pandas as pd
    import psycopg2
    import os.path
    import string
    from functools import reduce
    from psycopg2 import sql

    conn = None
    cursor = None

    # Database connection parameters
    db_params = {
            'dbname': 'oceanomics',
            'user': 'postgres',
            'password': 'oceanomics',
            'host': '115.146.85.41',  # Change to your host if necessary
            'port': 5432          # Default PostgreSQL port
        }

    master_species_path = "data/Master Species List - Marine Vertebrates.xlsx"

    print("Connecting to PostgreSQL database...")
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()


    
    # Import Master Species Data
    print(f"Importing data from {master_species_path}")

    # Load the Excel file
    species_df = pd.read_excel(master_species_path, sheet_name="Sheet1")

    # Clean the DataFrame to remove non-breaking spaces globally
    species_df = species_df.applymap(lambda x: x.replace('\xa0', '').strip() if isinstance(x, str) else x)

    # Iterate over the DataFrame
    for index, row in species_df.iterrows():
        species = row["species"]
        og_ids = row["OG_number"]

        # Check if the species exists in the Species table
        cursor.execute(
            sql.SQL('SELECT species FROM "Species" WHERE species = %s;'),
            (species,)
        )
        result = cursor.fetchone()

        if result:
            # Update the row if it exists
            queries.update_species(cursor, species, row)
        else:
            # Insert the row if it doesn't exist
            queries.insert_species(cursor, species, row)
    # Commit the transaction to save the changes
    conn.commit()

