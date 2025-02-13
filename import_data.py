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

    excel_path = "OceanGenomes_database_250213.xlsx"
    draft_path = ""
    mito_path = "mtdnastat.250211.tsv"
    protein_path = ""
    accession_path = "Genbank_accession_numbers_240605.tsv"
    lca_paths = "lca.tsv"
    tolid_path = "TOLID_250211.xlsx"
    master_species_path = "Master_species_list250211.xlsx"
    pacbio_qc_path = ""

    print("Connecting to PostgreSQL database...")
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()


    
    # Import the lab excel database
    print("Importing data from " + excel_path)
    xls = pd.ExcelFile(excel_path)
    # Check
    print(f"Sheets in the Excel file: {xls.sheet_names}")

    
    # Loop through the database tables, and get the Excel sheets for that table
    for db_tble_name in name_convert.db_to_excel_tables.keys():
        excel_tble_name = name_convert.db_to_excel_tables[db_tble_name]
        primary_key = schema.primary_keys[db_tble_name]

        # Handle multiple Excel sheets for the database table
        if isinstance(excel_tble_name, list) and len(excel_tble_name) > 1:
            df_list = []
            for curr_table in excel_tble_name:
                primary_key_excel_name = name_convert.db_to_excel_cols[db_tble_name][primary_key]
                df_list.append(pd.read_excel(xls, curr_table))

                # Check dataframes
                print(f"Loaded data from sheet: {curr_table}")
                print(f"Columns: {df_list[-1].columns}")
                print(f"First few rows:\n{df_list[-1].head()}")

                # Rename the primary key column for easy merging
                if isinstance(primary_key_excel_name, list) and len(primary_key_excel_name) > 1:
                    for curr_pk_excel_name in primary_key_excel_name:
                        df_list[-1].columns = [primary_key if x == curr_pk_excel_name else x for x in df_list[-1].columns]
                else:
                    df_list[-1].columns = [primary_key if x == primary_key_excel_name else x for x in df_list[-1].columns]

            df = reduce(lambda left, right: pd.merge(left, right, on=[primary_key], how='outer'), df_list)
            #Check dataframes
            print(f"Merged DataFrame for table {db_tble_name}:")
            print(f"Columns: {df.columns}")
            print(f"First few rows:\n{df.head()}")

            # Fix duplication of some columns
            df = df.loc[:, ~df.columns.str.endswith("_x")]
            df.columns = df.columns.str.removesuffix("_y")

        else:
            df = pd.read_excel(xls, excel_tble_name)
            #Checks
            print(f"Loaded data from sheet {excel_tble_name}:")
            print(f"Columns: {df.columns}")
            print(f"First few rows:\n{df.head()}")

        try:
            df = df.sort_values(by=["OG"])
            # Check
            print(f"Data sorted by 'OG' column for table {db_tble_name}")
        except KeyError:
            print(f"KeyError: 'OG' column not found in table {db_tble_name}. Skipping sort.")

        # Rename columns to match database
        for col in df.columns:
            try:
                db_col_name = name_convert.excel_to_db_cols[db_tble_name][col]
                df.columns = [db_col_name if x == col else x for x in df.columns]
            except KeyError:
                pass

        # Handle specific table logic
        if db_tble_name == "Sequencing":
            df["rna_library_tube_id"] = ""
            df["illumina_library_tube_id"] = ""
            df["ont_library_tube_id"] = ""
            df["pacbio_library_tube_id"] = ""
            df["hic_library_tube_id"] = ""
            df["seq_type"] = ""

            # Figure out the foreign key based on the technology used
            for index, row in df.iterrows():
                tech = df.at[index, "technology"]
                if tech == "Hi-C":
                    df.at[index, "hic_library_tube_id"] = df.at[index, "Library Tube ID"]
                    df.at[index, "seq_type"] = "HiC"
                elif tech == "Illumina":
                    if "_D" in df.at[index, "Library Tube ID"]:
                        df.at[index, "illumina_library_tube_id"] = df.at[index, "Library Tube ID"]
                        df.at[index, "seq_type"] = "Illumina"
                    elif "_R" in df.at[index, "Library Tube ID"]:
                        df.at[index, "rna_library_tube_id"] = df.at[index, "Library Tube ID"]
                        df.at[index, "seq_type"] = "RNA"
                    else:
                        raise Exception(df.at[index, "Library Tube ID"] + " is not a valid ID for Illumina.")
                elif tech == "ONT":
                    df.at[index, "ont_library_tube_id"] = df.at[index, "Library Tube ID"]
                    df.at[index, "seq_type"] = "ONT"
                elif tech == "PacBio":
                    df.at[index, "pacbio_library_tube_id"] = df.at[index, "Library Tube ID"]
                    df.at[index, "seq_type"] = "PacBio"
                else:
                    raise Exception(tech + " in Sequencing table is not a valid technology.")

    # The sample table needs the tolid data added
        if db_tble_name == "Sample":
            tolid_df = pd.ExcelFile(tolid_path)
            tolid_df = pd.read_excel(tolid_df, "Sheet1")
            tolid_df = tolid_df.rename(columns={
                "Specimen ID": "og_id",
                "Assigned Species": "assigned_species",
                "Eschmeyer ID": "eschmeyer_id",
                "NCBI_Sample_Name": "ncbi_sample_name",
                "NCBI_BioSample_ID": "ncbi_biosample_id",
                "HIFI LCA outcome": "hifi_lca_outcome",
                "NCBI ID": "ncbi_id",
                "ToLID": "tol_id",
                "NCBI BioProject ID lvl 3 - HiFi": "ncbi_bioproject_id_lvl_3_hifi",
                "BioProject ID_Haplotype 1": "bioproject_id_haplotype_1",
                "BioProject ID_Haplotype 2": "bioproject_id_haplotype_2",
                "BioProject ID_Sequencing data": "bioproject_sequencing_data",
                "NCBI_Assembly_Upload": "ncbi_assembly_upload",
                "NCBI_Raw_Reads_Upload": "ncbi_raw_reads_upload",
                "HiFi - Public": "hifi_public"
            })
            df = pd.merge(df, tolid_df, on = "og_id", how = "left")

        # Start building a string to be used in the SQL INSERT/UPDATE query
        cols_to_update = list(name_convert.db_to_excel_cols[db_tble_name].keys())
        if db_tble_name == "Sample":
            cols_to_update.append("assigned_species")
            cols_to_update.append("eschmeyer_id")
            cols_to_update.append("ncbi_sample_name")
            cols_to_update.append("ncbi_biosample_id")
            cols_to_update.append("hifi_lca_outcome")
            cols_to_update.append("ncbi_id")
            cols_to_update.append("tol_id")
            cols_to_update.append("ncbi_bioproject_id_lvl_3_hifi")
            cols_to_update.append("bioproject_id_haplotype_1")
            cols_to_update.append("bioproject_id_haplotype_2")
            cols_to_update.append("bioproject_sequencing_data")
            cols_to_update.append("ncbi_assembly_upload")
            cols_to_update.append("ncbi_raw_reads_upload")
            cols_to_update.append("hifi_public")
            #cols_to_update.append("illumina_lca")
            #cols_to_update.append("ncbi_bioproject_id_draft")
            #cols_to_update.append("illumina_public")
            #cols_to_update.append("draft_sra_accessions")
            #cols_to_update.append("draft_assembly_accession")

        insert_query_string = ""
        insert_query_string_suffix = ""
        update_query_string = " SET "
        for col in cols_to_update:
            insert_query_string += col + ", "
            insert_query_string_suffix += "%s, "
            update_query_string += col + " = %s, "
        insert_query_string = insert_query_string[0:-2]
        insert_query_string_suffix = insert_query_string_suffix[0:-2]
        update_query_string = update_query_string[0:-2]
        insert_query_string = "(" + insert_query_string + ") VALUES (" + insert_query_string_suffix + ");"
        update_query_string = update_query_string + " WHERE " + primary_key + " = %s;"

        #Loop through rows of df
        for index, row in df.iterrows():
            key_value = df.at[index, primary_key]

            # Debugging: Print the query details
            print(f"Preparing to execute SELECT query for table {db_tble_name}:")
            print(f"Query: SELECT {primary_key} FROM {db_tble_name} WHERE {primary_key} = %s;")
            print(f"Parameter: {str(key_value)}")

            # Check if the row exists
            cursor.execute(
                sql.SQL("SELECT {primary_key} FROM {table} WHERE {primary_key} = %s").format(
                    primary_key=sql.Identifier(primary_key),
                    table=sql.Identifier(db_tble_name)
                ),
                (str(key_value),)
            )

            result = cursor.fetchone()

            if result:
                print(f"Row with primary key {key_value} exists in table {db_tble_name}. Preparing to update.")
                # Row exists, update it
                try:
                    update_funcs = {
                        "Sample": queries.update_sample,
                        "Tissue": queries.update_tissue,
                        "RNA_Extraction": queries.update_rnaextraction,
                        "DNA_Extraction": queries.update_dnaextraction,
                        "HiC_Lysate": queries.update_hiclysate,
                        "RNA_Library": queries.update_rnalibrary,
                        "Illumina_Library": queries.update_illuminalibrary,
                        "ONT_Library": queries.update_ontlibrary,
                        "PacBio_Library": queries.update_pacbiolibrary,
                        "HiC_Library": queries.update_hiclibrary,
                        "Sequencing": queries.update_sequencing,
                    }
                    update_funcs[db_tble_name](cursor, row)
                    conn.commit()  # Save changes after update
                except Exception as e:
                    print(f"Error in update function for {db_tble_name}: {e}")

            else:
                print(f"Row with primary key {key_value} does not exist in table {db_tble_name}. Preparing to insert.")
                # Row does not exist, insert it
                try:
                    insert_funcs = {
                        "Sample": queries.insert_sample,
                        "Tissue": queries.insert_tissue,
                        "RNA_Extraction": queries.insert_rnaextraction,
                        "DNA_Extraction": queries.insert_dnaextraction,
                        "HiC_Lysate": queries.insert_hiclysate,
                        "RNA_Library": queries.insert_rnalibrary,
                        "Illumina_Library": queries.insert_illuminalibrary,
                        "ONT_Library": queries.insert_ontlibrary,
                        "PacBio_Library": queries.insert_pacbiolibrary,
                        "HiC_Library": queries.insert_hiclibrary,
                        "Sequencing": queries.insert_sequencing,
                    }
                    insert_funcs[db_tble_name](cursor, row)
                    conn.commit()  # Save changes after insert
                except Exception as e:
                    print(f"Error in insert function for {db_tble_name}: {e}")

    # Get data from the Sequencing table that's needed for adding the genome data
    draft_sequencing_dict = {}
    mito_sequencing_dict = {}

    # Execute the query to fetch data from the Sequencing table
    cursor.execute(
        sql.SQL('SELECT sequencing_id, technology, run_date, run_id FROM "Sequencing";')
    )
    query_results = cursor.fetchall()

    # Process the query results
    for row in query_results:
        sequencing_id = row[0]
        technology = row[1]
        
        # Map technology to shorthand keys
        if technology == "Illumina":
            tech = "ilmn"
        elif technology == "PacBio":
            tech = "hifi"
        elif technology == "ONT":
            tech = "ont"
        elif technology == "Hi-C":
            tech = "hic"
        else:
            raise ValueError(f"Unknown technology: {technology}")

        # Extract and process values
        run_date = row[3].split("_")[1]  # Assumes `run_date` follows this specific format
        seq_id_split = sequencing_id.split("_")
        og_id = seq_id_split[0].split("-")[0].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        mach = seq_id_split[-3]

        # Create keys for the dictionaries
        draft_key = f"{og_id}_{mach}_{run_date}"
        mito_key = f"{og_id}_{tech}_{run_date}"

        # Populate the dictionaries
        draft_sequencing_dict[draft_key] = sequencing_id
        mito_sequencing_dict[mito_key] = sequencing_id

    # The `draft_sequencing_dict` and `mito_sequencing_dict` are now populated

# Import Draft Genome Data
    print("Importing data from " + draft_path)
    draft = pd.read_csv(draft_path)

    # Adjust column names to match database conventions
    draft.columns = draft.columns.str.replace(".", "_")
    draft.columns = draft.columns.str.lower()

    # Iterate over each row in the DataFrame
    for index, row in draft.iterrows():
        # Construct the primary key
        primary_key = row["sample"].rstrip(string.ascii_uppercase) + "_" + row["mach"] + "_" + str(row["seq_date"])

        # Retrieve sequencing_id using the dictionary
        sequencing_id = draft_sequencing_dict.get(primary_key)
        # Print or log the primary key and its corresponding sequencing_id
        print(f"Primary Key: {primary_key}, Sequencing ID: {sequencing_id}")

        if sequencing_id is None:
            raise ValueError(f"Sequencing ID not found for primary key: {primary_key}")

        # Check if the primary key already exists in Draft_Genome
        cursor.execute(
            sql.SQL("SELECT draft_id FROM Draft_Genome WHERE draft_id = %s;"),
            (primary_key,)
        )
        result = cursor.fetchone()

        if result:
            # Update row if it exists
            queries.update_draft(cursor, primary_key, sequencing_id, row)
        else:
            # Insert row if it does not exist
            queries.insert_draft(cursor, primary_key, sequencing_id, row)

# Import Mitogenome Data
    print(f"Importing data from {mito_path}, {protein_path}, and {accession_path}")

    # Load and preprocess mito data
    mito = pd.read_csv(mito_path, sep="\t")
    mito["code"] = mito["code"].str.split(".").str[0]

    # Load and preprocess protein data
    protein = pd.read_csv(protein_path, sep="\t")
    protein["code"] = protein["code"].str.split(".").str[0]

    # Load and preprocess accession data
    accession = pd.read_csv(accession_path, sep="\t", names=["bankit_num_key", "accession_num"])
    accession["bankit_num"] = accession["bankit_num_key"].str.split(" ").str[0]
    accession["key"] = accession["bankit_num_key"].str.split(" ").str[1]
    accession["OG_id"] = accession["key"].str.split(".").str[0]
    accession["tech"] = accession["key"].str.split(".").str[1]
    accession["seq_date"] = accession["key"].str.split(".").str[2].astype(int)
    accession["code"] = accession["key"].str.split(".").str[3]

    # Merge the datasets
    merged_df = pd.merge(mito, protein, how="outer", on=["OG_id", "tech", "seq_date", "code"])
    merged_df = pd.merge(merged_df, accession, how="outer", on=["OG_id", "tech", "seq_date", "code"])

    # Iterate over the merged DataFrame
    for index, row in merged_df.iterrows():
        # Construct the primary key
        primary_key = f"{row['OG_id']}_{row['tech']}_{row['seq_date']}"

        # Retrieve sequencing_id from mito_sequencing_dict
        sequencing_id = mito_sequencing_dict.get(primary_key)
        if sequencing_id is None:
            raise ValueError(f"Sequencing ID not found for primary key: {primary_key}")

        # Check if the primary key exists in the Mito_Genome table
        cursor.execute(
            sql.SQL("SELECT mito_id FROM Mito_Genome WHERE mito_id = %s;"),
            (primary_key,)
        )
        result = cursor.fetchone()

        if result:
            # Update the row if it exists
            queries.update_mito(cursor, primary_key, sequencing_id, row)
        else:
            # Insert the row if it doesn't exist
            queries.insert_mito(cursor, primary_key, sequencing_id, row)
    
    # Add in the accession numbers
    

# Import LCA Data
    print(f"Importing data from {str(lca_paths)}")

    lca_list = []

    # Process each LCA file
    for lca_path in lca_paths:
        curr_lca = pd.read_csv(lca_path, sep="\t")
        lca_path_split = lca_path.split("_")
        assay = lca_path_split[1]
        database = lca_path_split[2].split(".")[0]
        curr_lca["assay"] = assay
        curr_lca["database"] = database
        lca_list.append(curr_lca)

    # Combine all LCA data into one DataFrame
    lca = pd.concat(lca_list)

    # Iterate through the combined LCA DataFrame
    for index, row in lca.iterrows():
        # Construct the `mito_id` and `primary_key`
        mito_id = f"{row['Sample']}_{row['tech']}_{row['seq_date']}"
        primary_key = f"{mito_id}_{row['assay']}_{row['database']}"

        # Check if the `primary_key` exists in the LCA table
        cursor.execute(
            sql.SQL("SELECT lca_id FROM LCA WHERE lca_id = %s;"),
            (primary_key,)
        )
        result = cursor.fetchone()

        if result:
            # Update the row if it exists
            queries.update_lca(cursor, primary_key, mito_id, row)
        else:
            # Insert the row if it does not exist
            queries.insert_lca(cursor, primary_key, mito_id, row)
    
    
    # I don't have Reference Genome data yet, but Emma said it will be similar to Lauren's draft genome data
    #pacbio_qc_path
    #pacbio_qc_df = pd.ExcelFile(pacbio_qc_path)
    #pacbio_qc_df = pd.read_excel(pacbio_qc_df, "Sheet1")

    # Save changes to database
    conn.commit()
    print("Finished updating database")

    # Close the connection to the database
    cursor.close()
    conn.close()

except Exception as e:
    # Close the connection to the database
    if cursor is not None:
        cursor.close()
    if conn is not None:
        conn.close()

    raise e
