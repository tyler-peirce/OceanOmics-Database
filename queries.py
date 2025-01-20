import numpy as np
import math

def insert_sample(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None
    
    
    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Sample" (
        og_id,
        field_id,
        nominal_species_id,
        common_name,
        collector,
        contact,
        date_collected,
        sex,
        weight,
        lengthTL_and_lengthFL,
        country,
        state,
        location,
        latitude_collection,
        longitude_collection,
        depth_collection,
        collection_method,
        preservation_method,
        sample_condition,
        photo_voucher,
        photo_id,
        specimen_voucher,
        voucher_id,
        comments,
        priority,
        tissues,
        extracted,
        extraction_queue,
        ilmn,
        il_status,
        hifi,
        pb_status,
        hic,
        hic_status,
        nano,
        ont_num,
        rna,
        rna_status,
        ilrna,
        ilrna_status,
        assigned_species,
        eschmeyer_id,
        ncbi_sample_name,
        ncbi_biosample_id,
        hifi_lca_outcome,
        ncbi_id,
        tol_id,
        ncbi_bioproject_id_lvl_3_hifi,
        bioproject_id_haplotype_1,
        bioproject_id_haplotype_2,
        bioproject_sequencing_data,
        ncbi_assembly_upload,
        ncbi_raw_reads_upload,
        hifi_public,
        illumina_lca,
        ncbi_bioproject_id_draft,
        illumina_public,
        draft_sra_accessions,
        draft_assembly_accession
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["og_id"],
        row["field_id"],
        row["nominal_species_id"],
        row["common_name"],
        row["collector"],
        row["contact"],
        str(row["date_collected"]) if row["date_collected"] else None,  # Handle NULL for date
        str(row["sex"]) if row["sex"] else None,
        str(row["weight"]) if row["weight"] else None,
        row["lengthTL_and_lengthFL"],
        row["country"],
        row["state"],
        row["location"],
        row["latitude_collection"],
        row["longitude_collection"],
        str(row["depth_collection"]) if row["depth_collection"] else None,
        row["collection_method"],
        row["preservation_method"],
        row["sample_condition"],
        str(row["photo_voucher"]) if row["photo_voucher"] else None,
        row["photo_id"],
        row["specimen_voucher"],
        str(row["voucher_id"]) if row["voucher_id"] else None,
        str(row["comments"]) if row["comments"] else None,
        row["priority"],
        str(row["tissues"]) if row["tissues"] else None,
        str(row["extracted"]) if row["extracted"] else None,
        row["extraction_queue"],
        row["ilmn"],
        row["il_status"],
        row["hifi"],
        row["pb_status"],
        row["hic"],
        row["hic_status"],
        row["nano"],
        row["ont_num"],
        row["rna"],
        row["rna_status"],
        row["ilrna"],
        row["ilrna_status"],
        row["assigned_species"],
        str(row["eschmeyer_id"]) if row["eschmeyer_id"] else None,
        row["ncbi_sample_name"],
        row["ncbi_biosample_id"],
        row["hifi_lca_outcome"],
        str(row["ncbi_id"]) if row["ncbi_id"] else None,
        row["tol_id"],
        row["ncbi_bioproject_id_lvl_3_hifi"],
        row["bioproject_id_haplotype_1"],
        row["bioproject_id_haplotype_2"],
        row["bioproject_sequencing_data"],
        str(row["ncbi_assembly_upload"]) if row["ncbi_assembly_upload"] else None,
        row["ncbi_raw_reads_upload"],
        str(row["hifi_public"]) if row["hifi_public"] else None,
        row["illumina_lca"],
        row["ncbi_bioproject_id_draft"],
        str(row["illumina_public"]) if row["illumina_public"] else None,
        str(row["draft_sra_accessions"]) if row["draft_sra_accessions"] else None,
        str(row["draft_assembly_accession"]) if row["draft_assembly_accession"] else None,
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for Sample:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for Sample.")
    except Exception as e:
        print(f"Error during INSERT for Sample: {e}")







def insert_tissue(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None


    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Tissue" (
        tissue_id,
        og_id,
        field_id,
        alt_id,
        tissue,
        extracted,
        freezer,
        shelf,
        rack,
        level,
        box,
        comment
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["tissue_id"],
        row["og_id"],
        row["field_id"],
        row["alt_id"],
        row["tissue"],
        row["extracted"],
        row["freezer"],
        row["shelf"],
        row["rack"],
        row["level"],
        row["box"],
        row["comment"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for Tissue:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for Tissue.")
    except Exception as e:
        print(f"Error during INSERT for Tissue: {e}")







def insert_rnaextraction(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "RNA_Extraction" (
        rna_id,
        tissue_id,
        ext_num,
        status,
        extraction_method,
        extraction_date,
        extraction_batch_id,
        final_buffer,
        volume,
        qubit_conc,
        nano_drop_conc,
        ratio_260_280,
        ratio_260_230,
        total_yield,
        tapestation_id,
        rna_dv200,
        rin,
        extraction_qc,
        comment,
        rna_freezer,
        rna_shelf,
        rna_rack,
        rna_level,
        rna_box,
        rna_notes
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["rna_id"],
        row["tissue_id"],
        row["ext_num"],
        row["status"],
        row["extraction_method"],
        str(row["extraction_date"]) if row["extraction_date"] else None,  # Handle NULL for date
        row["extraction_batch_id"],
        row["final_buffer"],
        row["volume"],
        row["qubit_conc"],
        row["nano_drop_conc"],
        row["ratio_260_280"],
        row["ratio_260_230"],
        row["total_yield"],
        row["tapestation_id"],
        row["rna_dv200"],
        row["rin"],
        row["extraction_qc"],
        row["comment"],
        row["rna_freezer"],
        row["rna_shelf"],
        row["rna_rack"],
        row["rna_level"],
        row["rna_box"],
        row["rna_notes"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for RNA_Extraction:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for RNA_Extraction.")
    except Exception as e:
        print(f"Error during INSERT for RNA_Extraction: {e}")







def insert_dnaextraction(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "DNA_Extraction" (
        dna_id,
        tissue_id,
        ext_num,
        status,
        extraction_method,
        extraction_date,
        extraction_batch_id,
        final_buffer,
        volume,
        qubit_conc,
        nano_drop_conc,
        ratio_260_280,
        ratio_260_230,
        ratioqubit_nanodrop,
        total_yield,
        gdna_femtol_id,
        av_size,
        extraction_qc,
        comment,
        dna_freezer,
        dna_shelf,
        dna_rack,
        dna_level,
        dna_box,
        dna_notes
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["dna_id"],
        row["tissue_id"],
        row["ext_num"],
        row["status"],
        row["extraction_method"],
        str(row["extraction_date"]) if row["extraction_date"] else None,  # Handle NULL for date
        row["extraction_batch_id"],
        row["final_buffer"],
        row["volume"],
        row["qubit_conc"],
        row["nano_drop_conc"],
        row["ratio_260_280"],
        row["ratio_260_230"],
        row["ratioqubit_nanodrop"],
        row["total_yield"],
        row["gdna_femtol_id"],
        row["av_size"],
        row["extraction_qc"],
        row["comment"],
        row["dna_Freezer"],
        row["dna_Shelf"],
        row["dna_Rack"],
        row["dna_Level"],
        row["dna_Box"],
        row["dna_Notes"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for DNA_Extraction:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for DNA_Extraction.")
    except Exception as e:
        print(f"Error during INSERT for DNA_Extraction: {e}")








def insert_hiclysate(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "HiC_Lysate" (
        lysate_id,
        tissue_id,
        lysate_num,
        lysate_status,
        lysate_prep_date,
        lysate_batch_id,
        lysate_conc,
        total_lysate,
        lysate_cde,
        lysate_comments
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["lysate_id"],
        row["tissue_id"],
        row["lysate_num"],
        row["lysate_status"],
        str(row["lysate_prep_date"]) if row["lysate_prep_date"] else None,  # Handle NULL for date
        row["lysate_batch_id"],
        row["lysate_conc"],
        row["total_lysate"],
        row["lysate_cde"],
        row["lysate_comments"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for HiC_Lysate:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for HiC_Lysate.")
    except Exception as e:
        print(f"Error during INSERT for HiC_Lysate: {e}")









def insert_rnalibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "RNA_Library" (
        rna_library_tube_id,
        rna_id,
        rna_num,
        rna_status,
        library_method,
        library_date,
        library_id,
        library_size,
        perc_product,
        library_qubit_conc,
        library_molarity,
        index_set,
        index_well,
        index_inx,
        kinnex_primers,
        kinnex_barcode,
        il_comments
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["rna_library_tube_id"],
        row["rna_id"],
        row["rna_num"],
        row["rna_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for date
        row["library_id"],
        row["library_size"],
        row["perc_product"],
        row["library_qubit_conc"],
        row["library_molarity"],
        row["index_set"],
        row["index_well"],
        row["index_inx"],
        row["kinnex_primers"],
        row["kinnex_barcode"],
        row["il_comments"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for RNA_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for RNA_Library.")
    except Exception as e:
        print(f"Error during INSERT for RNA_Library: {e}")







def insert_illuminalibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Illumina_Library" (
        illumina_library_tube_id,
        dna_id,
        ilmn_num,
        ilmn_status,
        library_method,
        library_date,
        library_id,
        index_set,
        index_well,
        index_idx,
        library_qubit_conc,
        il_comments
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["illumina_library_tube_id"],
        row["dna_id"],
        row["ilmn_num"],
        row["ilmn_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for date
        row["library_id"],
        row["index_set"],
        row["index_well"],
        row["index_idx"],
        row["library_qubit_conc"],
        row["il_comments"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for Illumina_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for Illumina_Library.")
    except Exception as e:
        print(f"Error during INSERT for Illumina_Library: {e}")







def insert_ontlibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "ONT_Library" (
        ont_library_tube_id,
        dna_id,
        ont_num,
        ont_status,
        library_date,
        library_id,
        library_method,
        library_type,
        est_loading_size,
        ont_comments
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["ont_library_tube_id"],
        row["dna_id"],
        row["ont_num"],
        row["ont_status"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for date
        row["library_id"],
        row["library_method"],
        row["library_type"],
        row["est_loading_size"],
        row["ont_comments"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for ONT_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for ONT_Library.")
    except Exception as e:
        print(f"Error during INSERT for ONT_Library: {e}")








def insert_pacbiolibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "PacBio_Library" (
        pacbio_library_tube_id,
        dna_id,
        pacb_num,
        pacb_status,
        library_method,
        library_date,
        library_id,
        dna_treatment,
        index_well,
        barcode,
        shear_femtol_id,
        shear_av_size,
        seq_femto_id,
        seq_av_size,
        library_conc,
        comment
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["pacbio_library_tube_id"],
        row["dna_id"],
        row["pacb_num"],
        row["pacb_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for date
        row["library_id"],
        row["dna_treatment"],
        row["index_well"],
        row["barcode"],
        row["shear_femtol_id"],
        row["shear_av_size"],
        row["seq_femto_id"],
        row["seq_av_size"],
        row["library_conc"],
        row["comment"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for PacBio_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for PacBio_Library.")
    except Exception as e:
        print(f"Error during INSERT for PacBio_Library: {e}")








def insert_hiclibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "HiC_Library" (
        hic_library_tube_id, 
        lysate_id, 
        hic_num, 
        hic_status, 
        library_method, 
        library_date, 
        library_id, 
        prox_ligation_conc, 
        purified_dna_total, 
        index_set, 
        library_conc, 
        library_size, 
        hic_comments
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["hic_library_tube_id"],
        row["lysate_id"],
        row["hic_num"],
        row["hic_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle potential NULL for dates
        row["library_id"],
        row["prox_ligation_conc"],
        row["purified_dna_total"],
        row["index_set"],
        row["library_conc"],
        row["library_size"],
        row["hic_comments"],
    )

    # Debugging: Print the query and parameters
    print(f"Executing INSERT query for HiC_Library: {sql_query}")
    print(f"With parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for HiC_Library.")
    except Exception as e:
        print(f"Error during INSERT for HiC_Library: {e}")







def insert_sequencing(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None
    # Replace empty strings with None to handle NULL values in the database
    for key in ["rna_library_tube_id", "illumina_library_tube_id", "ont_library_tube_id", "pacbio_library_tube_id", "hic_library_tube_id"]:
        if row[key] == "":
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Sequencing" (
        sequencing_id, 
        rna_library_tube_id, 
        illumina_library_tube_id, 
        ont_library_tube_id, 
        pacbio_library_tube_id, 
        hic_library_tube_id, 
        technology, 
        instrument, 
        run_date, 
        run_id, 
        cell_id, 
        smrt_num, 
        seq_comments, 
        seq_type
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        row["sequencing_id"],
        row["rna_library_tube_id"],
        row["illumina_library_tube_id"],
        row["ont_library_tube_id"],
        row["pacbio_library_tube_id"],
        row["hic_library_tube_id"],
        row["technology"],
        row["instrument"],
        str(row["run_date"]),  # Ensure run_date is a string
        row["run_id"],
        row["cell_id"],
        row["smrt_num"],
        row["seq_comments"],
        row["seq_type"],
    )

    # Debugging: Print the query and parameters
    print(f"Executing INSERT query for Sequencing:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful.")
    except Exception as e:
        print(f"Error during INSERT: {e}")








def insert_draft(cursor, primary_key, sequencing_id, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Draft_Genome" (
        draft_id, 
        sequencing_id, 
        mach, 
        initial, 
        passed_filter_reads, 
        low_quality_reads, 
        too_many_N_reads, 
        too_short_reads, 
        too_long_reads, 
        raw_total_reads, 
        raw_total_bases, 
        raw_q20_bases, 
        raw_q30_bases, 
        raw_q20_rate, 
        raw_q30_rate, 
        raw_read1_mean_length, 
        raw_read2_mean_length, 
        raw_gc_content, 
        total_reads, 
        total_bases, 
        q20_bases, 
        q30_bases, 
        q20_rate, 
        q30_rate, 
        read1_mean_length, 
        read2_mean_length, 
        gc_content, 
        homozygosity, 
        heterozygosity, 
        genome_size, 
        repeat_size, 
        unique_size, 
        model_fit, 
        error_rate, 
        num_contigs, 
        num_contigs_mitochondrion, 
        num_contigs_plastid, 
        num_contigs_prokarya, 
        bp_mitochondrion, 
        bp_plastid, 
        bp_prokarya, 
        complete, 
        single_copy, 
        multi_copy, 
        fragmented, 
        missing, 
        n_markers, 
        domain, 
        number_of_scaffolds, 
        number_of_contigs, 
        total_length, 
        percent_gaps, 
        scaffold_N50, 
        contigs_N50, 
        unique_k_mers_assembly, 
        k_mers_total, 
        qv, 
        error, 
        k_mer_set, 
        solid_k_mers, 
        total_k_mers, 
        completeness, 
        dep_method, 
        adjust, 
        read_bp, 
        map_adjust, 
        sc_depth, 
        est_genome_size
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        primary_key,
        sequencing_id,
        row["mach"],
        row["initial"],
        row["passed_filter_reads"],
        row["low_quality_reads"],
        row["too_many_n_reads"],
        row["too_short_reads"],
        row["too_long_reads"],
        row["raw_total_reads"],
        row["raw_total_bases"],
        row["raw_q20_bases"],
        row["raw_q30_bases"],
        row["raw_q20_rate"],
        row["raw_q30_rate"],
        row["raw_read1_mean_length"],
        row["raw_read2_mean_length"],
        row["raw_gc_content"],
        row["total_reads"],
        row["total_bases"],
        row["q20_bases"],
        row["q30_bases"],
        row["q20_rate"],
        row["q30_rate"],
        row["read1_mean_length"],
        row["read2_mean_length"],
        row["gc_content"],
        row["homozygosity"],
        row["heterozygosity"],
        row["genomesize"],
        row["repeatsize"],
        row["uniquesize"],
        row["modelfit"],
        row["errorrate"],
        row["num_contigs"],
        row["num_contigs_mitochondrion"],
        row["num_contigs_plastid"],
        row["num_contigs_prokarya"],
        row["bp_mitochondrion"],
        row["bp_plastid"],
        row["bp_prokarya"],
        row["complete"],
        row["single_copy"],
        row["multi_copy"],
        row["fragmented"],
        row["missing"],
        row["n_markers"],
        row["domain"],
        row["number_of_scaffolds"],
        row["number_of_contigs"],
        row["total_length"],
        row["percent_gaps"],
        row["scaffold_n50"],
        row["contigs_n50"],
        row["unique_k_mers_assembly"],
        row["k_mers_total"],
        row["qv"],
        row["error"],
        row["k_mer_set"],
        row["solid_k_mers"],
        row["total_k_mers"],
        row["completeness"],
        row["depmethod"],
        row["adjust"],
        row["readbp"],
        row["mapadjust"],
        row["scdepth"],
        row["estgenomesize"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for Draft_Genome:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for Draft_Genome.")
    except Exception as e:
        print(f"Error during INSERT for Draft_Genome: {e}")







def insert_mito(cursor, primary_key, sequencing_id, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Mito_Genome" (
        mito_id,
        sequencing_id,
        tech,
        code,
        stats,
        seq_length,
        _12srna,
        _16srna,
        atp6,
        atp8,
        cox1,
        cox2,
        cox3,
        cytB,
        nd1,
        nd2,
        nd3,
        nd4,
        nd4l,
        nd5,
        nd6,
        bankit_num,
        accession_num
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        primary_key,
        sequencing_id,
        row["tech"],
        row["code"],
        row["stats"],
        row["seq_length"],
        row["12srna"],
        row["16srna"],
        row["ATP6"],
        row["ATP8"],
        row["COX1"],
        row["COX2"],
        row["COX3"],
        row["CYTB"],
        row["ND1"],
        row["ND2"],
        row["ND3"],
        row["ND4"],
        row["ND4L"],
        row["ND5"],
        row["ND6"],
        row["bankit_num"],
        row["accession_num"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for Mito_Genome:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for Mito_Genome.")
    except Exception as e:
        print(f"Error during INSERT for Mito_Genome: {e}")








def insert_lca(cursor, primary_key, mito_id, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "LCA" (
        lca_id,
        mito_id,
        assay,
        database,
        code,
        tech,
        seq_date,
        lca_tax,
        lca,
        percent_identity,
        length_bp
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        primary_key,
        mito_id,
        row["assay"],
        row["database"],
        row["code"],
        row["tech"],
        row["seq_date"],
        row["LCAtax"],
        row["LCA"],
        row["percent_identity"],
        row["length_bp"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for LCA:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for LCA.")
    except Exception as e:
        print(f"Error during INSERT for LCA: {e}")








def insert_species(cursor, primary_key, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Preprocess the row to replace NaN with None
    for key, value in row.items():
        if isinstance(value, float) and math.isnan(value):  # Check if the value is NaN
            row[key] = None

    # Construct the INSERT query
    sql_query = """
    INSERT INTO "Species" (
        species,
        class,
        ordr,
        family,
        genus,
        epithet,
        afd_common_name,
        family_common_name,
        ncbi_taxon_id,
        synonym,
        specimen_tol_id,
        sequencing_status,
        ont,
        hifi,
        hic,
        draft_sequencing_status,
        illumina,
        draft_genome_bioproject_id,
        genome_available,
        internal_aus_status_fishbase,
        cites_listing,
        iucn_code,
        iucn_assessment,
        iucn_dateassessed,
        epbc,
        internal_first_in_family,
        internal_first_in_genus,
        internal_conservation_value,
        internal_research,
        internal_endemic,
        sequencing_priority,
        collaboration,
        comments,
        lab_database_status
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
    """

    # Prepare parameters
    params = (
        primary_key,
        row["class"],
        row["order"],
        row["family"],
        row["genus"],
        row["epithet"],
        row["AFD_common_name"],
        row["Family_Common_Name"],
        row["ncbi_taxon_id"] if row["ncbi_taxon_id"] else None,
        row["synonym"],
        row["Specimen ToLID"],
        row["sequencing_status"],
        row["ONT"],
        row["HiFi"],
        row["HiC"],
        row["DRAFT_sequencing_status"],
        row["Illumina"],
        row["Draft_Genome_Bioproject_ID"],
        row["Genome_Available"],
        row["INTERNAL_Aus_Status_FishBase"],
        row["CITES Listing"],
        row["IUCN_Code"],
        row["IUCN_Assessment"],
        str(row["IUCN_DateAssessed"]) if row["IUCN_DateAssessed"] else None,  # Handle NULL for dates
        row["EPBC"],
        row["INTERNAL_First in Family"],
        row["INTERNAL_First in Genus"],
        row["INTERNAL_ Conservation value"],
        row["INTERNAL_ Research"],
        row["INTERNAL_Endemic"],
        row["Sequencing priority Immediate / High / Medium / Draft /Hold"],
        row["Collaboration"],
        row["Comments"],
        row["Lab_database_status"],
    )

    # Debugging: Log the query and parameters
    print(f"Executing INSERT query for Species:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Insert successful for Species.")
    except Exception as e:
        print(f"Error during INSERT for Species: {e}")




# UPDATE QUERIES

def update_sample(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the base SQL query
    sql_query = """
    UPDATE "Sample"
    SET
        og_id = %s,
        field_id = %s,
        nominal_species_id = %s,
        common_name = %s,
        collector = %s,
        contact = %s,
        date_collected = %s,
        sex = %s,
        weight = %s,
        lengthTL_and_lengthFL = %s,
        country = %s,
        state = %s,
        location = %s,
        latitude_collection = %s,
        longitude_collection = %s,
        depth_collection = %s,
        collection_method = %s,
        preservation_method = %s,
        sample_condition = %s,
        photo_voucher = %s,
        photo_id = %s,
        specimen_voucher = %s,
        voucher_id = %s,
        comments = %s,
        priority = %s,
        tissues = %s,
        extracted = %s,
        extraction_queue = %s,
        ilmn = %s,
        il_status = %s,
        hifi = %s,
        pb_status = %s,
        hic = %s,
        hic_status = %s,
        nano = %s,
        ont_num = %s,
        rna = %s,
        rna_status = %s,
        ilrna = %s,
        ilrna_status = %s,
        assigned_species = %s,
        eschmeyer_id = %s,
        ncbi_sample_name = %s,
        ncbi_biosample_id = %s,
        hifi_lca_outcome = %s,
        ncbi_id = %s,
        tol_id = %s,
        ncbi_bioproject_id_lvl_3_hifi = %s,
        bioproject_id_haplotype_1 = %s,
        bioproject_id_haplotype_2 = %s,
        bioproject_sequencing_data = %s,
        ncbi_assembly_upload = %s,
        ncbi_raw_reads_upload = %s,
        hifi_public = %s,
        illumina_lca = %s,
        ncbi_bioproject_id_draft = %s,
        illumina_public = %s,
        draft_sra_accessions = %s,
        draft_assembly_accession = %s
    WHERE og_id = %s;
    """

    # Prepare parameters
    params = (
        row["og_id"],
        row["field_id"],
        row["nominal_species_id"],
        row["common_name"],
        row["collector"],
        row["contact"],
        str(row["date_collected"]) if row["date_collected"] else None,  # Handle NULL for date
        str(row["sex"]) if row["sex"] else None,
        str(row["weight"]) if row["weight"] else None,
        row["lengthTL_and_lengthFL"],
        row["country"],
        row["state"],
        row["location"],
        row["latitude_collection"],
        row["longitude_collection"],
        str(row["depth_collection"]) if row["depth_collection"] else None,
        row["collection_method"],
        row["preservation_method"],
        row["sample_condition"],
        str(row["photo_voucher"]) if row["photo_voucher"] else None,
        row["photo_id"],
        row["specimen_voucher"],
        str(row["voucher_id"]) if row["voucher_id"] else None,
        str(row["comments"]) if row["comments"] else None,
        row["priority"],
        str(row["tissues"]) if row["tissues"] else None,
        str(row["extracted"]) if row["extracted"] else None,
        row["extraction_queue"],
        row["ilmn"],
        row["il_status"],
        row["hifi"],
        row["pb_status"],
        row["hic"],
        row["hic_status"],
        row["nano"],
        row["ont_num"],
        row["rna"],
        row["rna_status"],
        row["ilrna"],
        row["ilrna_status"],
        row["assigned_species"],
        str(row["eschmeyer_id"]) if row["eschmeyer_id"] else None,
        row["ncbi_sample_name"],
        row["ncbi_biosample_id"],
        row["hifi_lca_outcome"],
        str(row["ncbi_id"]) if row["ncbi_id"] else None,
        row["tol_id"],
        row["ncbi_bioproject_id_lvl_3_hifi"],
        row["bioproject_id_haplotype_1"],
        row["bioproject_id_haplotype_2"],
        row["bioproject_sequencing_data"],
        str(row["ncbi_assembly_upload"]) if row["ncbi_assembly_upload"] else None,
        row["ncbi_raw_reads_upload"],
        str(row["hifi_public"]) if row["hifi_public"] else None,
        row["illumina_lca"],
        row["ncbi_bioproject_id_draft"],
        str(row["illumina_public"]) if row["illumina_public"] else None,
        str(row["draft_sra_accessions"]) if row["draft_sra_accessions"] else None,
        str(row["draft_assembly_accession"]) if row["draft_assembly_accession"] else None,
        row["og_id"],  # This is the WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Sample:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    #print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Sample.")
    except Exception as e:
        print(f"Error during UPDATE for Sample: {e}")







def update_tissue(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = """
    UPDATE "Tissue"
    SET
        tissue_id = %s,
        og_id = %s,
        field_id = %s,
        alt_id = %s,
        tissue = %s,
        extracted = %s,
        freezer = %s,
        shelf = %s,
        rack = %s,
        level = %s,
        box = %s,
        comment = %s
    WHERE tissue_id = %s;
    """

    # Prepare parameters
    params = (
        row["tissue_id"],
        row["og_id"],
        row["field_id"],
        row["alt_id"],
        row["tissue"],
        row["extracted"],
        row["freezer"],
        row["shelf"],
        row["rack"],
        row["level"],
        row["box"],
        row["comment"],
        row["tissue_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Tissue:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    #print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Tissue.")
    except Exception as e:
        print(f"Error during UPDATE for Tissue: {e}")







def update_rnaextraction(cursor, update_query_string, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = f"""
    UPDATE "RNA_Extraction"
    {update_query_string}
    WHERE rna_id = %s;
    """

    # Prepare parameters
    params = (
        row["rna_id"],
        row["tissue_id"],
        row["ext_num"],
        row["status"],
        row["extraction_method"],
        str(row["extraction_date"]) if row["extraction_date"] else None,  # Handle NULL for dates
        row["extraction_batch_id"],
        row["final_buffer"],
        row["volume"],
        row["qubit_conc"],
        row["nano_drop_conc"],
        row["ratio_260_280"],
        row["ratio_260_230"],
        row["total_yield"],
        row["tapestation_id"],
        row["rna_dv200"],
        row["rin"],
        row["extraction_qc"],
        row["comment"],
        row["rna_freezer"],
        row["rna_shelf"],
        row["rna_rack"],
        row["rna_level"],
        row["rna_box"],
        row["rna_notes"],
        row["rna_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for RNA_Extraction:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for RNA_Extraction.")
    except Exception as e:
        print(f"Error during UPDATE for RNA_Extraction: {e}")







def update_dnaextraction(cursor, update_query_string, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = f"""
    UPDATE "DNA_Extraction"
    {update_query_string}
    WHERE dna_id = %s;
    """

    # Prepare parameters
    params = (
        row["dna_id"],
        row["tissue_id"],
        row["ext_num"],
        row["status"],
        row["extraction_method"],
        str(row["extraction_date"]) if row["extraction_date"] else None,  # Handle NULL for dates
        row["extraction_batch_id"],
        row["final_buffer"],
        row["volume"],
        row["qubit_conc"],
        row["nano_drop_conc"],
        row["ratio_260_280"],
        row["ratio_260_230"],
        row["ratioqubit_nanodrop"],
        row["total_yield"],
        row["gdna_femtol_id"],
        row["av_size"],
        row["extraction_qc"],
        row["comment"],
        row["dna_Freezer"],
        row["dna_Shelf"],
        row["dna_Rack"],
        row["dna_Level"],
        row["dna_Box"],
        row["dna_Notes"],
        row["dna_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for DNA_Extraction:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for DNA_Extraction.")
    except Exception as e:
        print(f"Error during UPDATE for DNA_Extraction: {e}")







def update_hiclysate(cursor, update_query_string, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = f"""
    UPDATE HiC_Lysate
    {update_query_string}
    WHERE lysate_id = %s;
    """

    # Prepare parameters
    params = (
        row["lysate_id"],
        row["tissue_id"],
        row["lysate_num"],
        row["lysate_status"],
        str(row["lysate_prep_date"]) if row["lysate_prep_date"] else None,  # Handle NULL for dates
        row["lysate_batch_id"],
        row["lysate_conc"],
        row["total_lysate"],
        row["lysate_cde"],
        row["lysate_comments"],
        row["lysate_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for HiC_Lysate:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for HiC_Lysate.")
    except Exception as e:
        print(f"Error during UPDATE for HiC_Lysate: {e}")







def update_rnalibrary(cursor, update_query_string, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = f"""
    UPDATE RNA_Library
    {update_query_string}
    WHERE rna_library_tube_id = %s;
    """

    # Prepare parameters
    params = (
        row["rna_library_tube_id"],
        row["rna_id"],
        row["rna_num"],
        row["rna_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for dates
        row["library_id"],
        row["library_size"],
        row["perc_product"],
        row["library_qubit_conc"],
        row["library_molarity"],
        row["index_set"],
        row["index_well"],
        row["index_inx"],
        row["kinnex_primers"],
        row["kinnex_barcode"],
        row["il_comments"],
        row["rna_library_tube_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for RNA_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for RNA_Library.")
    except Exception as e:
        print(f"Error during UPDATE for RNA_Library: {e}")








def update_illuminalibrary(cursor, update_query_string, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = f"""
    UPDATE Illumina_Library
    {update_query_string}
    WHERE illumina_library_tube_id = %s;
    """

    # Prepare parameters
    params = (
        row["illumina_library_tube_id"],
        row["dna_id"],
        row["ilmn_num"],
        row["ilmn_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for dates
        row["library_id"],
        row["index_set"],
        row["index_well"],
        row["index_idx"],
        row["library_qubit_conc"],
        row["il_comments"],
        row["illumina_library_tube_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Illumina_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Illumina_Library.")
    except Exception as e:
        print(f"Error during UPDATE for Illumina_Library: {e}")








def update_ontlibrary(cursor, update_query_string, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the SQL query
    sql_query = f"""
    UPDATE ONT_Library
    {update_query_string}
    WHERE ont_library_tube_id = %s;
    """

    # Prepare parameters
    params = (
        row["ont_library_tube_id"],
        row["dna_id"],
        row["ont_num"],
        row["ont_status"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for dates
        row["library_id"],
        row["library_method"],
        row["library_type"],
        row["est_loading_size"],
        row["ont_comments"],
        row["ont_library_tube_id"],  # WHERE clause parameter
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for ONT_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")
    print(f"Number of placeholders: {update_query_string.count('%s')}")
    print(f"Number of parameters: {len(params)}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for ONT_Library.")
    except Exception as e:
        print(f"Error during UPDATE for ONT_Library: {e}")







def update_pacbiolibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the UPDATE query
    sql_query = """
    UPDATE "PacBio_Library"
    SET
        pacbio_library_tube_id = %s,
        dna_id = %s,
        pacb_num = %s,
        pacb_status = %s,
        library_method = %s,
        library_date = %s,
        library_id = %s,
        dna_treatment = %s,
        index_well = %s,
        barcode = %s,
        shear_femtol_id = %s,
        shear_av_size = %s,
        seq_femto_id = %s,
        seq_av_size = %s,
        library_conc = %s,
        comment = %s
    WHERE pacbio_library_tube_id = %s;
    """

    # Prepare parameters
    params = (
        row["pacbio_library_tube_id"],
        row["dna_id"],
        row["pacb_num"],
        row["pacb_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for date
        row["library_id"],
        row["dna_treatment"],
        row["index_well"],
        row["barcode"],
        row["shear_femtol_id"],
        row["shear_av_size"],
        row["seq_femto_id"],
        row["seq_av_size"],
        row["library_conc"],
        row["comment"],
        row["pacbio_library_tube_id"],  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for PacBio_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for PacBio_Library.")
    except Exception as e:
        print(f"Error during UPDATE for PacBio_Library: {e}")







def update_hiclibrary(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the UPDATE query
    sql_query = """
    UPDATE "HiC_Library"
    SET
        hic_library_tube_id = %s,
        lysate_id = %s,
        hic_num = %s,
        hic_status = %s,
        library_method = %s,
        library_date = %s,
        library_id = %s,
        prox_ligation_conc = %s,
        purified_dna_total = %s,
        index_set = %s,
        library_conc = %s,
        library_size = %s,
        hic_comments = %s
    WHERE hic_library_tube_id = %s;
    """

    # Prepare parameters
    params = (
        row["hic_library_tube_id"],
        row["lysate_id"],
        row["hic_num"],
        row["hic_status"],
        row["library_method"],
        str(row["library_date"]) if row["library_date"] else None,  # Handle NULL for date
        row["library_id"],
        row["prox_ligation_conc"],
        row["purified_dna_total"],
        row["index_set"],
        row["library_conc"],
        row["library_size"],
        row["hic_comments"],
        row["hic_library_tube_id"],  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for HiC_Library:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for HiC_Library.")
    except Exception as e:
        print(f"Error during UPDATE for HiC_Library: {e}")







def update_sequencing(cursor, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None
    # Replace empty strings with None to handle NULL values in the database
    for key in ["rna_library_tube_id", "illumina_library_tube_id", "ont_library_tube_id", "pacbio_library_tube_id", "hic_library_tube_id"]:
        if row[key] == "":
            row[key] = None

    # Construct the UPDATE query
    sql_query = """
    UPDATE "Sequencing"
    SET
        sequencing_id = %s,
        rna_library_tube_id = %s,
        illumina_library_tube_id = %s,
        ont_library_tube_id = %s,
        pacbio_library_tube_id = %s,
        hic_library_tube_id = %s,
        technology = %s,
        instrument = %s,
        run_date = %s,
        run_id = %s,
        cell_id = %s,
        smrt_num = %s,
        seq_comments = %s,
        seq_type = %s
    WHERE sequencing_id = %s;
    """

    # Prepare parameters
    params = (
        row["sequencing_id"],
        row["rna_library_tube_id"],
        row["illumina_library_tube_id"],
        row["ont_library_tube_id"],
        row["pacbio_library_tube_id"],
        row["hic_library_tube_id"],
        row["technology"],
        row["instrument"],
        str(row["run_date"]) if row["run_date"] else None,  # Handle NULL for date
        row["run_id"],
        row["cell_id"],
        row["smrt_num"],
        row["seq_comments"],
        row["seq_type"],
        row["sequencing_id"],  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Sequencing:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Sequencing.")
    except Exception as e:
        print(f"Error during UPDATE for Sequencing: {e}")







def update_draft(cursor, primary_key, sequencing_id, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the UPDATE query
    sql_query = """
    UPDATE "Draft_genome"
    SET
        draft_id = %s,
        sequencing_id = %s,
        mach = %s,
        initial = %s,
        passed_filter_reads = %s,
        low_quality_reads = %s,
        too_many_N_reads = %s,
        too_short_reads = %s,
        too_long_reads = %s,
        raw_total_reads = %s,
        raw_total_bases = %s,
        raw_q20_bases = %s,
        raw_q30_bases = %s,
        raw_q20_rate = %s,
        raw_q30_rate = %s,
        raw_read1_mean_length = %s,
        raw_read2_mean_length = %s,
        raw_gc_content = %s,
        total_reads = %s,
        total_bases = %s,
        q20_bases = %s,
        q30_bases = %s,
        q20_rate = %s,
        q30_rate = %s,
        read1_mean_length = %s,
        read2_mean_length = %s,
        gc_content = %s,
        homozygosity = %s,
        heterozygosity = %s,
        genome_size = %s,
        repeat_size = %s,
        unique_size = %s,
        model_fit = %s,
        error_rate = %s,
        num_contigs = %s,
        num_contigs_mitochondrion = %s,
        num_contigs_plastid = %s,
        num_contigs_prokarya = %s,
        bp_mitochondrion = %s,
        bp_plastid = %s,
        bp_prokarya = %s,
        complete = %s,
        single_copy = %s,
        multi_copy = %s,
        fragmented = %s,
        missing = %s,
        n_markers = %s,
        domain = %s,
        number_of_scaffolds = %s,
        number_of_contigs = %s,
        total_length = %s,
        percent_gaps = %s,
        scaffold_N50 = %s,
        contigs_N50 = %s,
        unique_k_mers_assembly = %s,
        k_mers_total = %s,
        qv = %s,
        error = %s,
        k_mer_set = %s,
        solid_k_mers = %s,
        total_k_mers = %s,
        completeness = %s,
        dep_method = %s,
        adjust = %s,
        read_bp = %s,
        map_adjust = %s,
        sc_depth = %s,
        est_genome_size = %s
    WHERE draft_id = %s;
    """

    # Prepare parameters
    params = (
        primary_key,
        sequencing_id,
        row["mach"],
        row["initial"],
        row["passed_filter_reads"],
        row["low_quality_reads"],
        row["too_many_n_reads"],
        row["too_short_reads"],
        row["too_long_reads"],
        row["raw_total_reads"],
        row["raw_total_bases"],
        row["raw_q20_bases"],
        row["raw_q30_bases"],
        row["raw_q20_rate"],
        row["raw_q30_rate"],
        row["raw_read1_mean_length"],
        row["raw_read2_mean_length"],
        row["raw_gc_content"],
        row["total_reads"],
        row["total_bases"],
        row["q20_bases"],
        row["q30_bases"],
        row["q20_rate"],
        row["q30_rate"],
        row["read1_mean_length"],
        row["read2_mean_length"],
        row["gc_content"],
        row["homozygosity"],
        row["heterozygosity"],
        row["genome_size"],
        row["repeat_size"],
        row["unique_size"],
        row["model_fit"],
        row["error_rate"],
        row["num_contigs"],
        row["num_contigs_mitochondrion"],
        row["num_contigs_plastid"],
        row["num_contigs_prokarya"],
        row["bp_mitochondrion"],
        row["bp_plastid"],
        row["bp_prokarya"],
        row["complete"],
        row["single_copy"],
        row["multi_copy"],
        row["fragmented"],
        row["missing"],
        row["n_markers"],
        row["domain"],
        row["number_of_scaffolds"],
        row["number_of_contigs"],
        row["total_length"],
        row["percent_gaps"],
        row["scaffold_N50"],
        row["contigs_N50"],
        row["unique_k_mers_assembly"],
        row["k_mers_total"],
        row["qv"],
        row["error"],
        row["k_mer_set"],
        row["solid_k_mers"],
        row["total_k_mers"],
        row["completeness"],
        row["dep_method"],
        row["adjust"],
        row["read_bp"],
        row["map_adjust"],
        row["sc_depth"],
        row["est_genome_size"],
        primary_key,  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Draft_genome:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Draft_genome.")
    except Exception as e:
        print(f"Error during UPDATE for Draft_genome: {e}")








def update_mito(cursor, primary_key, sequencing_id, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the UPDATE query
    sql_query = """
    UPDATE "Mito_genome"
    SET
        mito_id = %s,
        sequencing_id = %s,
        tech = %s,
        code = %s,
        stats = %s,
        seq_length = %s,
        _12srna = %s,
        _16srna = %s,
        atp6 = %s,
        atp8 = %s,
        cox1 = %s,
        cox2 = %s,
        cox3 = %s,
        cytB = %s,
        nd1 = %s,
        nd2 = %s,
        nd3 = %s,
        nd4 = %s,
        nd4l = %s,
        nd5 = %s,
        nd6 = %s,
        bankit_num = %s,
        accession_num = %s
    WHERE mito_id = %s;
    """

    # Prepare parameters
    params = (
        primary_key,
        sequencing_id,
        row["tech"],
        row["code"],
        row["stats"],
        row["seq_length"],
        row["12srna"],
        row["16srna"],
        row["ATP6"],
        row["ATP8"],
        row["COX1"],
        row["COX2"],
        row["COX3"],
        row["CYTB"],
        row["ND1"],
        row["ND2"],
        row["ND3"],
        row["ND4"],
        row["ND4L"],
        row["ND5"],
        row["ND6"],
        row["bankit_num"],
        row["accession_num"],
        primary_key,  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Mito_genome:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Mito_genome.")
    except Exception as e:
        print(f"Error during UPDATE for Mito_genome: {e}")








def update_lca(cursor, primary_key, mito_id, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Construct the UPDATE query
    sql_query = """
    UPDATE "LCA"
    SET
        lca_id = %s,
        mito_id = %s,
        assay = %s,
        database = %s,
        code = %s,
        tech = %s,
        seq_date = %s,
        lca_tax = %s,
        lca = %s,
        percent_identity = %s,
        length_bp = %s
    WHERE lca_id = %s;
    """

    # Prepare parameters
    params = (
        primary_key,
        mito_id,
        row["assay"],
        row["database"],
        row["code"],
        row["tech"],
        str(row["seq_date"]) if row["seq_date"] else None,  # Handle NULL for dates
        row["LCAtax"],
        row["LCA"],
        row["percent_identity"],
        row["length_bp"],
        primary_key,  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for LCA:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for LCA.")
    except Exception as e:
        print(f"Error during UPDATE for LCA: {e}")







def update_species(cursor, primary_key, row):
    # Replace `nan` values with `None`
    for key, value in row.items():
        if value != value:  # NaN values are not equal to themselves
            row[key] = None

    # Preprocess the row to replace NaN with None
    for key, value in row.items():
        if isinstance(value, float) and math.isnan(value):  # Check if the value is NaN
            row[key] = None
    # Construct the UPDATE query
    sql_query = """
    UPDATE "Species"
    SET 
        species = %s,
        class = %s,
        ordr = %s,
        family = %s,
        genus = %s,
        epithet = %s,
        afd_common_name = %s,
        family_common_name = %s,
        ncbi_taxon_id = %s,
        synonym = %s,
        specimen_tol_id = %s,
        sequencing_status = %s,
        ont = %s,
        hifi = %s,
        hic = %s,
        draft_sequencing_status = %s,
        illumina = %s,
        draft_genome_bioproject_id = %s,
        genome_available = %s,
        internal_aus_status_fishbase = %s,
        cites_listing = %s,
        iucn_code = %s,
        iucn_assessment = %s,
        iucn_dateassessed = %s,
        epbc = %s,
        internal_first_in_family = %s,
        internal_first_in_genus = %s,
        internal_conservation_value = %s,
        internal_research = %s,
        internal_endemic = %s,
        sequencing_priority = %s,
        collaboration = %s,
        comments = %s,
        lab_database_status = %s
    WHERE species = %s;
    """

    # Prepare parameters, replacing empty strings with None
    params = (
        primary_key,
        row.get("class"),
        row.get("order"),
        row.get("family"),
        row.get("genus"),
        row.get("epithet"),
        row.get("AFD_common_name"),
        row.get("Family_Common_Name"),
        row.get("ncbi_taxon_id"),
        row.get("synonym"),
        row.get("Specimen ToLID"),
        row.get("sequencing_status"),
        row.get("ONT"),
        row.get("HiFi"),
        row.get("HiC"),
        row.get("DRAFT_sequencing_status"),
        row.get("Illumina"),
        row.get("Draft_Genome_Bioproject_ID"),
        row.get("Genome_Available"),
        row.get("INTERNAL_Aus_Status_FishBase"),
        row.get("CITES Listing"),
        row.get("IUCN_Code"),
        row.get("IUCN_Assessment"),
        str(row.get("IUCN_DateAssessed")) if row.get("IUCN_DateAssessed") else None,  # Handle NULL for dates
        row.get("EPBC"),
        row.get("INTERNAL_First in Family"),
        row.get("INTERNAL_First in Genus"),
        row.get("INTERNAL_ Conservation value"),
        row.get("INTERNAL_ Research"),
        row.get("INTERNAL_Endemic"),
        row.get("Sequencing priority Immediate / High / Medium / Draft /Hold"),
        row.get("Collaboration"),
        row.get("Comments"),
        row.get("Lab_database_status"),
        primary_key,  # WHERE clause
    )

    # Debugging: Log the query and parameters
    print(f"Executing UPDATE query for Species:")
    print(f"Query: {sql_query}")
    print(f"Parameters: {params}")
    print(f"Parameter types: {[type(p) for p in params]}")

    # Execute the query
    try:
        cursor.execute(sql_query, params)
        print("Update successful for Species.")
    except Exception as e:
        print(f"Error during UPDATE for Species: {e}")
