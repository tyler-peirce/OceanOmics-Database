import psycopg2
from psycopg2 import sql

def create_postgresql_db():
    try:
        # PostgreSQL connection parameters
        db_params = {
            'dbname': 'oceanomics',
            'user': 'postgres',
            'password': 'oceanomics',
            'host': '115.146.85.41',  # Change to your host if necessary
            'port': 5432          # Default PostgreSQL port
        }

        # Establish connection to PostgreSQL
        print("Connecting to PostgreSQL database...")
        conn = psycopg2.connect(**db_params)
        cursor = conn.cursor()

        # Species Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "Species" (
            species TEXT PRIMARY KEY,
            class TEXT,
            ordr TEXT,
            family TEXT,
            genus TEXT,
            epithet TEXT,
            afd_common_name TEXT,
            family_common_name TEXT,
            ncbi_taxon_id INTEGER,
            synonym TEXT,
            specimen_tol_id TEXT,
            sequencing_status TEXT,
            ont TEXT,
            hifi TEXT,
            hic TEXT,
            draft_sequencing_status TEXT,
            illumina TEXT,
            draft_genome_bioproject_id TEXT,
            genome_available TEXT,
            internal_aus_status_fishbase TEXT,
            cites_listing TEXT,
            iucn_code TEXT,
            iucn_assessment TEXT,
            iucn_dateassessed TEXT,
            epbc TEXT,
            internal_first_in_family TEXT,
            internal_first_in_genus TEXT,
            internal_conservation_value TEXT,
            internal_research TEXT,
            internal_endemic TEXT,
            sequencing_priority TEXT,
            collaboration TEXT,
            comments TEXT,
            lab_database_status TEXT
            );
            """
        )

        # Sample table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "Sample" (
            og_num INTEGER GENERATED ALWAYS AS (substring(og_id FROM 3)::INTEGER) STORED,
            og_id TEXT PRIMARY KEY,
            field_id TEXT,
            nominal_species_id TEXT,
            common_name TEXT,
            collector TEXT,
            contact TEXT,
            date_collected DATE,
            sex TEXT,
            weight REAL,
            lengthTL_and_lengthFL TEXT,
            country TEXT,
            state TEXT,
            location TEXT,
            latitude_collection TEXT,
            longitude_collection TEXT,
            depth_collection TEXT,
            collection_method TEXT,
            preservation_method TEXT,
            sample_condition TEXT,
            photo_voucher TEXT,
            photo_id TEXT,
            specimen_voucher TEXT,
            voucher_id TEXT,
            comments TEXT,
            priority TEXT,
            tissues TEXT,
            extracted TEXT,
            extraction_queue TEXT,
            ilmn TEXT,
            il_status TEXT,
            hifi TEXT,
            pb_status TEXT,
            hic TEXT,
            hic_status TEXT,
            nano TEXT,
            ont_num TEXT,
            rna TEXT,
            rna_status TEXT,
            ilrna TEXT,
            ilrna_status TEXT,
            assigned_species TEXT,
            eschmeyer_id TEXT,
            ncbi_sample_name TEXT,
            ncbi_biosample_id TEXT,
            hifi_lca_outcome TEXT,
            ncbi_id TEXT,
            tol_id TEXT,
            ncbi_bioproject_id_lvl_3_hifi TEXT,
            bioproject_id_haplotype_1 TEXT,
            bioproject_id_haplotype_2 TEXT,
            bioproject_sequencing_data TEXT,
            ncbi_assembly_upload TEXT,
            ncbi_raw_reads_upload TEXT,
            hifi_public TEXT,
            illumina_lca TEXT,
            ncbi_bioproject_id_draft TEXT,
            illumina_public TEXT,
            draft_sra_accessions TEXT,
            draft_assembly_accession TEXT
            );
            """
        )
    
        # Tissue Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "Tissue" (
            tissue_id TEXT PRIMARY KEY,
            og_id TEXT,
            field_id TEXT,
            alt_id TEXT,
            tissue TEXT,
            extracted INTEGER,
            freezer TEXT,
            shelf INTEGER,
            rack INTEGER,
            level TEXT,
            box TEXT,
            comment TEXT,
            CONSTRAINT fk_og_id FOREIGN KEY (og_id) REFERENCES "Sample" (og_id)
            );
            """
        )

        # RNA_Extraction Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "RNA_Extraction" (
            rna_id TEXT PRIMARY KEY,
            tissue_id TEXT,
            ext_num INTEGER,
            status TEXT,
            extraction_method TEXT,
            extraction_date TEXT,
            extraction_batch_id TEXT,
            final_buffer TEXT,
            volume INTEGER,
            qubit_conc REAL,
            nano_drop_conc REAL,
            ratio_260_280 REAL,
            ratio_260_230 REAL,
            total_yield INTEGER,
            tapestation_id TEXT,
            rna_dv200 REAL,
            rin INTEGER,
            extraction_qc TEXT,
            comment TEXT,
            rna_freezer TEXT,
            rna_shelf TEXT,
            rna_rack TEXT,
            rna_level TEXT,
            rna_box TEXT,
            rna_notes TEXT,
            CONSTRAINT fk_tissue_id FOREIGN KEY (tissue_id) REFERENCES "Tissue" (tissue_id)
            );
            """
        )

        # DNA_Extraction Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "DNA_Extraction" (
            dna_id TEXT PRIMARY KEY,
            tissue_id TEXT,
            ext_num INTEGER,
            status TEXT,
            extraction_method TEXT,
            extraction_date DATE,
            extraction_batch_id TEXT,
            final_buffer TEXT,
            volume INTEGER,
            qubit_conc REAL,
            nano_drop_conc REAL,
            ratio_260_280 REAL,
            ratio_260_230 REAL,
            ratioqubit_nanodrop REAL,
            total_yield INTEGER,
            gdna_femtol_id TEXT,
            av_size INTEGER,
            extraction_qc TEXT,
            comment TEXT,
            dna_freezer TEXT,
            dna_shelf INTEGER,
            dna_rack INTEGER,
            dna_level TEXT,
            dna_box TEXT,
            dna_notes TEXT,
            CONSTRAINT fk_tissue_id FOREIGN KEY (tissue_id) REFERENCES "Tissue" (tissue_id)
            );
            """
        )

        # HiC_Lysate Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "HiC_Lysate" (
            lysate_id TEXT PRIMARY KEY,
            tissue_id TEXT,
            lysate_num INTEGER,
            lysate_status TEXT,
            lysate_prep_date DATE,
            lysate_batch_id TEXT,
            lysate_conc REAL,
            total_lysate REAL,
            lysate_cde REAL,
            lysate_comments TEXT,
            CONSTRAINT fk_tissue_id FOREIGN KEY (tissue_id) REFERENCES "Tissue" (tissue_id)
            );
            """
        )

        # RNA_Library Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "RNA_Library" (
            rna_library_tube_id TEXT PRIMARY KEY,
            rna_id TEXT,
            rna_num INTEGER,
            rna_status TEXT,
            library_method TEXT,
            library_date TEXT,
            library_id TEXT,
            library_size INTEGER,
            perc_product REAL,
            library_qubit_conc REAL,
            library_molarity REAL,
            index_set TEXT,
            index_well TEXT,
            index_inx TEXT,
            kinnex_primers TEXT,
            kinnex_barcode TEXT,
            il_comments TEXT,
            CONSTRAINT fk_rna_id FOREIGN KEY (rna_id) REFERENCES "RNA_Extraction" (rna_id)
            );
            """
        )

        # Illumina_Library Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "Illumina_Library" (
            illumina_library_tube_id TEXT PRIMARY KEY,
            dna_id TEXT,
            ilmn_num INTEGER,
            ilmn_status TEXT,
            library_method TEXT,
            library_date TEXT,
            library_id TEXT,
            index_set TEXT,
            index_well TEXT,
            index_idx TEXT,
            library_qubit_conc REAL,
            il_comments TEXT,
            CONSTRAINT fk_dna_id FOREIGN KEY (dna_id) REFERENCES "DNA_Extraction" (dna_id)
            );
            """
        )

        # ONT_Library Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "ONT_Library" (
            ont_library_tube_id TEXT PRIMARY KEY,
            dna_id TEXT,
            ont_num INTEGER,
            ont_status TEXT,
            library_date TEXT,
            library_id TEXT,
            library_method TEXT,
            library_type TEXT,
            est_loading_size INTEGER,
            ont_comments TEXT,
            CONSTRAINT fk_dna_id FOREIGN KEY (dna_id) REFERENCES "DNA_Extraction" (dna_id)
            );
            """
        )

        # PacBio_Library Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "PacBio_Library" (
            pacbio_library_tube_id TEXT PRIMARY KEY,
            dna_id TEXT,
            pacb_num INTEGER,
            pacb_status TEXT,
            library_method TEXT,
            library_date TEXT,
            library_id TEXT,
            dna_treatment TEXT,
            index_well TEXT,
            barcode TEXT,
            shear_femtol_id TEXT,
            shear_av_size INTEGER,
            seq_femto_id TEXT,
            seq_av_size INTEGER,
            library_conc REAL,
            comment TEXT,
            CONSTRAINT fk_dna_id FOREIGN KEY (dna_id) REFERENCES "DNA_Extraction" (dna_id)
            );
            """
        )

        # HiC_Library Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "HiC_Library" (
            hic_library_tube_id TEXT PRIMARY KEY,
            lysate_id TEXT,
            hic_num INTEGER,
            hic_status TEXT,
            library_method TEXT,
            library_date DATE,
            library_id TEXT,
            prox_ligation_conc REAL,
            purified_dna_total REAL,
            index_set TEXT,
            library_conc REAL,
            library_size INTEGER,
            hic_comments TEXT,
            CONSTRAINT fk_lysate_id FOREIGN KEY (lysate_id) REFERENCES "HiC_Lysate" (lysate_id)
            );
            """
        )

        # Sequencing Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "Sequencing" (
            sequencing_id TEXT PRIMARY KEY,
            og_num TEXT GENERATED ALWAYS AS (substring(sequencing_id FROM '([A-Z]{2}[0-9]+)')) STORED,
            rna_library_tube_id TEXT,
            illumina_library_tube_id TEXT,
            ont_library_tube_id TEXT,
            pacbio_library_tube_id TEXT,
            hic_library_tube_id TEXT,
            technology TEXT,
            instrument TEXT,
            run_date DATE,
            run_id TEXT,
            seq_date TEXT GENERATED ALWAYS AS (split_part(run_id, '_', 2)) STORED,
            cell_id TEXT,
            smrt_num INTEGER,
            seq_comments TEXT,
            seq_type TEXT,
            CONSTRAINT fk_rna_library FOREIGN KEY (rna_library_tube_id) REFERENCES "RNA_Library" (rna_library_tube_id),
            CONSTRAINT fk_illumina_library FOREIGN KEY (illumina_library_tube_id) REFERENCES "Illumina_Library" (illumina_library_tube_id),
            CONSTRAINT fk_ont_library FOREIGN KEY (ont_library_tube_id) REFERENCES "ONT_Library" (ont_library_tube_id),
            CONSTRAINT fk_pacbio_library FOREIGN KEY (pacbio_library_tube_id) REFERENCES "PacBio_Library" (pacbio_library_tube_id),
            CONSTRAINT fk_hic_library FOREIGN KEY (hic_library_tube_id) REFERENCES "HiC_Library" (hic_library_tube_id)
            );
            """
        )



    # Commit changes
        conn.commit()
        print("Database and tables created successfully.")

    except Exception as e:
        print(f"Error: {e}")

    finally:
        # Close connections
        if cursor:
            cursor.close()
        if conn:
            conn.close()

if __name__ == "__main__":
    create_postgresql_db()
