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

        # ✅ Print a success message if connected
        print("✅ Successfully connected to PostgreSQL!")

        # Species Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS "Species" (
                species TEXT PRIMARY KEY,
                class TEXT,
                "ordr" TEXT,
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
            extraction_date DATE,
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
            av_size TEXT,
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

            -- Foreign Key Constraints for Libraries
            CONSTRAINT fk_rna_library FOREIGN KEY (rna_library_tube_id) REFERENCES "RNA_Library" (rna_library_tube_id),
            CONSTRAINT fk_illumina_library FOREIGN KEY (illumina_library_tube_id) REFERENCES "Illumina_Library" (illumina_library_tube_id),
            CONSTRAINT fk_ont_library FOREIGN KEY (ont_library_tube_id) REFERENCES "ONT_Library" (ont_library_tube_id),
            CONSTRAINT fk_pacbio_library FOREIGN KEY (pacbio_library_tube_id) REFERENCES "PacBio_Library" (pacbio_library_tube_id),
            CONSTRAINT fk_hic_library FOREIGN KEY (hic_library_tube_id) REFERENCES "HiC_Library" (hic_library_tube_id)
            );
            """
        )

        # Mitogenome Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS mitogenome_data (
                og_id TEXT NOT NULL,  -- Foreign Key
                tech TEXT NOT NULL,
                seq_date TEXT NOT NULL,  -- Foreign Key
                code TEXT NOT NULL,
                stats TEXT,
                length INTEGER,
                length_emma INTEGER,
                seqlength_12s INTEGER,
                seqlength_16s INTEGER,
                seqlength_CO1 INTEGER,
                cds_no INTEGER,
                trna_no INTEGER,
                rrna_no INTEGER,
                status TEXT,
                genbank TEXT,
                rrna12s INTEGER,
                rrna16s INTEGER,
                atp6 INTEGER,
                atp8 INTEGER,
                cox1 INTEGER,
                cox2 INTEGER,
                cox3 INTEGER,
                cytb INTEGER,
                nad1 INTEGER,
                nad2 INTEGER,
                nad3 INTEGER,
                nad4 INTEGER,
                nad4l INTEGER,
                nad5 INTEGER,
                mad6 INTEGER,
                tRNA_Phe INTEGER,
                tRNA_Val INTEGER,
                tRNA_LeuUAG INTEGER,
                tRNA_LeuUAA INTEGER,
                tRNA_Ile INTEGER,
                tRNA_Met INTEGER,
                tRNA_Thr INTEGER,
                tRNA_Pro INTEGER,
                tRNA_Lys INTEGER,
                tRNA_Asp INTEGER,
                tRNA_Glu INTEGER,
                tRNA_SerGCU INTEGER,
                tRNA_SerUGA INTEGER,
                tRNA_Tyr INTEGER,
                tRNA_Cys INTEGER,
                tRNA_Trp INTEGER,
                tRNA_Ala INTEGER,
                tRNA_Asn INTEGER,
                tRNA_Gly INTEGER,
                tRNA_Arg INTEGER,
                tRNA_His INTEGER,
                tRNA_Gln INTEGER,

                -- Composite Primary Key
                PRIMARY KEY (og_id, tech, seq_date, code)
            );
            """
        )


        # LCA Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS lca (
                og_id TEXT NOT NULL,
                tech TEXT NOT NULL,
                seq_date TEXT NOT NULL,
                code TEXT NOT NULL,
                annotation TEXT NOT NULL,
                taxonomy TEXT,
                lca TEXT,
                percent_match DECIMAL(5,2),
                length INTEGER,
                lca_run_date TEXT,
                region TEXT NOT NULL,

                -- Composite Primary Key
                PRIMARY KEY (og_id, tech, seq_date, code, annotation, region),

                -- Foreign Key Constraints (Referencing mitogenome_data)
                CONSTRAINT fk_mitogenome FOREIGN KEY (og_id, tech, seq_date, code) 
                    REFERENCES mitogenome_data(og_id, tech, seq_date, code)
            );
            """
        )
    
    # Draft Genome Table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS draft_genomes (
                sample TEXT,
                mach TEXT,
                seq_date TEXT,
                initial TEXT,
                passed_filter_reads INTEGER,
                low_quality_reads INTEGER,
                too_many_n_reads INTEGER,
                too_short_reads INTEGER,
                too_long_reads INTEGER,
                raw_total_reads BIGINT,
                raw_total_bases BIGINT,
                raw_q20_bases BIGINT,
                raw_q30_bases BIGINT,
                raw_q20_rate DECIMAL(7,6),
                raw_q30_rate DECIMAL(7,6),
                raw_read1_mean_length INTEGER,
                raw_read2_mean_length INTEGER,
                raw_gc_content DECIMAL(7,6),
                total_reads BIGINT,
                total_bases BIGINT,
                q20_bases BIGINT,
                q30_bases BIGINT,
                q20_rate DECIMAL(7,6),
                q30_rate DECIMAL(7,6),
                read1_mean_length INTEGER,
                read2_mean_length INTEGER,
                gc_content DECIMAL(7,6),
                homozygosity DECIMAL(4,2),
                heterozygosity DECIMAL(4,2),
                genomesize BIGINT,
                repeatsize BIGINT,
                uniquesize BIGINT,
                modelfit DECIMAL(5,2),
                errorrate DECIMAL(5,2),
                num_contigs INTEGER,
                num_contigs_mitochondrion INTEGER,
                num_contigs_plastid INTEGER,
                num_contigs_prokarya INTEGER,
                bp_mitochondrion BIGINT,
                bp_plastid BIGINT,
                bp_prokarya BIGINT,
                complete DECIMAL(4,1),
                single_copy DECIMAL(4,1),
                multi_copy DECIMAL(4,1),
                fragmented DECIMAL(4,1),
                missing DECIMAL(4,1),
                n_markers INTEGER,
                domain TEXT,
                number_of_scaffolds INTEGER,
                number_of_contigs INTEGER,
                total_length BIGINT,
                percent_gaps DECIMAL(5,2),
                scaffold_n50 INTEGER,
                contigs_n50 INTEGER,
                unique_k_mers_assembly BIGINT,
                k_mers_total BIGINT,
                qv DECIMAL(6,4),
                error DECIMAL(12,11),
                k_mer_set TEXT,
                solid_k_mers BIGINT,
                total_k_mers BIGINT,
                completeness DECIMAL(7,4),
                depmethod TEXT,
                adjust TEXT,
                readbp BIGINT,
                mapadjust DECIMAL(7,6),
                scdepth DECIMAL(4,2),
                estgenomesize BIGINT,
                PRIMARY KEY (sample, seq_date)  -- Composite Primary Key
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
