# OceanOmics-Database
Postgresql database code, adapted from Adams sqlite code.

## Updating data
To update the data in the database to the most recent you need to first run the import_sharepoint.py, this will download the excel lab database and the TOLID spreadsheet with the days date to log when it was last downloaded.

Use the import_sharepoint_species_data.py script to download the master species list for updating the species data. Use the import_species_data.py to update the species table in the database.