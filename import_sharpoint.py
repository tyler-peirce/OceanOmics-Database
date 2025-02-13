import requests

from office365.runtime.auth.authentication_context import AuthenticationContext
from office365.sharepoint.client_context import ClientContext

# Define your SharePoint site URL
site_url = 'https://uniwa.sharepoint.com/teams/EXT-OceanOmicsLab'

# Authentication cookies (replace with actual values)
# These will likely need updating. To update you go to sharepoint on chrome browser then press F12 to get developer tools open.
# Then in developer tools navigate to application tab and copy the values over for the two cookies.
cookies = {
    "FedAuth": "77u/PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0idXRmLTgiPz48U1A+VjEzLDBoLmZ8bWVtYmVyc2hpcHwxMDAzMjAwMjVjODFhNGIyQGxpdmUuY29tLDAjLmZ8bWVtYmVyc2hpcHwwMDExMTIwMkB1d2EuZWR1LmF1LDEzMzgzODgxMjQ3MDAwMDAwMCwxMzM4MTIyMTI1OTAwMDAwMDAsMTMzODM5Njc2NDgwMTEyMzU1LDEzMC45NS4xODkuMTg0LDIsMDU4OTRhZjAtY2IyOC00NmQ4LTg3MTYtNzRjZGI0NmUyMjI2LCwwMDIwODgzOS0yYzg2LWE0MTYtMzQ2Yy0xOTcxY2NlZjAzYzksMzNlZjgwYTEtYzBjZC00MDAwLThmNmQtZDZjMjNlOGZlMTBkLDMzZWY4MGExLWMwY2QtNDAwMC04ZjZkLWQ2YzIzZThmZTEwZCwsMCwxMzM4Mzg4NDg0Nzk5NTU4NjAsMTMzODQxNDA0NDc5OTU1ODYwLCwsZXlKNGJYTmZZMk1pT2lKYlhDSkRVREZjSWwwaUxDSjRiWE5mYzNOdElqb2lNU0lzSW5CeVpXWmxjbkpsWkY5MWMyVnlibUZ0WlNJNklqQXdNVEV4TWpBeVFIVjNZUzVsWkhVdVlYVWlMQ0oxZEdraU9pSmxMVXhWV1VkMVdGRlZaVzlVTm1Ka1VIRm5ha0ZCSW4wPSwyNjUwNDY3NzQzOTk5OTk5OTk5LDEzMzgzODgxMjQ3MDAwMDAwMCw0YTkwYmRlMS0wMDEzLTQwMmItOTQ3ZC1mZTk5NTBiNDBhYmUsLCwsLCwxMTUyOTIxNTA0NjA2ODQ2OTc2LCwxOTU1ODAsdVhlaFFKUGxlVmpOQ2Jha1VoR0Q2SXlGUVFrLFpnb0g4alhQeTJuVktMeXlhWWNkak9tZjhNVW8zR1V5dFdSdjlVc213WmZRckxaWExOUGFwdE9PZWErclBXSkFaRGk3K3hCZUREdWlOQ3F6SkF4N0I0QXcxQVRWaDJXKzh6VXN2ZWRGVm5nR2hjQmp1cmVrVFk5ZUdQOUViNUFpYmdzemxGa3VCSjUxN3ZsUkdQK3hFNjhySVgvRDVVYWF5RDFKSWJhY3JiUGZMc2VsVHNMNTg0Q3ZtZUVGdUVWTVlSVmY0Z1B2QUxqUVE5MGJTVkJFVXdOc3c1Y3NvbDFhQm13UGZKbFoxZVREYUtWMHhJSTlGWHlSQmkwS1EwZ1R5S25RZSswWThJa3J4NklkTEVuTXNYSFJhMitaZkxwc3I4eElGTElCdFpaclNYTGdjN2RjL0dZZG1JYmhoZ3oyV3B3Mm5vKy9WYloxUGQ2V1hxUnNMUT09PC9TUD4=",
    "rtFa": "nLwlbYMsJ1OWQktF6LF1quA75JR6RH+S0ypWu98A+YwmMDU4OTRhZjAtY2IyOC00NmQ4LTg3MTYtNzRjZGI0NmUyMjI2IzEzMzgzODgxMjQ4MDI2ODUwNiMzM2VmODBhMS1jMGNkLTQwMDAtOGY2ZC1kNmMyM2U4ZmUxMGQjMDAxMTEyMDIlNDB1d2EuZWR1LmF1IzE5NTU4MCNiM0lRaDQ1a3lhZnh2ZkpUTlNQM1JzdW9PRk0jUDNtWnBmRU4yeV9rVjhJbG5lUzJXTTJBWkhZoh6LlFtHwku8/83F8GHe5D8vHPIe25HO1NpUTRR6AIbsP22N6I0NM4ENe8V2wPBGZbjCL7DcMx0W6zAtUjlgb0MCf8N9Glq7cxWA/8fD10XQ12pZZ8Bl8K1vzzc9Lbdn8aZzPSaeufryNW5LpnzR/KUM64BqX6aoauwD+lVVUbb/4rCz7ztyx4o5RTsT+HNXUYIg53IuvWTp07R7cilHtEjr7tIh5cW4B9nH4/g5T2V9wtKPANr8UZxf/z8jpkX3wAgvlz+8AYEwAsNhCZ/F/eV/XsTOZVv0yGbrvBKvUr0eGb6t/QCTIEz4guMRjM9iZG7/sxjFkRpcWX5N0kP49NIAAAA=",
}

### Download Lab Database exel code
# File URL (relative to the SharePoint site)
file_path = '/teams/EXT-OceanOmicsLab/Shared%20Documents/Lab%20-%20Sequencing%20Runs/OceanGenomes%20Database%20v5.Feb24.xlsx'
file_url = f"{site_url}/_api/web/GetFileByServerRelativeUrl('{file_path}')/$value"

print("Final URL:", file_url)

# Make the request
response = requests.get(file_url, cookies=cookies)

from datetime import datetime

# Get the current date in YYYY-MM-DD format
current_date = datetime.now().strftime("%y%m%d")

# Construct the filename with the date
filename = f"OceanGenomes_database_{current_date}.xlsx"

# Save the file with the new name
if response.status_code == 200:
    with open(filename, "wb") as file:
        file.write(response.content)
    print(f"File downloaded successfully as {filename}")
else:
    print(f"Failed to download file: {response.status_code}")
    print("Response content:", response.text)

## Download Toll ID spreadsheet
# File URL (relative to the SharePoint site)
file_path2 = '/teams/EXT-OceanOmicsLab/Shared%20Documents/Lab%20-%20Sequencing%20Runs/TOLID.xlsx'
file_url2 = f"{site_url}/_api/web/GetFileByServerRelativeUrl('{file_path2}')/$value"

print("Final URL:", file_url2)

# Make the request
response = requests.get(file_url2, cookies=cookies)

from datetime import datetime

# Get the current date in YYYY-MM-DD format
current_date = datetime.now().strftime("%y%m%d")

# Construct the filename with the date
filename2 = f"TOLID_{current_date}.xlsx"

# Save the file with the new name
if response.status_code == 200:
    with open(filename2, "wb") as file:
        file.write(response.content)
    print(f"File downloaded successfully as {filename2}")
else:
    print(f"Failed to download file: {response.status_code}")
    print("Response content:", response.text)