import requests

from office365.runtime.auth.authentication_context import AuthenticationContext
from office365.sharepoint.client_context import ClientContext

# Define your SharePoint site URL
site_url = 'https://uniwa.sharepoint.com/teams/EXT-OceanOmicsLab'

# Authentication cookies (replace with actual values)
# These will likely need updating. To update you go to sharepoint on chrome browser then press F12 to get developer tools open.
# Then in developer tools navigate to application tab and copy the values over for the two cookies.
cookies = {
    "FedAuth": "77u/PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0idXRmLTgiPz48U1A+VjEzLDBoLmZ8bWVtYmVyc2hpcHwxMDAzMjAwMjVjODFhNGIyQGxpdmUuY29tLDAjLmZ8bWVtYmVyc2hpcHwwMDExMTIwMkB1d2EuZWR1LmF1LDEzMzgzNzExNjQxMDAwMDAwMCwxMzM4MTIyMTI1OTAwMDAwMDAsMTMzODM4MTA1NDg4NzAxOTIzLDEzMC45NS4xODkuMTg0LDIsMDU4OTRhZjAtY2IyOC00NmQ4LTg3MTYtNzRjZGI0NmUyMjI2LCwwMDIwM2ExOS04NjRlLThkMWYtOTM1NS0yYmY1NDUyODlmZWMsNjE1OTgwYTEtNzA4NS00MDAwLThmNmQtZDU3NjMyMzEyZWU4LDYxNTk4MGExLTcwODUtNDAwMC04ZjZkLWQ1NzYzMjMxMmVlOCwsMCwxMzM4MzcyNzc0ODgzODkzMDEsMTMzODM5ODMzNDg4Mzg5MzAxLCwsZXlKNGJYTmZZMk1pT2lKYlhDSkRVREZjSWwwaUxDSjRiWE5mYzNOdElqb2lNU0lzSW5CeVpXWmxjbkpsWkY5MWMyVnlibUZ0WlNJNklqQXdNVEV4TWpBeVFIVjNZUzVsWkhVdVlYVWlMQ0oxZEdraU9pSnFRMVZ2YUVGSFkwbHJMVXBGT0RsbVZqRnpRMEZCSW4wPSwyNjUwNDY3NzQzOTk5OTk5OTk5LDEzMzgzNzI0MTQ3MDAwMDAwMCw0YTkwYmRlMS0wMDEzLTQwMmItOTQ3ZC1mZTk5NTBiNDBhYmUsLCwsLCwxMTUyOTIxNTA0NjA2ODQ2OTc2LCwxOTU1ODAsdVhlaFFKUGxlVmpOQ2Jha1VoR0Q2SXlGUVFrLGpLUzdsbS93V3c2NzFZRGVJUWFGWTl1d09xemdleGhjVmMvWUFYaEpUa2tnQnRkcXdHeDV0Ym9NVGxSMm9qbVlwSkI4alN1RkhqRE9XLzJ0dHhuYnJ4WlVyVzRzWC9zelFHbWpsVDJlc0dKbTNyNGp5ZklPM0V5Mm5PNkhJMDhHTHlzaWJ2cDBPMEhPdlRkWkpYL3JpeGFaNjh6enc1RkZRVUtheEhiSmplRTlMS1hJYmZiOXdlWE9TVjFYQmJlbndFYzZCdHhWTjhaak1IZURsaXNGWndzTFVKUU96Qlk4eDhZTVhWZkFFNDlYQ2hmYkRlc1RGWENpNkc2L3FOR0xJdnRtZEpXTGM2VnBtcmV0YzFleDNhMmVIUUFGWCsrclM5NGwrRXJHN1JRcXFjbkNIUXB0VDlhY2hIM2NvL3pFUE5ub1JmWG9Rb2NTTitUbThMYXVoQT09PC9TUD4=",
    "rtFa": "Cy1qXI2pbpSxggEcZpzjqPEXzDGcDel1ZduZdErFN/MmMDU4OTRhZjAtY2IyOC00NmQ4LTg3MTYtNzRjZGI0NmUyMjI2IzEzMzgzNzI0MTQ4ODcwMTgzNiM2MTU5ODBhMS03MDg1LTQwMDAtOGY2ZC1kNTc2MzIzMTJlZTgjMDAxMTEyMDIlNDB1d2EuZWR1LmF1IzE5NTU4MCNiM0lRaDQ1a3lhZnh2ZkpUTlNQM1JzdW9PRk0jUDNtWnBmRU4yeV9rVjhJbG5lUzJXTTJBWkhZkkAerr6b6mNtOGvlnvPWViBdFkl2Sn+DICdYc+MsT+ccuLA6N7xlRwNCJ1Z4G7GVvw5fA6OJpe6wiZcGJmB/JrAt8tu1IXTSzMBb+HmV2pr6GmiYnZ03a4cWC+FsmUMrYz2LsOXLZ8fanXtU+hOxGjF39uoZ9Jtn6TLjM66sa+gbsKPBbE5SLncr+HJOKCN+olG4i0grg8rNypxHeieDJK/xB8QgIpcBJxC8Ur5oAZ30Xs0wimHt/ASlPw6/HrAfaxmwn9FNhfezb7DPetCj+ry760gMyf6IvUqAKuK3F2l7iHq+XlnIutVOpJuA72uNut0UVumQAy00gtXF5j6urdIAAAA=",
}

### Download species data spreadsheet code
# File URL (relative to the SharePoint site)
file_path = '/teams/EXT-OceanOmicsLab/Shared%20Documents/Database/Master%20Species%20List%20-%20Marine%20Vertebrates.xlsx'
file_url = f"{site_url}/_api/web/GetFileByServerRelativeUrl('{file_path}')/$value"

print("Final URL:", file_url)

# Make the request
response = requests.get(file_url, cookies=cookies)

from datetime import datetime

# Get the current date in YYYY-MM-DD format
current_date = datetime.now().strftime("%y%m%d")

# Construct the filename with the date
filename = f"Master_species_list_{current_date}.xlsx"

# Save the file with the new name
if response.status_code == 200:
    with open(filename, "wb") as file:
        file.write(response.content)
    print(f"File downloaded successfully as {filename}")
else:
    print(f"Failed to download file: {response.status_code}")
    print("Response content:", response.text)

