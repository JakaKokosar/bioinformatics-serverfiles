import sys
import requests

token = str(sys.argv[1])
run_id = str(sys.argv[2])
version = str(sys.argv[3])

resp = requests.get('https://api.github.com/repos/JakaKokosar/bioinformatics-serverfiles/actions/artifacts')
resp2 = requests.get(f'https://api.github.com/repos/JakaKokosar/bioinformatics-serverfiles/actions/runs/{run_id}/artifacts')

print(resp.text)
print()
print(resp2.text)

header = {'Authorization': f'Token {token}'}
url = f'https://service.biolab.si/serverfiles_update/bioinformatics/{run_id}?version={version}'
response = requests.get(url, headers=header)

print(url)
print(response.content)

# Fail build if service.biolab.si unavailable
response.raise_for_status()
