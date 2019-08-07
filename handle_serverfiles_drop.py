import sys
import requests

header = {'Authorization': f'Token {str(sys.argv[1])}'}
url = f'http://service.biolab.si/serverfiles_update/bioinformatics/{str(sys.argv[2])}'
response = requests.post(url, headers=header)

# Fail build if service.biolab.si unavailable
response.raise_for_status()
