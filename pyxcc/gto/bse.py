from . import basis

def load(basis_name, bse="http://basissetexchange.org"):
  import requests, os
  # This allows for overriding the URL via an environment variable
  # Feel free to just use the base_url below
  base_url = os.environ.get('BSE_API_URL', bse)
  headers = {}
  #params = {'elements': ['H','C','O']}
  r = requests.get(
    base_url + '/api/basis/%s/format/json'% (basis_name),
    headers=headers,
    #params=params
  )
  return basis.load_json(r.json())
