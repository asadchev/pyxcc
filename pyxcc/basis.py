from collections import namedtuple
import os.path

def load_json(basis):
  if isinstance(basis, str):
    from json import loads as load
    basis = load(basis)
  elif hasattr(basis, "read"):
    from json import load as load
    basis = load(basis)
  else:
    pass
  elements = basis['elements']
  basis = {}
  for Z in map(int,elements):
    basis[Z] = []
    # print(Z)
    for f in (elements[str(Z)]['electron_shells']):
      angular_momentum = f['angular_momentum']
      exponents = list(map(float, f['exponents']))
      coefficients = [list(map(float,c)) for c in f['coefficients']]
      for i,L in enumerate(angular_momentum):
        primitives = list(zip(exponents, coefficients[i]))
      basis[Z].append((L, primitives))
  return basis


def load(basis, format="json", path=None):
  basis = basis + "." + format
  file = None
  if not path:
    import pyxcc.data.basis
    from importlib_resources import files
    file = files(pyxcc.data.basis).joinpath(basis).read_text()
  else:
    file = open(os.path.join(path, basis))
  return load_json(file)

def load_from_bse(basis_name, bse="http://basissetexchange.org"):
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
  return loads(r.json())

Primitive = namedtuple("Primitive", ["exp", "C"])    

class Shell:  
  def __init__(self, L=None, primitives=[], r=None, Z=None):
    self.L = int(L)
    self.primitives = tuple(
      Primitive(float(exp), float(C)) for exp,C in primitives
    )
    self.r = r
    self.Z = Z
  def __repr__(self):
    return "Shell(L=%i, primitives=%s)" % (self.L, self.primitives)

class Basis(dict):
  def __init__(self, basis, name=None):
    if isinstance(basis,str):
      name = basis
      basis = load(name.lower())
    if isinstance(basis,dict):
      basis = list(basis.items())
    self.name = name
    for (z,b) in basis:
      if z not in self: self[z] = []
      self[z].extend([Shell(L,C) for (L,C) in b])
  def __repr__(self):
    return "\n".join("%i : %r" % (k,b) for (k,b) in self.items())
