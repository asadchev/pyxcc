from collections import namedtuple
import os.path
import pyxcc.atom as atom

def parse(basis, format="json"):
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


def load(name, file=None, format="json", key=None):
  data = None
  if not file:
    import pyxcc.data.basis
    from importlib_resources import files
    file = files(pyxcc.data.basis).joinpath("%s.json" % name.lower())
    data = file.read_text()
  else:
    data = open(path).read()
  return parse(data,format)

def get(basis,key):
  if not isinstance(basis,BasisSet): basis = BasisSet(basis)
  return basis.get(key)

def atom_basis(basis, atom, *args):
  if len(args) == 1:
    [x,y,z] = args[0]
  else:
    x,y,z = args
  return (atom,(x,y,z))

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

def basis_set(basis):
  if isinstance(basis,str):
    basis = load(basis)
  if isinstance(basis,dict):
    basis = basis.items()
  basis_set = dict()
  for (k,v) in basis:
    k = atom.key(k)
    Z = atom.Z(k)
    basis_set[k] = [ Shell(L,p) for (L,p) in v ]
  return basis_set

class BasisSet:
  ## k is Atom
  def __getitem__(self,a):
    a = Atom(a)
    basis = self.get(a.name)
    (name,symbol,Z) = (a.name, a.symbol, a.Z)
    if not basis and name: basis = self.get(name)
    if not basis and symbol: basis = self.get(symbol)
    if not basis and Z: basis = self.get(Z)
    #assert (basis o)
    return (Atom,basis)


class Basis(list):
  def __init__(self, basis, *args, name=None, pure=True):
    if isinstance(basis,str):
      basis_set = load(basis)
    # for (a,r) in args:
    #   self.extend(Shell(L,p,r=r,Z=Z) for (L,p,r,Z) in args)
    #   for s in get(basis,a):
    # )

def make_basis(basis, *args, pure=True):
  name = None
  if isinstance(basis,str):
    name = basis
    basis = BasisSet(basis)
  if isinstance(basis,dict):
    basis = BasisSet(basis)
  if not isinstance(basis,BasisSet):
    assert(not args)
    return Basis(basis, name=name, pure=pure)
  return Basis([ basis[a] for a in args ], name=name, pure=pure)
