from . import Basis, basis

class Wavefunction:
  def __init__(self, basis, charge=0, spin=0):
    if not isinstance(basis, Basis):
      basis = Basis(basis)
    self._basis = basis
      
def wfn(basis, atoms, charge=0, spin=0):
  if isinstance(basis, str):
    basis = gto.basis.load(basis)
  assert(isinstance(dict,basis))
  for a,b in basis.items():
    if isinstance(str,b):
      basis[a] = gto.basis.load(b).get(a)
  if isinstance(str,atoms):
    atoms = make_atoms(atoms)
  return Wavefunction(
    basis.make_basis(basis, atoms),
    charge,
    spin
  )
    
