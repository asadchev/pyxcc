from pyxcc.gto import Basis, wfn
import pyxcc.gto.basis as basis
import pyxcc.atom as atom

# def test_atom():
#   assert(atom.key('h_0') == 'H_0')
#   assert(atom.Z('h') == 1)
#   assert(atom.Z('h_0') == 1)
#   assert(atom.symbol('h_0') == 'H')
#   assert(atom.symbol(1) == 'H')

def test_basis():
  # basis.atom_basis(1, '6-31g')
  # basis.atom_basis('h', '6-31g')
  # basis.atom_basis((1,0,0,0), '6-31g')
  # basis.atom_basis((1,[0,0,0]), '6-31g')
  # basis.atom_basis((1,[0,0,0]), '6-31g')
  # basis.atom_basis(('h',[0,0,0]), '6-31g')
  # basis.atom_basis(('h',[0,0,0]), '6-31g')

  basis.load('6-31g')

  bs = Basis('6-31g')
  #assert(bs[1] == bs['h'] == bs['H'])

#   # b = basis.make_basis('6-31g', 1, 'o')
#   # print('basis', b)
#   # Basis([])
#   # Basis({
#   #   2 : [
#   #     (0, [[0,1]]), # s
#   #     (1, [[0,1]]), # p
#   #   ]
#   # })

# # def test_bse():
# #   import pyxcc.gto.bse as bse
# #   bse.load('6-31g')
