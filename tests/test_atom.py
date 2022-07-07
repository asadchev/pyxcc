import pyxcc.atom as atom

def test_atom():

  assert(atom.key('h_0') == 'H_0')
  assert(atom.Z('h') == 1)
  assert(atom.Z('h_0') == 1)
  assert(atom.symbol('h_0') == 'H')
  assert(atom.symbol(1) == 'H')

  assert(atom.Z(100) == 100)
  assert(atom.key(100) == '100')
  assert(atom.symbol(100) == '100')
  
