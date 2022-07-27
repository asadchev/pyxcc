import pyxcc.cartesian as cartesian

def test_cartesian():

  assert (cartesian.orbitals(0) == [(0,0,0)])

  p = [
     (1,0,0),
     (0,1,0),
     (0,0,1),
   ]

  d = [
     (2,0,0),
     (1,1,0),
     (1,0,1),
     (0,2,0),
     (0,1,1),
     (0,0,2)
   ]

  assert (cartesian.orbitals(1) == p)
  assert (cartesian.orbitals(2) == d)
  assert (cartesian.orbitals(3))
