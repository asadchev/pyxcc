from pyxcc.gto.gaussian import *

def test_gaussian():

  A = Shell(1, [(.05, 3.4)], R=[0,0,0])

  bra,ket = braket([A,A],[A,A])

  for p in primitives2(bra,ket):
    print (p)
