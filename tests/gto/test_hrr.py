import pyxcc.gto.reference.hrr as hrr
from pyxcc.cartesian import orbitals,nbf,axis

from pytest import approx
import numpy as np

def test_hrr_hgp():

  HGP = hrr.HGP

  A = 4
  B = 4

  R = (3.14, -2.71, 0.5)

  G = np.random.rand(
    sum(nbf(L) for L in range(A,A+B+1)),
    1
  )

  for t in range(0,B):
    As = orbitals(*range(A,A+B-t+1))
    At = orbitals(*range(A,A+B-t))
    Bs = orbitals(t)
    Bt = orbitals(t+1)
    H = HGP.hrr1(As,Bs,At,Bt,R,G)
    for k in (0,1,2):
      for ia,a in enumerate(At):
        for ib,b in enumerate(Bs):
          ia_plus_1 = As.index(a + axis(k))
          ib_plus_1 = Bt.index(b + axis(k))
          assert(H[ia,ib_plus_1] == approx(R[k]*G[ia,ib] + G[ia_plus_1,ib]))
    G = H
