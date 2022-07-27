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

  H = G
  for t in range(0,B):
    As = orbitals(*range(A,A+B-t+1))
    At = orbitals(*range(A,A+B-t))
    Bs = orbitals(t)
    Bt = orbitals(t+1)
    I = HGP.hrr1(As,Bs,At,Bt,R,H)
    for k in (0,1,2):
      for ia,a in enumerate(At):
        for ib,b in enumerate(Bs):
          ia_plus_1 = As.index(a + axis(k))
          ib_plus_1 = Bt.index(b + axis(k))
          assert(I[ia,ib_plus_1] == approx(R[k]*H[ia,ib] + H[ia_plus_1,ib]))
    H = I

  assert(np.allclose(H, HGP.hrr(A,B,R,G)))


def test_hrr_direct():

  print()

  Direct = hrr.Direct
  HGP = hrr.HGP

  A = 4
  B = 4

  R = (3.14, -2.71, 0.5)

  G = np.random.rand(
    sum(nbf(L) for L in range(A,A+B+1)),
    1
  )

  H = Direct.hrr(A,B,R,G)
  Ref = HGP.hrr(A,B,R,G)

  # print(H)
  # print(Ref)

  assert(np.allclose(Direct.hrr(A,B,R,G), Ref))
