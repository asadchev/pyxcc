import pyxcc.cartesian as cart
import pyxcc.hermitian as herm
import pyxcc.gto.reference.md as md

from pytest import approx
import numpy as np
import itertools

def test_md_reference():

  a,b = 3.14,137.0
  Xab = [ -2.71, 0.0, 6.63 ]
  A,B,CD = 3,2,4

  r1 = md.R1(A+B+CD)
  E = herm.E(A,B,a,b,Xab)

  for q in herm.orbitals(CD):
    PQ_q = r1[ [ herm.index(p+q) for p in herm.orbitals(A+B) ] ]
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        assert(
          np.dot(Eij, PQ_q) == approx(md.recursive(i,j,None,PQ_q,a,b,Xab))
        )


  PQ = r1.tensor(herm.orbitals(A+B), herm.orbitals(CD))
  for q in herm.orbitals(CD):
    PQ_q = PQ[:,herm.index(q)]
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        assert(
          np.dot(Eij, PQ_q) == approx(md.recursive(i,j,None,PQ_q,a,b,Xab))
        )

  # MD 1-Factorisation
  AR = r1.tensor(cart.orbitals(A), herm.orbitals(B+CD))
  for q in herm.orbitals(CD):
    PQ_q = PQ[:,herm.index(q)]
    s = 1/((2*(a+b))**(A+B))
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        T = AR[cart.index(i),herm.index(j+q)]
        nr = herm.nbf2(A+B-1) # sans highest L=A+B
        assert(
          np.dot(Eij[:nr], PQ_q[:nr]) + s*T ==
          approx(np.dot(Eij, PQ_q))
        )
