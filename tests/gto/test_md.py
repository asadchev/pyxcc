import pyxcc.cartesian as cart
import pyxcc.hermitian as herm
import pyxcc.gto.reference.md as md

from pytest import approx
import numpy as np
import itertools

def test_md_reference():

  a,b = 3.14,3.70
  Xab = [ -2.71, 0.0, 6.63 ]
  A,B,CD = 3,2,3

  r1 = md.R1(A+B+CD)
  PQ = r1.tensor(herm.orbitals(A+B), herm.orbitals(CD))
  E = herm.E(A,B,a,b,Xab)

  for q in herm.orbitals(CD):
    PQ_q = r1[ [ herm.index(p+q) for p in herm.orbitals(A+B) ] ]
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        assert(
          np.dot(Eij, PQ_q) == approx(md.recursive(i,j,None,PQ_q,a,b,Xab))
        )

  for q in herm.orbitals(CD):
    PQ_q = PQ[:,herm.index(q)]
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        assert(
          np.dot(Eij, PQ_q) == approx(md.recursive(i,j,None,PQ_q,a,b,Xab))
        )

  # MD 1-Factorisation
  #t1 = r1.tensor(cart.orbitals(A), herm.orbitals(B+CD))
  for q in herm.orbitals(CD):
    PQ_q = PQ[:,herm.index(q)]
    s = 0.5/(a+b)
    s = s**(A+B)
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        T = r1[herm.index(i+j+q)]
        nr = herm.nbf2(A+B-1) # sans highest L=A+B
        result = np.dot(Eij[:nr], PQ_q[:nr]) + s*T
        assert(result == approx(np.dot(Eij, PQ_q)))

  # MD 2-Factorisation
  for q in herm.orbitals(CD):
    s = 0.5/(a+b)
    PQ_q = r1[ [herm.index(i+q) for i in herm.orbitals(A+B) ] ]
    s0 = s**(A+B)
    s1 = -b/(a+b)*s**(A+B-1)
    s2 = +a/(a+b)*s**(A+B-1)
    for i in cart.orbitals(A):
      for j in cart.orbitals(B):
        Eij = E[cart.index(i),cart.index(j),:]
        T0,T1,T2 = 0,0,0
        # [i,j+q] -> [i,j,q]
        T0 = r1[herm.index(i+j+q)]
        # [i-1,j+q] -> [i,j,q]
        for k in [0,1,2]:
          if not i[k]: continue
          ik = cart.Orbital('xyz'[k])
          idx = herm.index((i+j+q)-ik)
          T1 += Xab[k]*i[k]*r1[idx]
        # [j-1,i+q] -> [i,j,q]
        for k in [0,1,2]:
          if not j[k]: continue
          jk = cart.Orbital('xyz'[k])
          idx = herm.index((i+j+q)-jk)
          T2 += Xab[k]*j[k]*r1[idx]
        n2 = herm.nbf2(A+B-2) # sans highest L=A+B-1,A+B
        result = np.dot(Eij[:n2], PQ_q[:n2]) + s0*T0 + s1*T1 + s2*T2
        u = np.dot(Eij, PQ_q)
        assert(
          np.dot(Eij[:n2], PQ_q[:n2]) + s0*T0 + s1*T1 + s2*T2 ==
          approx(np.dot(Eij, PQ_q))
        )
