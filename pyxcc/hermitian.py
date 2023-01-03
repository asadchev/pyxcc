import itertools,math,numpy
from . import cartesian as cart

def nbf2(L):
  return cart.nbf(0,L)

def index(p):
  return cart.index(p,0)

def orbitals(L):
  return cart.orbitals(*range(0,L+1))

def E1(A,B,p,a,b,Xab):
  if (A == B == p == 0): return 1
  if (p < 0 or p > A + B): return 0
  if B:
    B = B-1
    X = +a*Xab/(a+b)
  else:
    assert(A)
    A = A-1
    X = -b*Xab/(a+b)
  return (
    0.5/(a+b)*E1(A,B,p-1,a,b,Xab) +
    X*E1(A,B,p,a,b,Xab) +
    (p+1)*E1(A,B,p+1,a,b,Xab)
  )

def E(A,B,a,b,Xab):
  E = numpy.zeros([cart.nbf(A),cart.nbf(A),nbf2(A+B)])
  for i,(ix,iy,iz) in enumerate(cart.orbitals(A)):
    for j,(jx,jy,jz) in enumerate(cart.orbitals(B)):
      for p,(px,py,pz) in enumerate(orbitals(A+B)):
        Ex = E1(ix,jx,px,a,b,Xab[0])
        Ey = E1(iy,jy,py,a,b,Xab[1])
        Ez = E1(iz,jz,pz,a,b,Xab[2])
        E[i,j,p] = Ex*Ey*Ez
  return E
