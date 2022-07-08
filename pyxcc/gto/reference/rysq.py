import numpy as np
import math
import pyxcc.gto.boys as boys

class Roots:

  # ROUTINE RETURNS AN N BY N TRIANGULAR MATRIX C SUCH THAT
  # C(TRANSPOSE)SC=I,  WHERE I IS AN N BY N IDENTITY MATRIX.
  @staticmethod
  def gram_schmidt(N,F):
    C = np.zeros([N,N])
    for j in range(N):
      Fj = F[j:j+j+1]
      v = np.zeros(j+1)
      v[j] = 1
      norm = 0
      for k in range(0,j):
        Ck = C[0:k+1,k]
        cf = np.dot(Ck, Fj[:k+1])
        v[0:k+1] -= cf*Ck
        norm += cf**2
      C[0:j+1,j] = v/math.sqrt(Fj[j] - norm)
    return C

  def __init__(self,N,X):

    K = max(2,N)

    F = np.array([ boys.reference(i,X) for i in range(2*K+1) ])
    #print(F)
    P = Roots.gram_schmidt(K+1,F)

    R = np.zeros([K,K])
    W = np.zeros([K,K])

    # W[0,0]=Fm[0]
    # R[0,0]=Fm[1]/Fm[0]
    # DUM = math.sqrt(P[1,2]**2 - 4*P[0,2]*P[2,2])
    # R[0,1] = 0.5*(-P[1,2]-DUM)/P[2,2]
    # R[1,1] = 0.5*(-P[1,2]+DUM)/P[2,2]

    for k in range(1,K+1):
      p = np.polynomial.Polynomial(P[0:k+1,k])
      R[0:k,k-1] = p.roots()

      for k in range(0,K):
        for i in range(0,k+1):
          w = 0
          r = R[i,k]
          for j in range(k+1):
            p = np.polynomial.Polynomial(P[:,j])
            w += p(r)**2
          W[i,k] = 1/w

      #print (W)

    self._t2 = np.array(R[:,-1])
    self._W = np.array(W[:,-1])

  @property
  def t2(self): return self._t2

  @property
  def W(self): return self._W

  @property
  def t(self):
    return np.sqrt(self.t2)

  @property
  def U(self): return self.t2/(1-self.t2)

def roots2_weights(N,X):
  r = Roots(N,X)
  return (r.t2, r.W)


from ... import cartesian
from .. import gaussian
from math import prod

class integrals2d:

  @staticmethod
  def x_m_n(bra, ket, p, P, q, Q, t2):
    N = len(t2)
    m = sum(bra.L)
    n = sum(ket.L)
    t2 = t2.reshape(N,1) # implicit broadcasting
    PQ = (P-Q)
    alpha = (p+q)
    G = np.zeros([m+1,n+1,N,3])
    G[0,0] = 1
    if m:
      PA = (P-bra[0].R)
      g0 = (PA - alpha/p*PQ*t2)
      gm = 1/(2*p)*(1 - alpha/p*t2)
      G[1,0] = g0*G[0,0]
      for i in range(1,m):
        G[i+1,0] = g0*G[i,0] + gm*i*G[i-1,0]
    if n:
      PC = (P-ket[0].R)
      g0 = (PC - alpha/q*PQ*t2)
      gn = 1/(2*q)*(1 - alpha/q*t2)
      gm = alpha/(2*p*q)
      G[:,1] = g0*G[:,0]
      # G[i,j+1] = g0*G[i,j] + j*gn*G[i,j-1] + i*gm*G[i-1,j]
      for j in range(1,n):
        G[0,j+1] = g0*G[0,j] + gn*j*G[i,j+1]
        for i in range(1,m+1):
          G[i,j+1] = g0*G[i,j] + gn*j*G[i,j-1] + gm*i*G[i-1,j]
    return np.moveaxis(G, -1, 0)

  @staticmethod
  def x_m_cd(I, *ket): return integrals2d.transfer(I,2,*ket)

  @staticmethod
  def x_ab_cd(I, *bra): return integrals2d.transfer(I,1,*bra)

  @staticmethod
  def transfer(G, index, first, second=None):
    if second is None: return G
    G = np.moveaxis(G,index,0)
    shape = list(G.shape)
    if not second.L:
      G = G.reshape(G.shape[0], 1, G.shape[2:])
      G = np.moveaxis(G, [0,1], [index,index+1])
      return G
    I = np.zeros((first.L+1, second.L+1) + G.shape[1:])
    R = (first.R - second.R)
    # implicit broadcast
    R = R.reshape([1,3] + [1 for i in range(G.ndim-2)])
    m = first.L+1
    i = slice(0,first.L+1)
    I[i,0] = G[i]
    for j in range(1,second.L+1):
      G[0:m+1-j] = R*G[0:m+1-j] + G[1:1+m+1-j]
      I[:,j] = G[i]
      #print (I.shape)
    I = np.moveaxis(I, [0,1], [index,index+1])
    return I


def eri(bra,ket):

  bra,ket = gaussian.braket(bra,ket)
  L = sum(bra.L+ket.L)

  V = np.zeros(bra.shape+ket.shape)
  #print (V.shape)

  V = V.reshape([ V.size ])

  for p,P,q,Q,pq in gaussian.primitives(bra,ket):
    PQ = P-Q
    X = (p+q)*np.linalg.norm(P-Q)
    t2,W = roots2_weights(L//2+1, X)
    I = integrals2d.x_m_n(bra, ket, p, P, q, Q, t2)
    #print(I.shape)
    I = integrals2d.x_m_cd(I, *ket)
    #print(I.shape)
    I = integrals2d.x_ab_cd(I, *bra)
    #print(I.shape)

    for i,orbitals in enumerate(cartesian.orbitals(bra.L+ket.L)):
      x,y,z = zip(*orbitals)
      Ix = I[0][x]
      Iy = I[1][y]
      Iz = I[2][z]
      #print('x',x,Ix,Ix.shape)
      v = 0
      for a,w in enumerate(W):
        v += w*Ix[a]*Iy[a]*Iz[a]
      V[i] += pq*v

  V = V.reshape(bra.shape+ket.shape)
  return V
