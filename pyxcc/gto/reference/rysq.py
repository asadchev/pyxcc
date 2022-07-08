import numpy as np
import math
import pyxcc.gto.boys as boys

class Q:

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

    K = N
    if K < 2: K=2

    F = np.array([ boys.reference(i,X) for i in range(2*K+1) ])
    #print(F)
    P = Q.gram_schmidt(K+1,F)

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
