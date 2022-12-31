import pyxcc.cartesian as cart
import pyxcc.hermitian as herm
import numpy,itertools

def axis(A):
  return cart.imax(A)

class Circulant(numpy.ndarray):
  def __new__(cls, *args, **kwargs):
    return super().__new__(cls, *args, **kwargs)
  def tensor(self,*args):
    idx = [ [ (i,v) for i,v in enumerate(k) ] for k in args ]
    v = numpy.ndarray([len(i) for i in idx])
    for ip in itertools.product(*idx):
      i = [ i[0] for i in ip ]
      p = [ p[1] for p in ip ]
      k = sum(p, start=cart.Orbital([0,0,0]))
      v[tuple(i)] = self[herm.index(k)]
    return v


def R1(L):
  return numpy.random.rand(herm.nbf2(L)).view(Circulant)

def recursive(A,B,p,R,a,b,Xab):
  if p is None: p = cart.Orbital([0,0,0])
  if A.L == B.L == 0:
    if min(p) < 0: return 0
    return R[herm.index(p)]
  if B.L:
    i = axis(B)
    t = cart.Orbital("xyz"[i])
    xi = +a*Xab[i]/(a+b)
    B = B-t
  else:
    assert(A.L)
    i = axis(A)
    t = cart.Orbital("xyz"[i])
    xi = -b*Xab[i]/(a+b)
    A = A-t
  pi = p[i]
  return (
    pi*recursive(A,B,p-t,R,a,b,Xab) +
    xi*recursive(A,B,p,  R,a,b,Xab) +
    0.5/(a+b)*recursive(A,B,p+t,R,a,b,Xab)
  )
