import itertools,numpy

def nbf(L):
  return ((L+1)*(L+2))//2

class Orbital(tuple):
  def __new__(cls, i, j, k):
    return super(Orbital,cls).__new__(Orbital, map(int, (i,j,k)))
  def __add__(self,other):
    return Orbital(*(a+b for a,b in zip(self,other)))
  def __sub__(self,other):
    return Orbital(*(a-b for a,b in zip(self,other)))
  @property
  def L(self): return sum(self)
  @property
  def array(self): return numpy.array(self)
  @property
  def nonzero(self): return tuple(i for i,v in enumerate(self) if v)

def axis(j):
  return Orbital(*(int(i==j) for i in range(3)))

def orbitals(*args):
  f = list()
  for L in args:
    for l in range(0,L+1):
      i = L-l
      for k in range(0,l+1):
        j = l-k
        f.append(Orbital(i,j,k))
  return tuple(f)
