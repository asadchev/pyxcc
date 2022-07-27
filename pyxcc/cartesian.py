import itertools,math,numpy

def nbf(L):
  return ((L+1)*(L+2))//2

class Orbital(tuple):
  def __new__(cls, i, j, k):
    return super(Orbital,cls).__new__(Orbital, map(int, (i,j,k)))
  def __add__(self,other):
    return Orbital(*(a+b for a,b in zip(self,other)))
  def __sub__(self,other):
    return Orbital(*(a-b for a,b in zip(self,other)))
  def __ge__(self,other):
    return all(a>=b for a,b in zip(self,other))
  def comb(self,k):
    return math.prod(math.comb(a,b) for a,b in zip(self,k))
  @property
  def L(self): return sum(self)
  @property
  def array(self): return numpy.array(self)
  @property
  def nonzero(self): return tuple(i for i,v in enumerate(self) if v)
  @property
  def index(self):
    i,j,k = self
    jk = j+k;
    return (jk*(jk+1))//2 + k

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

def index(A):
  return Orbital(*A).index

def index0(A):
  i,j,k = A
  L = i+j+k
  first = (L*(L+1)*(L+2))//6
  return first + index(A)
