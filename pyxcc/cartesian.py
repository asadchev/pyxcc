import itertools,math,numpy

def triangular(L):
  return (L*(L+1))//2

def tetrahedral(L):
  return (L*(L+1)*(L+2))//6

def nbf(L, M=None):
  if M is None:
    return triangular(L+1)
  return (tetrahedral(M+1) - tetrahedral(L))

class Orbital(tuple):

  def __new__(cls, ijk):
    if isinstance(ijk,str):
      ijk = [ ijk.lower().count(x) for x in 'xyz' ]
    i,j,k = map(int, ijk)
    return super(Orbital,cls).__new__(Orbital, map(int, (i,j,k)))

  def __add__(self,other):
    return Orbital((a+b for a,b in zip(self,other)))

  def __sub__(self,other):
    return Orbital((a-b for a,b in zip(self,other)))

  def __lt__(self,other):
    if self.L != other.L: return self.L < other.L
    for a,b in zip(self,other):
      if a != b: return (a > b)
    return False

  def comb(self,k):
    return [math.comb(a,b) for a,b in zip(self,k)]

  def __str__(self):
    return "(%s)" % "".join(map(str,self))

  def __getitem__(self,k):
    if isinstance(k,str):
      k = 'xyz'.index(k)
    return super().__getitem__(k)

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
        f.append(Orbital([i,j,k]))
  return tuple(f)

def index(A,L=None):
  A = Orbital(A)
  idx = A.index
  if L is None: return idx
  return idx + (tetrahedral(A.L) - tetrahedral(L))

def index0(A):
  if isinstance(A,int): return tetrahedral(A+1)
  i,j,k = A
  return index0(i+j+k) + index(A)

def imax(A):
  return tuple.index(A,max(A))
