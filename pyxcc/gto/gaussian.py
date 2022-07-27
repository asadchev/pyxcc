import itertools, math, numpy
from pyxcc import cartesian

class Shell:
  def __init__(self, L, primitives, R=None):
    self._L = L
    self._primitives = [ (a,C) for a,C in primitives ]
    if R: self._R = numpy.array(R)
  @property
  def L(self): return self._L
  @property
  def nbf(self): return cartesian.nbf(self.L)
  @property
  def primitives(self): return self._primitives
  @property
  def R(self): return self._R

# def primitives(arg, *args):
#   return arg.primitives

class Tuple(tuple):

  class Primitive(tuple):
    @property
    def exp(self): return [ a for (a,C) in self ]
    @property
    def C(self): return [ C for (a,C) in self ]

  def __new__(cls, *args):
    return super(Tuple,cls).__new__(cls, args)

  @property
  def L(self):
    return [s.L for s in self]

  @property
  def primitives(self):
    primitives = ( s.primitives for s in self)
    return (Tuple.Primitive(ab) for ab in itertools.product(*primitives))

  @property
  def shape(self):
    return tuple(s.nbf for s in self)

  # @property
  # def exp(self):
  #   return sum(map(lambda s: s.primitives, self))

  @property
  def R(self):
    return tuple(s.R for s in self)

def braket(bra, ket):
  return (Tuple(*bra),Tuple(*ket))

def overlap_center(*args):
  #print(list(args))
  R = sum(q*numpy.array(R) for (q,R) in args)
  q = sum(q for (q,R) in args)
  return R/q

def exp(first, second=None):
  if not (first and second): return 1
  (a,A),(b,B) = first,second
  norm = numpy.linalg.norm
  return math.exp(-a*b*norm(A-B)/(a+b))

class Primitive2:
  def __init__(self, ab, AB, cd, CD):
    self.P = overlap_center(*zip(ab.exp,AB))
    self.Q = overlap_center(*zip(cd.exp,CD))
    p = sum(ab.exp)
    q = sum(cd.exp)
    self.p = p
    self.q = q
    self.pq = (
      (math.prod(ab.C)*exp(*zip(ab.exp,AB)))*
      (math.prod(cd.C)*exp(*zip(cd.exp,CD)))/
      (2*math.pi**(2.5))/(p*q*math.sqrt(p+q))
    )
  def __str__(self):
    return ("(p=%s, P=%s, q=%s, Q=%s, pq=%s)" % tuple(self))
  def __iter__(self):
    return iter((self.p, self.P, self.q, self.Q, self.pq))

def primitives(bra,ket):
  return (
    Primitive2(p, bra.R, q, ket.R)
    for p,q in itertools.product(bra.primitives, ket.primitives)
  )
