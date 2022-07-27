import pyxcc.cartesian

import numpy
import itertools

class HGP:

  @staticmethod
  def recurrence(t,R):
    u = numpy.array([0,0,0])
    r = 0.0
    i = next(iter(t.nonzero), None)
    if i is not None:
      u[i] = 1
      r = R[i]
    return (u,r)

  @staticmethod
  def hrr1(As,Bs,A,B,R,G):
    G_ = G
    G = lambda a,b: G_[As.index(a),Bs.index(b)]
    H = numpy.ndarray([len(A),len(B)])
    for ia,a in enumerate(A):
      for ib,b in enumerate(B):
        i,Ri = HGP.recurrence(b,R)
        H[ia,ib] = Ri*G(a,b-i) + G(a+i,b-i)
    return H

  @staticmethod
  def hrr(bra,ket,I):
    if bra and bra.second:
      I = hrr1(bra.first, bra.second, I)
    if ket and ket.second:
      I = numpy.moveaxis(I,[-1],0)
      I = hrr1(ket.first, ket.second, I)
      I = hrr(numpy.moveaxis(I,[0,1],-1))
    return I
