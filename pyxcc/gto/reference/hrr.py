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
  def hrr(A,B,R,G):
    from pyxcc.cartesian import orbitals,axis
    #print (G.shape)
    for t in range(0,B):
      As = orbitals(*range(A,A+B-t+1))
      At = orbitals(*range(A,A+B-t))
      Bs = orbitals(t)
      Bt = orbitals(t+1)
      H = HGP.hrr1(As,Bs,At,Bt,R,G)
      G = H
    return H


class Direct:

  @staticmethod
  def hrr(A,B,R,G):
    from pyxcc.cartesian import nbf,orbitals,index0
    from numpy import prod
    index = lambda t: index0(t) - index0([A,0,0])

    As = []
    for L in range(0,B+1):
      for t in orbitals(L):
        idx = [index(a) for a in orbitals(A+L) if a >= t]
        As.append(G[idx,0])

    #print (As)

    R = numpy.array(R)
    I = numpy.zeros([nbf(A),nbf(B)])

    R_ = R
    R = lambda i,j,k: prod(R_**(i,j,k))

    # if B == 1:
    #   I[:,0] = R(1,0,0)*As[0] + As[1]
    #   I[:,1] = R(0,1,0)*As[0] + As[2]
    #   I[:,2] = R(0,0,1)*As[0] + As[3]
    #   return I

    # if B == 2:
    #   I[:,0] = R(2,0,0)*As[0] + R(1,0,0)*As[1] + R(1,0,0)*As[1] + As[4]
    #   I[:,1] = R(1,1,0)*As[0] + R(1,0,0)*As[2] + R(0,1,0)*As[1] + As[5]
    #   I[:,2] = R(1,0,1)*As[0] + R(1,0,0)*As[3] + R(0,0,1)*As[1] + As[6]
    #   I[:,3] = R(0,2,0)*As[0] + R(0,1,0)*As[2] + R(0,1,0)*As[2] + As[7]
    #   I[:,4] = R(0,1,1)*As[0] + R(0,1,0)*As[3] + R(0,0,1)*As[2] + As[8]
    #   I[:,5] = R(0,0,2)*As[0] + R(0,0,1)*As[3] + R(0,0,1)*As[3] + As[9]
    #   return I

    for L in range(0,B+1):
      for t in orbitals(L):
        for x in orbitals(B-L):
          #print ("b=",x+t,(x+t).index,"=",'x',x,'t',t,'@',index0(t))
          M = (x+t).comb(t)
          I[:,(x+t).index] += M*R(*x)*As[index0(t)]

    return I
