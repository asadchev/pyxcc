from pyxcc import Basis
import pyxcc.basis as basis

def test_basis():
  basis.load('6-31g', path="/home/andrey/github/pyxcc/pyxcc/data/basis")
  Basis('6-31g')
  Basis([])
  Basis({
    2 : [
      (0, [[0,1]]), # s
      (1, [[0,1]]), # p
    ]
  })
