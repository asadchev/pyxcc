import pyxcc.gto.boys as boys

def test_boys():
  assert boys.Boys.reference(0,0) == 1
