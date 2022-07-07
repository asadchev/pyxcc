atoms = (
  None,
  "H",
  "He",
  "Li",
  "Be",
  "B",
  "C",
  "N",
  "O"
)

def symbol(a):
  assert(isinstance(a, (int,str)))
  if isinstance(a,int):
    try: return atoms[a]
    except: return str(a)
  a = a.capitalize()
  for i in range(len(a)):
    if a[i].isalpha(): continue
    return a[0:i]
  return a

def Z(a):
  assert(isinstance(a, (int,str)))
  if isinstance(a,int): return a
  try: return atoms.index(symbol(a))
  except: return None

def key(a):
  assert(isinstance(a, (int,str)))
  if isinstance(a,int): return symbol(a)
  return a.capitalize()

# def parse(args):
#   return [ a for a in args ]
  
