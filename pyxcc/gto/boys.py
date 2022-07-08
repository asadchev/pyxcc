import math

# computes a single value of \f$ F_m(T) \f$ using MacLaurin series to full precision of @c double
def reference(m,T):
  assert(m < 100)
  denom = (m + 0.5)
  term = math.exp(-T)/(2*denom)
  old_term = 0
  term_sum = term
  epsilon = 1e-16 # get_epsilon(T)
  epsilon_divided_10 = epsilon / 10
  while True:
    denom += 1
    old_term = term
    term = old_term * T / denom
    term_sum += term
    #rel_error = term / term_sum , hence iterate until rel_error = epsilon
    # however, must ensure that contributions are decreasing to ensure that omitted contributions are smaller than epsilon
    #print (term,term_sum,old_term)
    if (term > term_sum * epsilon_divided_10 or old_term < term): continue
    break
  return term_sum

class Boys:
  @staticmethod
  def reference(m,T):
    return reference(m,T)
