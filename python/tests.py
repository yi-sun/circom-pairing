import time
from curve_field_elements import field_modulus, FQ, FQ2, FQ12
from curve import double, add, eq, multiply, is_on_curve, neg, b2, b12, curve_order, G1, G2, G12, is_inf
from pairing import pairing

test_field = False
if test_field:
  print('Starting bn128 tests')

  assert FQ(2) * FQ(2) == FQ(4)
  assert FQ(2) / FQ(7) + FQ(9) / FQ(7) == FQ(11) / FQ(7)
  assert FQ(2) * FQ(7) + FQ(9) * FQ(7) == FQ(11) * FQ(7)
  assert FQ(9) ** field_modulus == FQ(9)
  print('FQ works fine')

  x = FQ2([1, 0])
  f = FQ2([1, 2])
  fpx = FQ2([2, 2])
  one = FQ2.one()
  assert x + f == fpx
  assert f / f == one
  assert one / f + x / f == (one + x) / f
  assert one * f + x * f == (one + x) * f
  assert x ** (field_modulus ** 2 - 1) == one
  print('FQ2 works fine')

  x = FQ12([1] + [0] * 11)
  f = FQ12([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
  fpx = FQ12([2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
  one = FQ12.one()
  assert x + f == fpx
  assert f / f == one
  assert one / f + x / f == (one + x) / f
  assert one * f + x * f == (one + x) * f
  # This check takes too long
  # assert x ** (field_modulus ** 12 - 1) == one
  print('FQ12 works fine')

  # G1, G2, G12, b, b2, b12, is_inf, is_on_curve, eq, add, double, curve_order, multiply = \
  #   G1, G2, G12, b, b2, b12, is_inf, is_on_curve, eq, add, double, curve_order, multiply

  assert eq(add(add(double(G1), G1), G1), double(double(G1)))
  assert not eq(double(G1), G1)
  assert eq(add(multiply(G1, 9), multiply(G1, 5)), add(multiply(G1, 12), multiply(G1, 2)))
  assert is_inf(multiply(G1, curve_order))
  print('G1 works fine')

  assert eq(add(add(double(G2), G2), G2), double(double(G2)))
  assert not eq(double(G2), G2)
  assert eq(add(multiply(G2, 9), multiply(G2, 5)), add(multiply(G2, 12), multiply(G2, 2)))
  assert is_inf(multiply(G2, curve_order))
  assert not is_inf(multiply(G2, 2 * field_modulus - curve_order))
  assert is_on_curve(multiply(G2, 9), b2)
  print('G2 works fine')

  assert eq(add(add(double(G12), G12), G12), double(double(G12)))
  assert not eq(double(G12), G12)
  assert eq(add(multiply(G12, 9), multiply(G12, 5)), add(multiply(G12, 12), multiply(G12, 2)))
  assert is_on_curve(multiply(G12, 9), b12)
  assert is_inf(multiply(G12, curve_order))
  print('G12 works fine')

print('Starting pairing tests')
a = time.time()
p1 = pairing(G2, G1)
pn1 = pairing(G2, neg(G1))
assert p1 * pn1 == FQ12.one()
print('Pairing check against negative in G1 passed')
np1 = pairing(neg(G2), G1)
assert p1 * np1 == FQ12.one()
assert pn1 == np1
print('Pairing check against negative in G2 passed')
assert p1 ** curve_order == FQ12.one()
print('Pairing output has correct order')
p2 = pairing(G2, multiply(G1, 2))
assert p1 * p1 == p2
print('Pairing bilinearity in G1 passed')
assert p1 != p2 and p1 != np1 and p2 != np1
print('Pairing is non-degenerate')
po2 = pairing(multiply(G2, 2), G1)
assert p1 * p1 == po2
print('Pairing bilinearity in G2 passed')
p3 = pairing(multiply(G2, 27), multiply(G1, 37))
po3 = pairing(G2, multiply(G1, 999))
assert p3 == po3
print('Composite check passed')
print('Total time for pairings: %.3f' % (time.time() - a))