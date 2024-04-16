import polynomial
import rational

print(
    polynomial.polynomial_factorization(
        polynomial.Polynomial([rational.Rational(x) for x in input("輸入多項式降冪排列的係數：").split()])))
