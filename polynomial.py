from enum import Enum

import util
from rational import Rational


class Monomial:
    """ 任意整係數單項式 """

    def __init__(self, coefficient: int | float | Rational, degree: int = 0):
        """
        任意整係數單項式

        :param coefficient: 係數
        :param degree: 次數
        """
        self.coefficient = Rational(coefficient)
        self.degree = degree

    def get_coefficient(self) -> Rational:
        """ 回傳單項式的係數 """
        return self.coefficient

    def get_degree(self) -> int:
        """ 回傳單項式的次數 """
        return self.degree

    def __str__(self) -> str:
        _str = str(self.coefficient) if self.coefficient != 1 or self.degree == 0 else ""
        if self.degree != 0:
            _str += "x"
            if self.degree != 1:
                _str += f"^{self.degree}"
        return _str

    def __repr__(self) -> str:
        return self.__str__()

    def __mul__(self, other):
        if not isinstance(other, Monomial):
            other = Monomial(other)
        return Monomial(self.coefficient * other.coefficient, self.degree + other.degree)

    def __imul__(self, other):
        return self.__mul__(other)

    def __eq__(self, other):
        return self.coefficient == other.coefficient and self.degree == other.degree


class ArrangementEnum(Enum):
    """ 多項式的排列方式 """
    DESCENDING = 0
    """ 降冪排列 """
    ASCENDING = 1
    """ 升冪排列 """


class Polynomial:
    """ 任意整係數多項式 """

    def __init__(self, coefficients: tuple[int | float | Rational, ...] | list[
        int | float | Rational] | int | float | Rational | None = None, monomial: Monomial = None,
                 arrangement: ArrangementEnum = ArrangementEnum.DESCENDING):
        """
        任意整係數多項式

        預設為降冪排列
        """
        if coefficients is None and monomial is None:
            raise ValueError("At least one parameter is required.")
        if isinstance(coefficients, (tuple, list)):
            if arrangement == ArrangementEnum.DESCENDING:
                coefficients = coefficients[::-1]
            while coefficients[-1] == 0:
                coefficients = coefficients[:-1]
            coefficients = list(map(Rational, coefficients))
        elif monomial is not None:
            coefficients = [Rational(0)] * monomial.get_degree()
            coefficients.append(monomial.get_coefficient())
        else:
            coefficients = [Rational(coefficients)]
        self.coefficients = coefficients

    def __len__(self) -> int:
        return len(self.coefficients)

    def __str__(self) -> str:
        _str = "("
        for i in range(len(self))[::-1]:
            if self.coefficients[i] != 0:
                if self.coefficients[i] > 0 and i != self.get_degree():
                    _str += "+"
                _str += str(self.get_monomial_by_degree(i))
        _str += ")"
        return _str

    def __repr__(self) -> str:
        return self.__str__()

    def get_coefficients(self) -> tuple[Rational, ...]:
        """ 回傳多項式的所有係數 """
        return tuple(self.coefficients)

    def get_degree(self) -> int:
        """ 回傳多項式的最高次數 """
        return len(self) - 1

    def lowest_degree_coefficient(self) -> Rational:
        """ 回傳多項式的常數項 """
        return self.coefficients[0]

    def highest_degree_coefficient(self) -> Rational:
        """ 回傳多項式的最高次項係數 """
        return self.coefficients[-1]

    def get_monomial_by_degree(self, degree: int) -> Monomial:
        """ 回傳多項式的 int degree 次項 """
        return Monomial(self.coefficients[degree], degree)

    def test(self, value: Rational | int | float) -> bool:
        """ 測試 int value 代入後是否為 0 """
        result = Rational(0)
        for i, coefficient in enumerate(self.coefficients):
            result += coefficient * (value ** i)
        return result == 0

    def __floordiv__(self, other):
        """
        進行綜合除法，回傳商式
        :return: 商式
        """
        if not isinstance(other, (Polynomial, Monomial)):
            raise TypeError("Unsupported type.")
        if isinstance(other, Monomial):
            other = Polynomial(monomial=other)
        if self.get_degree() < other.get_degree():
            return Polynomial(0)
        dividend_coefficients = self.coefficients[::-1]
        leading_coefficient = other.highest_degree_coefficient()
        divisor_coefficients = [-Rational(coefficient, leading_coefficient) for coefficient in
                                other.coefficients[-2::-1]]
        coefficients_difference = len(dividend_coefficients) - len(divisor_coefficients)
        table: list[list[Rational]] = [[Rational(0)] * coefficients_difference for _ in
                                       range(len(divisor_coefficients))]
        quotient_coefficients: list[Rational] = [dividend_coefficients[i] for i in range(coefficients_difference)]
        for i in range(coefficients_difference):
            for j in range(len(divisor_coefficients)):
                if i - j - 1 >= 0:
                    quotient_coefficients[i] += table[j][i - j - 1]
            for j, divisor_coefficient in enumerate(divisor_coefficients):
                table[j][i] = divisor_coefficient * quotient_coefficients[i]
        return Polynomial(
            [quotient_coefficient / leading_coefficient for quotient_coefficient in quotient_coefficients])

    def __ifloordiv__(self, other):
        return self.__floordiv__(other)

    def is_monomial(self) -> bool:
        """ 判斷是否為單項式 """
        for coefficient in self.coefficients[1:]:
            if coefficient != 0:
                return False
        return True

    def to_monomial(self) -> Monomial:
        """ 轉換為單項式 """
        return Monomial(self.highest_degree_coefficient(), self.get_degree())

    def __eq__(self, other):
        return self.coefficients == other.coefficients


class Quadratic(Polynomial):
    """ 二次整係數多項式 """

    class DiscriminantEnum(Enum):
        """ 判別式的結果 """
        TWO_SAME_ROOTS = 2
        """ 重實根 """
        TWO_DIFFERENT_ROOTS = 1
        """ 兩相異實根 """
        NO_REAL_ROOT = 0
        """ 無實根 """

    def __init__(self, a: int | float | Rational, b: int | float | Rational, c: int | float | Rational,
                 arrangement: ArrangementEnum = ArrangementEnum.DESCENDING):
        super().__init__((a, b, c), arrangement=arrangement)

    def get_a(self) -> Rational:
        """ 回傳二次多項式的二次項係數 """
        return self.coefficients[2]

    def get_b(self) -> Rational:
        """ 回傳二次多項式的一次項係數 """
        return self.coefficients[1]

    def get_c(self) -> Rational:
        """ 回傳二次多項式的常數項 """
        return self.coefficients[0]

    def discriminant(self) -> DiscriminantEnum:
        """
        回傳判別式的結果

        :return: 判別式的結果
        """
        result = self.get_b() ** 2 - Rational(4) * self.get_a() * self.get_c()
        if result < 0:
            return self.DiscriminantEnum.NO_REAL_ROOT
        if result == 0:
            return self.DiscriminantEnum.TWO_SAME_ROOTS
        if result > 0:
            return self.DiscriminantEnum.TWO_DIFFERENT_ROOTS


class Polynomials:
    def __init__(self, polynomials: any = (),
                 coefficient: Rational | Monomial | int | float = Rational(1)):
        if not isinstance(coefficient, Monomial):
            coefficient = Monomial(coefficient)
        if isinstance(polynomials, (list, tuple)):
            coefficient = coefficient
            polynomials = list(polynomials)
        elif isinstance(polynomials, Polynomial):
            coefficient = coefficient
            polynomials = [polynomials]
        elif isinstance(polynomials, Monomial):
            coefficient = coefficient * Monomial(polynomials.coefficient)
            polynomials = [polynomials]
        elif isinstance(polynomials, Polynomials):
            coefficient = polynomials.get_coefficient() * coefficient
            polynomials = list(polynomials.get_polynomials())
        else:
            raise ValueError("Unsupported type.")
        polynomials.sort(key=lambda polynomial: (polynomial.get_degree(), polynomial.coefficients))
        self.polynomials = polynomials
        self.coefficient = coefficient

    def __str__(self) -> str:
        _str = str(self.coefficient) if self.coefficient != Monomial(1) else ""
        count = 1
        i = 0
        while i < len(self.polynomials):
            _str += str(self.polynomials[i])
            while i < len(self.polynomials):
                i += 1
                if i >= len(self.polynomials) or self.polynomials[i] != self.polynomials[i - 1]:
                    break
                count += 1
            if count > 1:
                _str += f"^{count}"
                count = 1
        return _str

    def __repr__(self) -> str:
        return self.__str__()

    def append(self, polynomial: Polynomial | Monomial):
        self.polynomials.append(polynomial)
        self.polynomials = sorted(self.polynomials,
                                  key=lambda _polynomial: (_polynomial.get_degree(), _polynomial.coefficients))

    def get_first_polynomial(self) -> Polynomial | None:
        return self.polynomials[0] if len(self.polynomials) > 0 else None

    def __add__(self, other):
        if not isinstance(other, Polynomials):
            other = Polynomials(other)
        return Polynomials(list(self.polynomials) + list(other.polynomials), self.coefficient * other.coefficient)

    def __iadd__(self, other):
        return self.__add__(other)

    def get_polynomials(self) -> tuple[Polynomial, ...]:
        return tuple(self.polynomials)

    def get_coefficient(self) -> Monomial:
        return self.coefficient


def polynomial_factorization(polynomial: Polynomial) -> Polynomials:
    """
    因式分解多項式

    :param polynomial: 欲因式分解的多項式
    :return: 因式分解的結果
    """
    denominators: set[int] = {coefficient.get_denominator() for coefficient in polynomial.get_coefficients()}
    if denominators != {1}:
        if len(denominators) == 1:
            return Polynomials(polynomial_factorization(
                Polynomial([coefficient * denominators.pop() for coefficient in polynomial.get_coefficients()],
                           arrangement=ArrangementEnum.ASCENDING)), Monomial(Rational(1, denominators.pop())))
        denominators_lcm: int = util.lcm(tuple(denominators))
        return Polynomials(polynomial_factorization(
            Polynomial([coefficient * denominators_lcm for coefficient in polynomial.get_coefficients()],
                       arrangement=ArrangementEnum.ASCENDING)), Monomial(Rational(1, denominators_lcm)))
    degree = 0
    while polynomial.coefficients[degree] == 0:
        degree += 1
    if degree != 0:
        polynomial.coefficients = polynomial.coefficients[degree:]
    if len(polynomial) != 1:
        _gcd = util.gcd([coefficient.get_numerator() for coefficient in polynomial.get_coefficients()])
    else:
        _gcd = 1
    if polynomial.highest_degree_coefficient() < 0:
        _gcd = -_gcd
    for coefficient in polynomial.coefficients:
        if _gcd != 1:
            coefficient //= _gcd
    if _gcd != 1 or degree != 0:
        return Polynomials(
            polynomial_factorization(
                Polynomial([coefficient // _gcd for coefficient in polynomial.get_coefficients()],
                           arrangement=ArrangementEnum.ASCENDING)), Monomial(Rational(1, _gcd), degree))
    highest_degree_coefficient_factors = util.int_factorization(polynomial.highest_degree_coefficient().get_numerator(),
                                                                True)
    lowest_degree_coefficient_factors = util.int_factorization(polynomial.lowest_degree_coefficient().get_numerator())
    polynomial_factors = Polynomials()
    _found = False
    for factor1 in highest_degree_coefficient_factors:
        for factor2 in lowest_degree_coefficient_factors:
            if polynomial.test(-Rational(factor2, factor1)) and util.is_coprime(factor1, factor2):
                polynomial_factors.append(Polynomial((factor1, factor2)))
                _found = True
                break
        if _found:
            break
    if polynomial_factors.get_first_polynomial() is None:
        if polynomial.is_monomial():
            return Polynomials(coefficient=polynomial.to_monomial())
        return Polynomials(polynomial)
    polynomial_factors += polynomial_factorization(polynomial // polynomial_factors.get_first_polynomial())
    return polynomial_factors
