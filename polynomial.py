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
        return Monomial(self.coefficient * other.coefficient, self.degree + other.degree)

    def __imul__(self, other):
        return self.__mul__(other)


class Polynomial:
    """ 任意整係數多項式 """

    class ArrangementEnum(Enum):
        """ 多項式的排列方式 """
        DESCENDING = 0
        """ 降冪排列 """
        ASCENDING = 1
        """ 升冪排列 """

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
            if arrangement == self.ArrangementEnum.DESCENDING:
                coefficients = coefficients[::-1]
            while coefficients[-1] == 0:
                coefficients = coefficients[:-1]
            coefficients = [Rational(coefficient) for coefficient in coefficients]
        elif monomial is not None:
            coefficients = [Rational(0)] * monomial.get_degree()
            coefficients.append(monomial.get_coefficient())
        else:
            coefficients = [Rational(coefficients)]
        self.coefficients = coefficients

    def __len__(self) -> int:
        return len(self.coefficients)

    def __str__(self) -> str:
        _str = ""
        for i in range(len(self))[::-1]:
            if self.coefficients[i] != 0:
                if self.coefficients[i] > 0 and i != self.get_degree():
                    _str += "+"
                _str += str(self.get_monomial_by_degree(i))
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

    def __init__(self, a: int | float | Rational, b: int | float | Rational, c: int | float | Rational):
        super().__init__((a, b, c))

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


def polynomial_factorization(polynomial: Polynomial) -> list[Monomial | Polynomial]:
    """
    因式分解多項式

    :param polynomial: 欲因式分解的多項式
    :return: 因式分解的結果
    """
    denominators: set[int] = {coefficient.get_denominator() for coefficient in polynomial.get_coefficients()}
    if denominators != {1}:
        if len(denominators) == 1:
            return [Monomial(Rational(1, denominators.pop()))] + polynomial_factorization(
                Polynomial([coefficient * denominators.pop() for coefficient in polynomial.get_coefficients()]))
        denominators_lcm: int = util.lcm(denominators)
        return [Monomial(Rational(1, denominators_lcm))] + polynomial_factorization(
            Polynomial([coefficient * denominators_lcm for coefficient in polynomial.get_coefficients()]))

    highest_degree_coefficient_factors = util.int_factorization(polynomial.highest_degree_coefficient().get_numerator(),
                                                                True)
    lowest_degree_coefficient_factors = util.int_factorization(polynomial.lowest_degree_coefficient().get_numerator())
    polynomial_factors: list[Polynomial | Monomial] = []
    for factor1 in highest_degree_coefficient_factors:
        for factor2 in lowest_degree_coefficient_factors:
            if polynomial.test(-Rational(factor2, factor1)) and util.is_coprime(factor1, factor2):
                polynomial_factors.append(Polynomial((factor1, factor2)))
                break
        else:
            continue
        break
    if not polynomial_factors:
        if polynomial.is_monomial():
            return [polynomial.to_monomial()]
        return [polynomial]
    polynomial_factors += polynomial_factorization(polynomial // polynomial_factors[0])
    monomial = Monomial(1)
    i = 0
    while i < len(polynomial_factors):
        if isinstance(polynomial_factors[i], Monomial):
            monomial *= polynomial_factors.pop(i)
        else:
            i += 1
    if monomial.get_coefficient() != 1 or monomial.get_degree() != 0:
        polynomial_factors.insert(0, monomial)
    return polynomial_factors
