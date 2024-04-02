from enum import Enum

from rational import Rational


class Monomial:
    """ 任意整係數單項式 """

    def __init__(self, coefficient: int | float | Rational, degree: int):
        """
        任意整係數單項式

        :param coefficient: 係數
        :param degree: 次數
        """
        self.coefficient = Rational(coefficient)
        self.power = degree

    def get_coefficient(self) -> Rational:
        """ 回傳單項式的係數 """
        return self.coefficient

    def get_degree(self) -> int:
        """ 回傳單項式的次數 """
        return self.power


class Polynomial:
    """ 任意整係數多項式 """

    def __init__(self, coefficients: tuple[int | float | Rational, ...] | list[int | float | Rational]):
        """
        任意整係數多項式

        coefficients[i] = 第 i 次的係數
        """
        self.coefficients = [Rational(coefficient) for coefficient in coefficients]

    def __len__(self) -> int:
        return len(self.coefficients)

    def __str__(self) -> str:
        return str(self.coefficients)

    def __repr__(self) -> str:
        return self.__str__()

    def get_coefficients(self) -> tuple[Rational, ...]:
        """ 回傳多項式的所有係數 """
        return tuple(self.coefficients)

    def get_degree(self) -> int | Rational:
        """ 回傳多項式的最高次數 """
        return len(self) - 1

    def lower_degree_coefficient(self) -> Rational:
        """ 回傳多項式的常數項 """
        return self.coefficients[0]

    def highest_degree_coefficient(self) -> Rational:
        """ 回傳多項式的最高次項係數 """
        return self.coefficients[-1]

    def get_monomial_by_degree(self, degree: int) -> Monomial:
        """ 回傳多項式的 int degree 次項 """
        return Monomial(self.coefficients[degree], degree)

    def test(self, value: int) -> bool:
        """ 測試 int value 代入後是否為 0 """
        result = 0
        for i, coefficient in enumerate(self.coefficients):
            result += coefficient * (value ** i)
        return result == 0


class Quadratic(Polynomial):
    """ 二次整係數多項式 """

    class discriminant_enum(Enum):
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

    def discriminant(self) -> discriminant_enum:
        """
        回傳判別式的結果

        :return: 判別式的結果
        """
        result = self.get_b() ** 2 - Rational(4) * self.get_a() * self.get_c()
        if result < 0:
            return self.discriminant_enum.NO_REAL_ROOT
        if result == 0:
            return self.discriminant_enum.TWO_SAME_ROOTS
        if result > 0:
            return self.discriminant_enum.TWO_DIFFERENT_ROOTS


def synthetic_division(dividend: Polynomial, divisor: Polynomial) -> Polynomial:
    """
    進行綜合除法，回傳商式
    :param dividend: 被除式
    :param divisor: 除式
    :return: 商式
    """
    dividend_coefficients = dividend.get_coefficients()[::-1]
    leading_coefficient = divisor.highest_degree_coefficient()
    divisor_coefficients = [-Rational(coefficient, leading_coefficient) for coefficient in
                            divisor.get_coefficients()[:-1]]
    coefficients_difference = len(dividend_coefficients) - len(divisor_coefficients)
    table: list[list[Rational]] = [[Rational(0)] * coefficients_difference for _ in
                                   range(len(divisor_coefficients))]  # row為直，column為橫
    quotient_coefficients: list[int | Rational] = [dividend_coefficients[i] for i in range(coefficients_difference)]
    for i in range(coefficients_difference):
        for j in range(len(divisor_coefficients)):
            if i - j - 1 >= 0:
                quotient_coefficients[i] = table[j][i - j - 1] + quotient_coefficients[i]
        for j, divisor_coefficient in enumerate(divisor_coefficients):
            table[j][i] = divisor_coefficient * quotient_coefficients[i]
    print(table)
    return Polynomial([(quotient_coefficient if isinstance(quotient_coefficient, Rational) else Rational(
        quotient_coefficient)) / (leading_coefficient if isinstance(leading_coefficient, Rational) else Rational(
        leading_coefficient)) for quotient_coefficient in
                       quotient_coefficients])
