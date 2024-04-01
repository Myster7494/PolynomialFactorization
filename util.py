from fraction import Fraction
from polynomial import Polynomial


def int_factorization(value: int, only_nature=False) -> set[int]:
    """
    傳入一個 int value，回傳 value 的所有因數

    :param value: 要因數分解分解的正整數
    :param only_nature: 只回傳正因數
    """
    if value == 0:
        raise ValueError("0 has infinite factors.")
    factors: set[int] = set()
    for i in range(1, int(value ** 0.5) + 1):
        if value % i == 0:
            factors.add(i)
            factors.add(abs(value // i))
            if not only_nature:
                factors.add(-i)
                factors.add(-abs(value // i))
    return factors


def gcd(value1: int, value2: int) -> int:
    """
    求兩個 int value 的最大公因數

    :param value1: 第一個整數
    :param value2: 第二個整數
    """
    if value2 == 0:
        if value1 == 0:
            raise ValueError("Both values are 0.")
        return abs(value1)
    return gcd(value2, value1 % value2)


def lcm(value1: int, value2: int) -> int:
    """
    求兩個 int value 的最小公倍數

    :param value1: 第一個整數
    :param value2: 第二個整數
    """
    return abs(value1 * value2) // gcd(value1, value2)


def is_coprime(value1: int, value2: int) -> bool:
    """
    測試兩個 int value 是否互質

    :param value1: 第一個整數
    :param value2: 第二個整數
    """
    return gcd(value1, value2) == 1


def synthetic_division(dividend: Polynomial, divisor: Polynomial) -> Polynomial:
    """
    進行綜合除法，回傳商式

    :param dividend: 被除式
    :param divisor: 除式
    :return: 商式
    """
    dividend_coefficients = dividend.get_coefficients()
    leading_coefficient = divisor.highest_degree_coefficient()
    divisor_coefficients = [-Fraction(coefficient, leading_coefficient) for coefficient in
                            list(divisor.get_coefficients())[:-1]]
    print(divisor_coefficients)
    coefficients_difference = len(dividend_coefficients) - len(divisor_coefficients)
    table: list[list[int | Fraction]] = [[0] * coefficients_difference for _ in
                                         range(len(divisor_coefficients))]  # row為直，column為橫
    quotient_coefficients: list[int | Fraction] = [dividend_coefficients[i] for i in range(coefficients_difference)]
    for i in range(coefficients_difference):
        for j in range(len(divisor_coefficients)):
            quotient_coefficients[i] = table[j][i-j-1] + quotient_coefficients[i]
        for j, divisor_coefficient in enumerate(divisor_coefficients):
            table[j][i] = divisor_coefficient * quotient_coefficients[i]
    print(table)
    print([quotient_coefficient / leading_coefficient for quotient_coefficient in quotient_coefficients])
