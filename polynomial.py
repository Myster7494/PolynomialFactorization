from enum import Enum

from rational import Rational
from util import int_factorization, gcd, lcm


class Monomial:
    """ 任意有理數單項式 """

    def __init__(self, coefficient: int | float | Rational = 0, degree: int = 0):
        """
        任意有理係數單項式

        :param coefficient: 係數
        :param degree: 次數
        """
        self.coefficient = Rational(coefficient)
        self.degree = degree

    def __str__(self):
        """ 輸出單項式 """
        _str = str(self.coefficient)
        if self.degree != 0:
            if self.coefficient == 1:
                _str = ""
            elif self.coefficient == -1:
                _str = "-"
            _str += "x"
            if self.degree != 1:
                _str += f"^{self.degree}"
        return _str

    def __repr__(self):
        """ 輸出單項式 """
        return self.__str__()

    def __mul__(self, other):
        """ 單項式乘法 """
        if not isinstance(other, Monomial):
            other = Monomial(other)
        return Monomial(self.coefficient * other.coefficient, self.degree + other.degree)

    def __imul__(self, other):
        """ 單項式乘法 """
        return self.__mul__(other)

    def __eq__(self, other):
        """ 判斷單項式是否相等 """
        return self.coefficient == other.coefficient and self.degree == other.degree

    def get_coefficient(self) -> Rational:
        """
        回傳單項式的係數

        :return: 係數
        """
        return self.coefficient

    def get_degree(self) -> int:
        """
        回傳單項式的次數

        :return: 次數
        """
        return self.degree


class ArrangementEnum(Enum):
    """ 多項式的排列方式 """

    DESCENDING = 0
    """ 降冪排列 """
    ASCENDING = 1
    """ 升冪排列 """


class Polynomial:
    """ 有理數多項式 """

    def __init__(self, coefficients: tuple[int | float | Rational, ...] | list[
        int | float | Rational] | int | float | Rational | Monomial = 0,
                 arrangement: ArrangementEnum = ArrangementEnum.DESCENDING):
        """
        任意整係數多項式

        :param coefficients: 係數
        :param arrangement: 排列方式
        """
        if isinstance(coefficients, (tuple, list)):
            if arrangement == ArrangementEnum.DESCENDING:
                coefficients = coefficients[::-1]
            while coefficients[-1] == 0:
                if len(coefficients) == 1:
                    break
                coefficients = coefficients[:-1]
            coefficients = list(map(Rational, coefficients))
        elif isinstance(coefficients, Monomial):
            monomial = coefficients
            coefficients = [Rational(0)] * monomial.get_degree()
            coefficients.append(monomial.get_coefficient())
        else:
            coefficients = [Rational(coefficients)]
        self.coefficients = coefficients

    def __len__(self):
        """ 回傳多項式的項數 """
        return len(self.coefficients)

    def __str__(self):
        """ 輸出多項式 """
        _str = "("
        for i in range(len(self))[::-1]:
            if self.coefficients[i] != 0:
                if self.coefficients[i] > 0 and i != self.get_degree():
                    _str += "+"
                _str += str(self.get_monomial_by_degree(i))
        _str += ")"
        return _str

    def __repr__(self):
        """ 輸出多項式 """
        return self.__str__()

    def __floordiv__(self, other):
        """ 進行綜合除法，回傳商式 """
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
        """ 進行綜合除法，回傳商式 """
        return self.__floordiv__(other)

    def __mul__(self, other):
        """ 多項式乘法 """
        result = Polynomial()
        for i, coefficient1 in enumerate(self.coefficients):
            for j, coefficient2 in enumerate(other.coefficients):
                result += Monomial(coefficient1 * coefficient2, i + j)
        return result

    def __imul__(self, other):
        """ 多項式乘法 """
        return self.__mul__(other)

    def __add__(self, other):
        """ 多項式加法 """
        result = Polynomial()
        result.coefficients = [Rational(0)] * max(len(self), len(other))
        for i in range(len(result.coefficients)):
            if i <= self.get_degree():
                result.coefficients[i] += self.coefficients[i]
            if i <= other.get_degree():
                result.coefficients[i] += other.coefficients[i]
        while result.coefficients[-1] == 0:
            if len(result.coefficients) == 1:
                break
            result.coefficients = result.coefficients[:-1]
        return result

    def __iadd__(self, other):
        """ 多項式加法 """
        return self.__add__(other)

    def __sub__(self, other):
        """ 多項式減法 """
        return self.__add__(-other)

    def __isub__(self, other):
        """ 多項式減法 """
        return self.__sub__(other)

    def __truediv__(self, other):
        """ 多項式係數真除法 """
        if not isinstance(other, Rational):
            other = Rational(other)
        return Polynomial([coefficient / other for coefficient in self.coefficients],
                          arrangement=ArrangementEnum.ASCENDING)

    def __neg__(self):
        """ 多項式變號 """
        return self * Polynomial(-1)

    def __eq__(self, other):
        """ 判斷多項式是否相等 """
        return self.coefficients == other.coefficients

    def get_coefficients(self) -> tuple[Rational, ...]:
        """
        回傳多項式的所有係數

        :return: 係數"""
        return tuple(self.coefficients)

    def get_degree(self) -> int:
        """
        回傳多項式的最高次數

        :return: 最高次數
        """
        return len(self) - 1

    def lowest_degree_coefficient(self) -> Rational:
        """
        回傳多項式的常數項

        :return: 常數項
        """
        return self.coefficients[0]

    def highest_degree_coefficient(self) -> Rational:
        """
        回傳多項式的最高次項係數

        :return: 最高次項係數
        """
        return self.coefficients[-1]

    def get_monomial_by_degree(self, degree: int) -> Monomial:
        """ 回傳多項式的某次項

        :param degree: 次數
        :return: 單項式
        """
        return Monomial(self.coefficients[degree], degree)

    def substitute(self, value: Rational | int | float) -> Rational:
        """
        代入值

        :param value: 代入的值
        :return: 代入後的值
        """
        result = Rational(0)
        for i, coefficient in enumerate(self.coefficients):
            result += coefficient * (value ** i)
        return result

    def is_divisible(self, divisor: any) -> bool:
        """
        判斷是否能被整除

        :param divisor: 除數
        :return: 是否能被整除
        """
        if isinstance(divisor, (int, float, Rational)) or (
                isinstance(divisor, (Monomial, Polynomial)) and divisor.get_degree() == 0):
            return True
        if (self - divisor * (self // divisor)) == 0:
            return True
        return False

    def is_monomial(self) -> bool:
        """
        判斷是否為單項式

        :return: 是否為單項式
        """
        for coefficient in self.coefficients[1:]:
            if coefficient != 0:
                return False
        return True

    def to_monomial(self) -> Monomial:
        """
        轉換為單項式

        :return: 單項式
        """
        return Monomial(self.highest_degree_coefficient(), self.get_degree())


class Polynomials:
    """ 多個有理數多項式 """

    def __init__(self, polynomials: any = (),
                 coefficient: Rational | Monomial | int | float = Rational(1)):
        """
        多個有理數多項式

        :param polynomials: 多項式
        :param coefficient: 係數
        """
        if not isinstance(coefficient, Monomial):
            coefficient = Monomial(coefficient)
        if isinstance(polynomials, (list, tuple)):
            coefficient = coefficient
            polynomials = list(polynomials)
        elif isinstance(polynomials, Polynomial):
            coefficient = coefficient
            polynomials = [polynomials]
        elif isinstance(polynomials, Polynomials):
            coefficient = polynomials.get_coefficient() * coefficient
            polynomials = list(polynomials.get_polynomials())
        elif isinstance(polynomials, Monomial):
            coefficient = coefficient * polynomials
            polynomials = []
        else:
            raise ValueError("Unsupported type.")
        polynomials.sort(key=lambda polynomial: (polynomial.get_degree(), polynomial.coefficients))
        self.polynomials = polynomials
        self.coefficient = coefficient

    def __str__(self):
        """ 輸出多個多項式 """
        if len(self.polynomials) == 0:
            return str(self.coefficient)
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

    def __repr__(self):
        """ 輸出多個多項式 """
        return self.__str__()

    def __mul__(self, other):
        """ 多項式乘法 """
        if not isinstance(other, Polynomials):
            other = Polynomials(other)
        return Polynomials(self.polynomials + other.polynomials, self.coefficient * other.coefficient)

    def __imul__(self, other):
        """ 多項式乘法 """
        return self.__mul__(other)

    def append(self, polynomial: Polynomial):
        """
        加入多項式

        :param polynomial: 多項式
        """
        self.polynomials.append(polynomial)
        self.polynomials.sort(key=lambda _polynomial: (_polynomial.get_degree(), _polynomial.coefficients))
        return polynomial

    def get_only_polynomial(self) -> Polynomial | None:
        """
        回傳唯一的多項式

        :return: 多項式
        """
        if len(self.polynomials) == 1:
            return self.polynomials[0]
        raise ValueError("There are not only one polynomial in the polynomials.")

    def get_polynomials(self) -> tuple[Polynomial, ...]:
        """
        回傳所有多項式

        :return: 多項式
        """
        return tuple(self.polynomials)

    def get_coefficient(self) -> Monomial:
        """
        回傳係數

        :return: 係數
        """
        return self.coefficient


def grouping_by_common_factor(polynomial: Polynomial) -> Polynomials:
    """
    提公因式

    :param polynomial: 欲提取公因式的多項式
    :return: 公因式和提取後的多項式
    """
    denominators: set[int] = {coefficient.get_denominator() for coefficient in polynomial.get_coefficients()}
    if denominators != {1}:  # 檢查是否為整係數多項式
        if len(denominators) == 1:
            return Polynomials(grouping_by_common_factor(polynomial * Polynomial(denominators.pop())),
                               Monomial(Rational(1, denominators.pop())))
        denominators_lcm: int = lcm(tuple(denominators))
        return Polynomials(grouping_by_common_factor(polynomial * Polynomial(denominators_lcm)),
                           Monomial(Rational(1, denominators_lcm)))
    degree = 0
    while polynomial.coefficients[degree] == 0:
        degree += 1
    polynomial.coefficients = polynomial.coefficients[degree:]
    _gcd = gcd([coefficient.to_int() for coefficient in polynomial.get_coefficients()])
    if polynomial.highest_degree_coefficient() < 0:  # 提出領導係數的負號
        _gcd = -_gcd
    if _gcd != 1 or degree != 0:
        return Polynomials(polynomial / _gcd, Monomial(Rational(_gcd), degree))
    return Polynomials(polynomial)


def lagrange_interpolation(
        points: tuple[tuple[int | float | Rational, int | float | Rational] | list[int | float | Rational], ...] | list[
            tuple[int | float | Rational, int | float | Rational] | list[int | float | Rational]]) -> Polynomial:
    """
    生成拉格朗日插值多項式

    :param points: 所有點的座標
    :return: 插值多項式
    """
    result = Polynomial()
    for i in range(len(points)):
        # 生成基礎多項式
        basis_polynomial = Polynomial(1)
        for j in range(len(points)):
            if i == j:
                continue
            if points[i][0] == points[j][0]:
                raise ValueError("The x positions of the points must be different.")
            basis_polynomial *= Polynomial((1, -Rational(points[j][0]))) / Rational(points[i][0] - points[j][0])
        result += basis_polynomial * Polynomial(points[i][1])
    return result


def polynomial_factor_test(polynomial: Polynomial) -> Polynomials:
    """
    使用一次因式檢驗法和拉格朗日插值法尋找因式

    :param polynomial: 欲因式分解的多項式
    :return: 因式分解的結果
    """
    if polynomial.is_monomial():
        return Polynomials(coefficient=polynomial.to_monomial())
    polynomial_factors = Polynomials()

    # 嘗試一次因式檢驗法
    highest_degree_coefficient_factors = int_factorization(polynomial.highest_degree_coefficient().to_int(),
                                                           True)
    lowest_degree_coefficient_factors = int_factorization(polynomial.lowest_degree_coefficient().to_int())

    for factor1 in highest_degree_coefficient_factors:
        for factor2 in lowest_degree_coefficient_factors:
            if polynomial.substitute(-Rational(factor2, factor1)) == 0:
                _gcd = gcd((factor1, factor2))
                factor = Polynomial((factor1 // _gcd, factor2 // _gcd))
                polynomial_factors.append(factor)
                polynomial_factors *= polynomial_factor_test(polynomial // factor)
                return polynomial_factors

    # 若多項式次數小於 4，則其因式必有一個是一次因式
    if polynomial.get_degree() < 4:
        return Polynomials(polynomial)

    # 嘗試拉格朗日插值法

    # 代入 0,1 並因數分解結果
    test_list: list[tuple[int, list[int]]] = [
        (0, int_factorization(polynomial.lowest_degree_coefficient().to_int())),
        (1, int_factorization(polynomial.substitute(1).to_int()))]

    # 必有因式次數小於等於多項式次數的一半
    for i in range(2, polynomial.get_degree() // 2 + 1):

        # 代入 -1,2,-2,3,-3,...,n,-n 並因數分解結果
        if i % 2 == 0:
            test_list.append((-(i // 2), int_factorization(polynomial.substitute(-(i // 2)).to_int())))
        else:
            test_list.append((i // 2 + 1, int_factorization(polynomial.substitute(i // 2 + 1).to_int())))

        test_point_index: list[int] = [0] * len(test_list)  # 紀錄測試點的索引值
        tested_result: list[tuple[Polynomial, bool]] = []  # 紀錄已測試過的結果，減少運算次數

        while test_point_index[-1] != len(test_list[-1][1]):
            # 準備欲生成多項式的點
            test_points = []
            for j in range(len(test_list)):
                test_points.append((test_list[j][0], test_list[j][1][test_point_index[j]]))
            # 生成拉格朗日插值多項式
            lagrange_polynomial = lagrange_interpolation(test_points)
            if lagrange_polynomial.get_degree() == i:  # 檢測多項式次數是否正確
                # 部分多項式互為倍數關係，將其視為同一多項式
                lagrange_polynomial = grouping_by_common_factor(lagrange_polynomial).get_only_polynomial()
                # 將測試結果存入列表
                _found = False
                divisible = False
                for tested_polynomial, result in tested_result:
                    if tested_polynomial == lagrange_polynomial:
                        _found = True
                        divisible = result
                        break
                if not _found:
                    divisible = polynomial.is_divisible(lagrange_polynomial)
                    tested_result.append((lagrange_polynomial, divisible))
                if divisible:
                    polynomial_factors.append(lagrange_polynomial)
                    polynomial_factors *= polynomial_factor_test(polynomial // lagrange_polynomial)
                    return polynomial_factors
            # 更新測試點索引值
            test_point_index[0] += 1
            for j in range(len(test_list) - 1):
                if test_point_index[j] == len(test_list[j][1]):
                    test_point_index[j] = 0
                    test_point_index[j + 1] += 1
    return Polynomials(polynomial)


def polynomial_factorization(polynomial: Polynomial) -> Polynomials:
    """
    因式分解多項式

    :param polynomial: 欲因式分解的多項式
    :return: 因式分解的結果
    """

    # 嘗試提公因式
    polynomial_factors = grouping_by_common_factor(polynomial)

    # 嘗試一次因式檢驗法和拉格朗日插值法
    if polynomial_factors.get_only_polynomial().get_degree() > 1:
        polynomial_factors = polynomial_factor_test(
            polynomial_factors.get_only_polynomial()) * Polynomials(
            coefficient=polynomial_factors.get_coefficient())

    return polynomial_factors
