from enum import Enum


class Monomial:
    """ 任意整係數單項式 """

    def __init__(self, coefficient: int, degree: int):
        """
        任意整係數單項式

        :param coefficient: 係數
        :param degree: 次數
        """
        self.coefficient = coefficient
        self.power = degree

    def get_coefficient(self) -> int:
        """ 回傳單項式的係數 """
        return self.coefficient

    def get_degree(self) -> int:
        """ 回傳單項式的次數 """
        return self.power


class Polynomial:
    """ 任意整係數多項式 """

    def __init__(self, coefficients: tuple[int, ...]):
        """
        任意整係數多項式

        coefficients[i] = 第 i 次的係數
        """
        self.coefficients = coefficients

    def __len__(self) -> int:
        return len(self.coefficients)

    def get_degree(self) -> int:
        """ 回傳多項式的最高次數 """
        return len(self) - 1

    def lower_degree_coefficient(self) -> int:
        """ 回傳多項式的常數項 """
        return self.coefficients[0]

    def highest_degree_coefficient(self) -> int:
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

    def __init__(self, a: int, b: int, c: int):
        super().__init__((a, b, c))

    def get_a(self) -> int:
        """ 回傳二次多項式的二次項係數 """
        return self.coefficients[2]

    def get_b(self) -> int:
        """ 回傳二次多項式的一次項係數 """
        return self.coefficients[1]

    def get_c(self) -> int:
        """ 回傳二次多項式的常數項 """
        return self.coefficients[0]

    def discriminant(self) -> discriminant_enum:
        """
        回傳判別式的結果

        :return: 判別式的結果
        """
        result = self.get_b() ** 2 - 4 * self.get_a() * self.get_c()
        if result < 0:
            return self.discriminant_enum.NO_REAL_ROOT
        if result == 0:
            return self.discriminant_enum.TWO_SAME_ROOTS
        if result > 0:
            return self.discriminant_enum.TWO_DIFFERENT_ROOTS
