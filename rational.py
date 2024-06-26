"""
Date: 2024/06
Author: 桃園高中 邱顯智
"""

import util


class Rational:
    """ 有理數 """

    def __init__(self, numerator: any, denominator: any = 1):
        """
        有理數

        :param numerator: 分子
        :param denominator: 分母，預設為 1
        """
        if isinstance(numerator, str):
            if "/" in numerator:
                numerator, denominator = map(Rational, numerator.split("/"))
            else:
                try:
                    numerator = float(numerator)
                except ValueError:
                    raise ValueError("Invalid input.")
        if not isinstance(numerator, (int, float, Rational)) or not isinstance(denominator, (int, float, Rational)):
            raise TypeError("Numerator or denominator type unsupported.Only support int, float and Rational.")
        if denominator == 0:
            raise ValueError("Denominator can't be 0.")
        if isinstance(numerator, float) or isinstance(denominator, float):
            numerator, denominator = (numerator / denominator).as_integer_ratio()
        elif isinstance(numerator, Rational) and denominator == 1:
            numerator, denominator = numerator.numerator_and_denominator()
        elif isinstance(numerator, Rational) or isinstance(denominator, Rational):
            numerator, denominator = Rational(Rational(numerator) / Rational(denominator)).numerator_and_denominator()
        if not util.is_coprime(numerator, denominator):
            numerator //= (_gcd := util.gcd((numerator, denominator)))
            denominator //= _gcd
        self.numerator = numerator if denominator > 0 else -numerator
        self.denominator = abs(denominator)

    def __str__(self):
        """ 輸出有理數 """
        if self.is_integer():
            if isinstance(self.numerator, float) and self.numerator.is_integer():
                return str(int(self.numerator))
            return str(self.numerator)
        return f"({self.numerator}/{self.denominator})" if self > 0 else f"-({-self.numerator}/{self.denominator})"

    def __repr__(self):
        """ 輸出有理數 """
        return self.__str__()

    def __add__(self, other):
        """ 有理數加法 """
        lcd = util.lcm((self.denominator, other.denominator))
        return Rational(self.numerator * (lcd // self.denominator) + other.numerator * (lcd // other.denominator), lcd)

    def __iadd__(self, other):
        """ 有理數加法 """
        return self.__add__(other)

    def __sub__(self, other):
        """ 有理數減法 """
        return self.__add__(-other)

    def __isub__(self, other):
        """ 有理數減法 """
        return self.__sub__(other)

    def __mul__(self, other):
        """ 有理數乘法 """
        return Rational((self.numerator // (gcd1 := util.gcd((self.numerator, other.denominator)))) * (
                other.numerator // (gcd2 := util.gcd((self.denominator, other.numerator)))),
                        (self.denominator // gcd2) * (other.denominator // gcd1))

    def __imul__(self, other):
        """ 有理數乘法 """
        return self.__mul__(other)

    def __truediv__(self, other):
        """ 有理數除法 """
        return self.__mul__(other ** -1)

    def __floordiv__(self, other):
        """ 有理數整除法 """
        if isinstance(other, int) and self.is_integer():
            return self.numerator // other
        numerator, denominator = self.__truediv__(other).numerator_and_denominator()
        return numerator // denominator

    def __pow__(self, power: int):
        """
        有理數的整數次方

        :param power: 次方數
        """
        if power < 0:
            return Rational(self.denominator ** abs(power), self.numerator ** abs(power))
        return Rational(self.numerator ** power, self.denominator ** power)

    def __neg__(self):
        """有理數變號 """
        return Rational(-self.numerator, self.denominator)

    def __eq__(self, other):
        """ 判斷是否相等 """
        return self.numerator == other.numerator and self.denominator == other.denominator

    def __ge__(self, other):
        """ 判斷是否大於等於 """
        return self.numerator * other.denominator >= other.numerator * self.denominator

    def __gt__(self, other):
        """ 判斷是否大於 """
        return self.numerator * other.denominator > other.numerator * self.denominator

    def __le__(self, other):
        """ 判斷是否小於等於 """
        return self.numerator * other.denominator <= other.numerator * self.denominator

    def __lt__(self, other):
        """ 判斷是否小於 """
        return self.numerator * other.denominator < other.numerator * self.denominator

    def numerator_and_denominator(self) -> tuple[int, int]:
        """
        回傳分子與分母

        :return: 分子與分母
        """
        return self.numerator, self.denominator

    def is_integer(self) -> bool:
        """
        判斷是否為整數

        :return: 是否為整數
        """
        return self.denominator == 1

    def get_numerator(self) -> int:
        """
        回傳分子

        :return: 分子
        """
        return self.numerator

    def get_denominator(self) -> int:
        """
        回傳分母

        :return: 分母
        """
        return self.denominator

    def to_int(self) -> int:
        """
        轉換為整數

        :return: 整數
        """
        if self.is_integer():
            return self.numerator
        raise ValueError("Can't convert to integer.")
