import util


class Fraction:
    """ 一般分數 """

    def __init__(self, numerator: int | float, denominator: int | float = 1):
        """
        一般分數

        :param numerator: 分子或小數
        :param denominator: 分母
        """
        if isinstance(numerator, float) or isinstance(denominator, float):
            numerator, denominator = (numerator / denominator).as_integer_ratio()
        reduction = util.gcd(numerator, denominator)
        self.numerator = numerator // reduction
        self.denominator = denominator // reduction

    def __str__(self):
        if self.denominator == 1:
            if isinstance(self.numerator, float) and self.numerator.is_integer():
                return str(int(self.numerator))
            return str(self.numerator)
        return f"{self.numerator}/{self.denominator}"

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        """ 分數加法 """
        if not isinstance(other, self.__class__):
            other = Fraction(other)
        lcd = util.lcm(self.denominator, other.denominator)
        return Fraction(self.numerator * (lcd // self.denominator) + other.numerator * (lcd // other.denominator), lcd)

    def __sub__(self, other):
        """ 分數減法 """
        return self.__add__(-other)

    def __mul__(self, other):
        """ 分數乘法 """
        if not isinstance(other, self.__class__):
            other = Fraction(other)
        reduction1 = util.gcd(self.numerator, other.denominator)
        reduction2 = util.gcd(self.denominator, other.numerator)
        return Fraction((self.numerator // reduction1) * (other.numerator // reduction2),
                        (self.denominator // reduction2) * (other.denominator // reduction1))

    def __truediv__(self, other):
        """ 分數除法 """
        return self.__mul__(other ** -1)

    def __pow__(self, power: int):
        """ 分數的整數次方 """
        if power < 0:
            return Fraction(self.denominator ** abs(power), self.numerator ** abs(power))
        return Fraction(self.numerator ** power, self.denominator ** power)

    def __neg__(self):
        """ 分數變號 """
        return Fraction(-self.numerator, self.denominator)
