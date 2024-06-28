"""
Date: 2024/06
Author: 桃園高中 邱顯智
"""
import random

from polynomial import Polynomial, polynomial_factorization
from rational import Rational

# return polynomial_factorization(Polynomial([Rational(x) for x in input(
#     "輸入欲因式分解之多項式之降冪排列係數，可將結果放入wolframalpha得到更好的顯示效果：").split()]))


if __name__ == "__main__":
    if input("手動輸入或自動測試(Y/n)：").strip().lower() != "n":
        origin_polynomial = Polynomial([Rational(x) for x in input(
            "輸入欲因式分解之多項式之降冪排列係數，可將結果放入wolframalpha得到更好的顯示效果：").split()])
        print(origin_polynomial, "=", sep="", end="")
        print(polynomial_factorization(origin_polynomial))
    else:
        print("自動生成隨機多項式，可放入wolframalpha驗證：")
        for i in range(3):
            origin_polynomial = Polynomial(1)
            for j in range(random.randint(1, 2)):
                origin_polynomial *= Polynomial(
                    [Rational(random.randint(-1000, 1000), random.randint(-1000, 1000)) for _ in
                     range(random.randint(2, 3))])
            print(origin_polynomial, "=", sep="", end="")
            print(polynomial_factorization(origin_polynomial))
