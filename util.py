def int_factorization(value: int, only_nature=False) -> set[int]:
    """ 傳入一個 int value，回傳 value 的所有因數 """
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


def is_coprime(value1: int, value2: int) -> bool:
    """ 測試兩個 int value 是否互質 """
    return int_factorization(abs(value1), only_nature=True).union(
        int_factorization(abs(value2), only_nature=True)) == {1}
