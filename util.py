def int_factorization(value: int, only_nature=False) -> set[int]:
    """
    傳入一個 int value，回傳 value 的所有因數

    :param value: 要因數分解分解的正整數
    :param only_nature: 只回傳正因數
    """
    if value == 0:
        raise ValueError("0 has infinite factors.")
    factors: set[int] = set()
    for i in range(1, int(abs(value) ** 0.5) + 1):
        if value % i == 0:
            factors.add(i)
            factors.add(abs(value // i))
            if not only_nature:
                factors.add(-i)
                factors.add(-abs(value // i))
    return factors


def gcd(values: tuple[int, ...] | list[int]) -> int:
    """
    求最大公因數

    :param values: 欲求最大公因數的所有整數
    """
    if len(values) < 2:
        raise ValueError("At least 2 values are required.")
    if len(values) == 2:
        if values[1] == 0:
            if values[0] == 0:
                raise ValueError("Both values are 0.")
            return abs(values[0])
        return gcd((values[1], values[0] % values[1]))
    return gcd((values[0], gcd(values[1:])))


def lcm(values: tuple[int, ...] | list[int]) -> int:
    """
    求最小公倍數

    :param values: 欲求最小公倍數的所有整數
    """
    if len(values) < 2:
        raise ValueError("At least 2 values are required.")

    if len(values) == 2:
        return abs(values[0] * values[1]) // gcd((values[0], values[1]))
    else:
        return lcm((values[0], lcm(values[1:])))


def is_coprime(value1: int, value2: int) -> bool:
    """
    測試兩個 int value 是否互質

    :param value1: 第一個整數
    :param value2: 第二個整數
    """
    return gcd((value1, value2)) == 1
