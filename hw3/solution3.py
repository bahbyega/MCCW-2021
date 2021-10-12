from solution2 import *
from solution1 import *


def f(x): return x * x / (1 + x * x)


def inverse_table(table):
    return [(fx, x) for x, fx in table]


def get_eps_as_10_times_neg_degree(n):
    return 10**(-n)


if __name__ == '__main__':
    # Вариант 6
    var = "f(x) = x^2 / (1+x^2)"

    print(f'''Задание 3. Задача обратного интерполирования
            ----------------------------------------------
            Вариант 6: {var}''')

    a = float(input("Введите a: "))
    b = float(input("Введите b: "))

    print(f'Отрезок: [{a}, {b}]')

    table_dim = int(input("\nВведите число значений в таблице: "))
    table = create_table(a, b, table_dim)
    print_table(table, True)

    while True:
        F = float(input(f'\nВведите f: '))

        print(f'\n-------------------- Способ #1 --------------------')
        inv_pol_deg = int(input("\nВведите степень многочлена: "))
        inv_table = inverse_table(table)
        print_table(inv_table, True)

        dist_table = dict()
        for root in inv_table:
            dist_table[abs(root[0] - F)] = root

        sorted_inv_table = dict(sorted(dist_table.items()))
        print_table(sorted_inv_table, False)

        y_lagrange = do_lagrange(
            list(sorted_inv_table.values()), F, inv_pol_deg)
        print(
            f'\nЗначение аргументов многочлена (Лагранж): {y_lagrange}')
        print(
            f'Модуль невязки: {abs(f(y_lagrange) - F)}')

        print(f'\n-------------------- Способ #2 --------------------')
        pol_deg = int(input("\nВведите степень многочлена: "))

        eps = get_eps_as_10_times_neg_degree(
            int(input("\nВведите степень eps: ")))
        #eps = 10**(-8)

        new_dist_table = dict()
        for root in table:
            new_dist_table[abs(root[0] - F)] = root

        sorted_table = dict(sorted(new_dist_table.items()))
        print_table(sorted_table, False)

        def P(x):
            return do_lagrange(list(sorted_table.values()), x, pol_deg) - F

        segments = root_finding(P, a, b, 10)
        y_secant = secant(P, a, b, eps)

        print(
            f'\nЗначение аргументов многочлена (Секущие): {y_secant}')
        print(
            f'Модуль невязки: {abs(f(y_secant) - F)}')

        if input("\nПродолжить с новыми F, n, eps? (y, n): ") == "y":
            continue
        else:
            break
