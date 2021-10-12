from tabulate import tabulate


def f(x): return x * x / (1 + x * x)


def create_table(left, right, table_dim):
    table = list()

    for i in range(table_dim):
        x_val = left + (right - left) * i / (table_dim - 1)
        table.append((x_val, f(x_val)))

    return table


def print_table(table, init_table_flag):
    if init_table_flag:
        print(f'\nИсходная таблица значений: ')
        print(tabulate(table, headers=["x", "f(x)"]))
    else:
        print(f'\nОтсортированная таблица: ')
        print(tabulate(table.items(), headers=["|x - x_0|", "( x, f(x) )"]))


def get_pol_deg_less_than_table_dim():
    while True:
        n = int(input(f'\nВведите степень многочлена: '))
        if (n >= table_dim):
            continue
        else:
            return n


def do_lagrange(table, x, pol_deg):
    y = 0

    for i in range(pol_deg):
        pol_koef = 1

        for j in range(pol_deg):
            x_j = table[j][0]
            x_i = table[i][0]

            if i != j:
                pol_koef = pol_koef * (x - x_j) / (x_i - x_j)

        y_i = table[i][1]
        y = y + pol_koef * y_i

    return y


def do_newton(table, x, pol_deg):
    y = 0

    diff_table = []

    for i in range(pol_deg + 1):
        diff_table.append([table[i][1]])

    for j in range(pol_deg):
        for i in range(pol_deg - j):
            diff_value = (diff_table[i + 1][-1] - diff_table[i]
                          [-1]) / (table[i + j + 1][0] - table[i][0])
            diff_table[i].append(diff_value)

    for k in range(pol_deg + 1):
        y_x_multiply = 1
        for j in range(k):
            y_x_multiply *= (x - table[j][0])
        y += diff_table[0][k] * y_x_multiply

    return y


def main_loop(table, x=None):
    while True:
        print_table(table, True)

        if x is None:
            x = float(input(f'\nВведите x: '))

        dist_table = dict()
        for root in table:
            dist_table[abs(root[0] - x)] = root

        sorted_table = dict(sorted(dist_table.items()))
        print_table(sorted_table, False)

        pol_deg = get_pol_deg_less_than_table_dim()

        y_lagrange = do_lagrange(list(sorted_table.values()), x, pol_deg)
        y_newton = do_newton(list(sorted_table.values()), x, pol_deg)

        print(
            f'\nЗначение интеполяционного многочлена (Лагранж): {y_lagrange}')
        print(
            f'Значение абсолютной факт. погрешности (Лагранж): {abs(f(x) - y_lagrange)}')

        print(f'\nЗначение интеполяционного многочлена (Ньютон): {y_newton}')
        print(
            f'Значение абсолютной факт. погрешности (Ньютон): {abs(f(x) - y_newton)}')

        if input("\nПродолжить с новыми x, n? (y, n): ") == "y":
            continue
        else:
            break


if __name__ == "__main__":
    # Вариант 6
    var = "f(x) = x^2 / (1+x^2)"

    a = 0.4
    b = 1

    x = 0.85

    print(f'''Задание 2. Задача алгебраического интерполирования
            ----------------------------------------------
            Вариант 6: {var}, [{a}, {b}]''')

    table_dim = int(input("\nВведите число значений в таблице: "))

    table = create_table(a, b, table_dim)
    main_loop(table)
