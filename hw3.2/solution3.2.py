from dataclasses import astuple, dataclass
from tabulate import tabulate


@dataclass
class Node:
    x: float
    y: float
    dx: float
    dx_Oh: float
    err_dx_Oh: float
    dx_Oh2: float
    err_dx_Oh2: float
    d2x: float
    d2x_Oh2: float
    err_d2x_Oh2: float


def create_table(left, right, table_dim):
    table = list()

    for i in range(table_dim):
        x_val = left + (right - left) * i / (table_dim - 1)
        table.append(Node(x_val, f(x_val), None, None,
                     None, None, None, None, None, None))

    return table


def print_table(table):
    print(f'\nТаблица: ')
    print(tabulate([astuple(x) for x in table], headers=["x", "f(x)", "f'(x)", "f'~ O(h)",
                                                         "погр. O(h)", "f'~~ O(h^2)",
                                                         "погр. O(h^2)", "f''(x)",
                                                         "f''~ O(h^2)", "погр. O(h^2)"]))


def compute_dx_with_error_Oh1(table, error):
    table[0].dx = dfx(table[0].x)
    table[0].dx_Oh = (table[1].y - table[0].y) / error
    table[0].err_dx_Oh = abs(table[0].dx - table[0].dx_Oh)

    for i in range(1, len(table)):
        table[i].dx = dfx(table[i].x)
        table[i].dx_Oh = (table[i].y - table[i - 1].y) / error
        table[i].err_dx_Oh = abs(table[i].dx - table[i].dx_Oh)


def compute_dx_with_error_Oh2(table, error):
    table[0].dx = dfx(table[0].x)
    table[0].dx_Oh2 = (-3*table[0].y + 4*table[1].y - table[2].y) / (2*error)
    table[0].err_dx_Oh2 = abs(table[0].dx - table[0].dx_Oh2)

    for i in range(1, len(table) - 1):
        table[i].dx = dfx(table[i].x)
        table[i].dx_Oh2 = (table[i+1].y - table[i-1].y) / (2*error)
        table[i].err_dx_Oh2 = abs(table[i].dx - table[i].dx_Oh2)

    table[-1].dx = dfx(table[-1].x)
    table[-1].dx_Oh2 = (3*table[-1].y - 4*table[-2].y +
                        table[-3].y) / (2*error)
    table[-1].err_dx_Oh2 = abs(table[-1].dx - table[-1].dx_Oh2)


def compute_d2x_with_error(table, error):
    for i in range(1, len(table) - 1):
        table[i].d2x = d2fx(table[i].x)
        table[i].d2x_Oh2 = (table[i+1].y - 2*table[i].y +
                            table[i-1].y) / error**2
        table[i].err_d2x_Oh2 = abs(table[i].d2x - table[i].d2x_Oh2)


def main_loop():
    while True:
        a = float(input(f'Введите a: '))
        b = float(input(f'Введите b: '))
        h = float(input(f'Введите h (шаг): '))

        table = create_table(a, b, int((b-a)/h) + 1)

        compute_dx_with_error_Oh1(table, h)
        compute_dx_with_error_Oh2(table, h)
        compute_d2x_with_error(table, h)

        print_table(table)

        if input("\nContinue?: (y, n): ") == "y":
            continue
        else:
            break


if __name__ == '__main__':

    # Вариант 6
    var = "f(x) = x^2 / (1+x^2)"
    def f(x): return x * x / (1 + x * x)
    def dfx(x): return 2 * x / (1 + x * x)**2
    def d2fx(x): return (2 - 6*x*x) / (1 + x * x)**3

    # def f(x): return x+3
    # def dfx(x): return 1
    # def d2fx(x): return 0

    print(f'''Задание 3.2. Численное дифференцирование.
            ----------------------------------------------
            Вариант 6: {var}''')

    main_loop()
