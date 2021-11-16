import math
import scipy.integrate as integrate
from tabulate import tabulate


# АСТ 0
def qf_left_rect(f, a, b):
    return (b - a) * f(a)


# АСТ 0
def qf_right_rect(f, a, b):
    return (b - a) * f(b)


# АСТ 1
def qf_inter_rect(f, a, b):
    return (b - a) * f((a + b) / 2)


# АСТ 1
def qf_trap(f, a, b):
    return (b - a) * (f(a) + f(b)) / 2


# АСТ 3
def qf_simp(f, a, b):
    return (b - a) * (
        f(a) +
        4 * f((a + b) / 2) +
        f(b)
    ) / 6


# АСТ 3
def qf_3_8(f, a, b):
    h = (b - a) / 3
    return (b - a) * (
        f(a) +
        3 * f(a + h) +
        3 * f(a + 2 * h) +
        f(b)
    ) / 8


def main_loop():

    def quad_f(f, a, b):
        integral, _ = integrate.quad(f, a, b)
        return integral

    functions = [
        (lambda x: x * x / (1 + x * x), "x^2 / (1 + x^2)"),
        (lambda x: x * 0 + 1, "1"),
        (lambda x: x, "x"),
        (lambda x: x ** 2, "x^2"),
        (lambda x: x ** 3, "x^3"),
        (lambda x: math.e ** ((-x)**2), "e^(-x^2)")
    ]

    quadr_formulas = [
        qf_left_rect,
        qf_right_rect,
        qf_inter_rect,
        qf_trap,
        qf_simp,
        qf_3_8
    ]

    while True:
        a = float(input("Введите a: "))
        b = float(input("Введите b: "))

        for f, f_name in functions:
            corr_integral = quad_f(f, a, b)
            appr_integrals = [qf(f, a, b) for qf in quadr_formulas]
            errors = [abs(y - corr_integral) for y in appr_integrals]

            qf_names = [qf.__name__ for qf in quadr_formulas]

            print(f'\nf(x) = {f_name}')
            print(f'Exact integral: {corr_integral}\n')
            print(tabulate([(qf_name, y, err) for qf_name, y, err in
                            zip(qf_names, appr_integrals, errors)],
                           headers=["quad_formula", "appr", "err"],
                           floatfmt=".10f"))
            print()

        if input("\nПродолжить с новыми a, b? (y, n): ") == "y":
            continue
        else:
            break


if __name__ == '__main__':
    print(f'''Задание 4. Приближенное вычисление интеграла по квадратурным формулам
---------------------------------------------------------------------''')
    main_loop()
