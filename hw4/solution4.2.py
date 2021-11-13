import math
import scipy.integrate as integrate
from tabulate import tabulate


# АСТ 0
def qf_left_rect(f, a, n, h):
    sum = 0

    for i in range(1, n + 1):
        sum += f(a + (i - 1) * h)

    return sum * h


# АСТ 0
def qf_right_rect(f, a, n, h):
    sum = 0

    for i in range(1, n + 1):
        sum += f(a + h + (i - 1) * h)

    return sum * h


# АСТ 1
def qf_inter_rect(f, a, n, h):
    sum = 0

    for i in range(1, n + 1):
        sum += f(a + h/2 + (i - 1) * h)

    return sum * h


# АСТ 1
def qf_trap(f, a, n, h):
    sum = 0

    for i in range(1, n + 1):
        sum += f(a + (i - 1) * h) \
            + f(a + h * i)

    return sum * h/2


# АСТ 3
def qf_simp(f, a, n, h):
    sum = 0

    for i in range(1, n + 1):
        sum += f(a + (i - 1) * h) \
            + 4 * f(a + h/2 * (2 * i - 1)) \
            + f(a + h * i)

    return sum * h/6


# АСТ 3
def qf_3_8(f, a, n, h):
    sum = 0

    for i in range(1, n + 1):
        sum += f(a + (i - 1) * h) \
            + 3 * f(a + i * h - 2 * h/3) \
            + 3 * f(a + i * h - h/3) \
            + f(a + h * i)

    return sum * h/8


def quad_f(f, a, b):
    integral, _ = integrate.quad(f, a, b)
    return integral


def find_max(f, a, n, h):
    m = 0
    for i in range(n + 1):
        m = max(m, abs(f(a + h * i)))
    return m


def theor_err(qf, f, a, b, n, h):
    temp = find_max(f, a, n, h) * (b - a)

    if qf.__name__ == "qf_left_rect" \
            or qf.__name__ == "qf_right_rect":
        return temp * h / 2
    elif qf.__name__ == "qf_inter_rect":
        return temp * h ** 2 / 24
    elif qf.__name__ == "qf_trap":
        return temp * h ** 2 / 12
    elif qf.__name__ == "qf_simp":
        return temp * h ** 4 / 2880
    else:
        return temp * h ** 4 / 6480


def main_loop():

    functions = [
        (lambda x: x * x / (1 + x * x), "x^2 / (1 + x^2)"),
        (lambda x: x * 0 + 1, "1"),
        (lambda x: x * 2, "2x"),
        (lambda x: x ** 2 * 3, "3x^2"),
        (lambda x: x ** 3 * 4, "4x^3"),
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
        n = int(input("Введите число разбиений промежутка интегрирования: "))

        h = (b - a) / n

        for f, f_name in functions:
            corr_integral = quad_f(f, a, b)
            appr_integrals = [qf(f, a, n, h) for qf in quadr_formulas]
            errs = [abs(y - corr_integral) for y in appr_integrals]
            theor_errs = [theor_err(qf, f, a, b, n, h)
                          for qf in quadr_formulas]

            qf_names = [qf.__name__ for qf in quadr_formulas]

            print(f'\nf(x) = {f_name}')
            print(f'Exact integral: {corr_integral}\n')
            print(tabulate([(qf_name, y, err, terr) for qf_name, y, err, terr in
                            zip(qf_names, appr_integrals, errs, theor_errs)],
                           headers=["quad_formula",
                                    "appr", "err", "theor_err"],
                           floatfmt=".10f"))
            print()

        if input("\nПродолжить с новыми a, b, m? (y, n): ") == "y":
            continue
        else:
            break


if __name__ == '__main__':
    print(f'''Задание 4. Приближенное вычисление интеграла по квадратурным формулам
---------------------------------------------------------------------''')
    main_loop()
