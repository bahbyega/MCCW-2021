import scipy.integrate as integrate
from tabulate import tabulate

from solution4_2 import *


def clarify_by_runge(l, h_l, h, qf, f, a, n):
    if qf.__name__ in ["qf_left_rect", "qf_right_rect"]:
        ast = 0
    elif qf.__name__ in ["qf_inter_rect", "qf_trap"]:
        ast = 1
    else:
        ast = 3

    Jh = qf(f, a, n, h)
    Jhl = qf(f, a, n * l, h_l)
    r = ast + 1

    return (l ** r * Jhl - Jh) / (l ** r - 1)


def main_loop():

    while True:
        a = float(input("Введите a: "))
        b = float(input("Введите b: "))
        n = int(input("Введите число разбиений промежутка интегрирования: "))
        l = int(input("Увеличить m в l раз: "))

        h = (b - a) / n
        h_l = (b - a) / (n * l)

        print()
        print(f"h = (b-a)/m = {h}")
        print(f"h/l = (b-a)/(m*l) = {h_l}")

        for f, f_name in functions:
            corr_integral = quad_f(f, a, b)
            appr_integrals = [qf(f, a, n, h) for qf in quadr_formulas]
            appr_integrals_l = [qf(f, a, n * l, h_l) for qf in quadr_formulas]

            errs = [abs(y - corr_integral) for y in appr_integrals]
            errs_l = [abs(y - corr_integral) for y in appr_integrals_l]

            appr_integrals_clar = [
                clarify_by_runge(l, h_l, h, qf, f, a, n) for qf in quadr_formulas
            ]
            errs_clar = [abs(y - corr_integral) for y in appr_integrals_clar]

            qf_names = [qf.__name__ for qf in quadr_formulas]

            print(f"\nf(x) = {f_name}")
            print(f"J = {corr_integral}\n")
            print(
                tabulate(
                    [
                        (qf_name, y, err, y_l, err_l, appr_clar, err_clar)
                        for qf_name, y, err, y_l, err_l, appr_clar, err_clar in zip(
                            qf_names,
                            appr_integrals,
                            errs,
                            appr_integrals_l,
                            errs_l,
                            appr_integrals_clar,
                            errs_clar,
                        )
                    ],
                    headers=[
                        "quad_formula",
                        "J(h)",
                        "|J - J(h)|",
                        "J(h/l)",
                        "|J - J(h/l)|",
                        "J_r",
                        "|J - J_r|",
                    ],
                    floatfmt=".2e",
                )
            )
            print("")

        if input("\nПродолжить с новыми a, b, m, l? (y, n): ") == "y":
            continue
        else:
            break


if __name__ == "__main__":
    print(
        f"""Задание 4_3. Приближенное вычисление интеграла по составным квадратурным формулам
---------------------------------------------------------------------"""
    )
    main_loop()
