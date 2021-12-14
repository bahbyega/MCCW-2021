"""
This module contains solutions for sixth laboratory task,
approximate evaluation of integrals with composite quadratic formulas
for Methods of Computation and Computational Workshop course.
"""

import math
import numpy as np
import scipy.integrate as integrate

from typing import Callable, Sequence, Tuple, Final

from hw1.solution1 import root_finding, secant
from hw5.solution5 import (
    find_gauss_integral,
    map_gauss_coef_pares,
    print_node_coef_pares,
    compute_gauss_node_coef_pares,
)


def find_gauss_integral_with_composite_quad_formula(
    node_coef_pares: Tuple[float, float],
    num_nodes: int,
    func: Callable,
    num_div_intervals: int,
    len_div_interval: float,
    seg_beg: float,
) -> float:
    """
    Calculates integral using composite quadratic formala for gauss's
    node-coefficients pares
    """
    integral_sum = 0

    for i in range(num_div_intervals):
        new_seg_beg = seg_beg + i * len_div_interval
        new_seg_end = new_seg_beg + len_div_interval

        for j in range(num_nodes):
            root = node_coef_pares[j][0]
            coef = node_coef_pares[j][1]
            new_root = len_div_interval * root + (new_seg_beg + new_seg_end) / 2

            integral_sum += coef * (func(new_root))

    return integral_sum * len_div_interval / 2


def calc_moments_of_weight_func(
    num_nodes: int, weight_func: Callable, seg_a: float, seg_b: float
) -> Sequence[float]:
    moments = []

    for i in range(2 * num_nodes):
        moment, _ = integrate.quad(lambda x: weight_func(x) * x ** i, seg_a, seg_b)
        moments.append(moment)

    return moments


def solve_eq_system(num_nodes: int, moments_wf: Sequence[float]) -> Sequence[float]:
    rows_list = []
    b_col = np.empty(shape=num_nodes)

    for i in range(num_nodes):
        row = [moments_wf[i + j] for j in range(num_nodes)]
        rows_list.append(row)

        b_col[i] = -moments_wf[num_nodes + i]

    return (np.linalg.inv(rows_list).dot(b_col)).tolist()


def solve_orth_pol_eq_0(
    coef_list: Sequence[float], seg_a: float, seg_b: float
) -> Sequence[float]:
    nodes = []
    segments = root_finding(_pol(coef_list), seg_a, seg_b)

    for seg in segments:
        node = secant(
            _pol(coef_list),
            seg[0],
            seg[1],
            10 ** (-12),
        )
        nodes.append(node)

    return nodes


def _pol(coef_list) -> float:
    return lambda x: sum(
        coef * x ** deg for deg, coef in enumerate(coef_list)
    ) + x ** len(coef_list)


def find_coefs(
    num_nodes: int, roots: Sequence[float], moments_qf: Sequence[float]
) -> Sequence[float]:

    rows_list = []

    for deg in range(num_nodes):
        row = [root ** deg for root in roots]
        rows_list.append(row)

    return (np.linalg.inv(rows_list).dot(moments_qf[:num_nodes])).tolist()


def do_task_1() -> None:
    """
    I/O interactions for task 1
    """
    print("\nЗадача 1.")
    seg_beg = input("\nВведите начало отрезка (default=0): ")
    seg_end = input("Введите конец отрезка (default=1): ")

    seg_beg = EXAMPLE_SEG_BEG if seg_beg == "" else float(seg_beg)
    seg_end = EXAMPLE_SEG_END if seg_end == "" else float(seg_end)

    num_nodes = int(input("Введите число узлов (N): "))
    num_div_intervals = int(input("Введите число промежутков деления (m): "))
    len_div_interval = (seg_end - seg_beg) / num_div_intervals

    node_coef_pares = compute_gauss_node_coef_pares(num_nodes)
    mapped_node_coef_pares = map_gauss_coef_pares(node_coef_pares, seg_beg, seg_end)

    appr_sol = find_gauss_integral_with_composite_quad_formula(
        node_coef_pares,
        num_nodes,
        func_x_weight,
        num_div_intervals,
        len_div_interval,
        seg_beg,
    )

    print("\nУзлы и коэфициенты на [-1, 1]:")
    print(f"Checksum: {print_node_coef_pares(node_coef_pares):.12f}")
    print(f"\nУзлы и коэфициенты на [{seg_beg}, {seg_end}]:")
    print(f"Checksum: {print_node_coef_pares(mapped_node_coef_pares):.12f}")
    print(f"\nApproximate solution: {appr_sol:.12f}")


def do_task_2() -> None:
    """
    I/O interactions for task 2
    """
    print("\nЗадача 2.")

    seg_beg = input("\nВведите начало отрезка (default=0): ")
    seg_end = input("Введите конец отрезка (default=1): ")

    seg_beg = EXAMPLE_SEG_BEG if seg_beg == "" else float(seg_beg)
    seg_end = EXAMPLE_SEG_END if seg_end == "" else float(seg_end)

    num_nodes = int(input("Введите число узлов (N): "))

    print("\nМоменты весовой функции:")
    moments_wf = calc_moments_of_weight_func(num_nodes, weight, seg_beg, seg_end)
    for i, moment in enumerate(moments_wf):
        print(f"{i}) {moment:.12f}")

    print("\nКоэффициенты a_n, ..., a_1 многочлена w(x):")
    a_coefs = solve_eq_system(num_nodes, moments_wf)
    for i, a in enumerate(a_coefs):
        print(f"{num_nodes - i}) {a:.12f}")

    print("\nОртогональный многочлен w(x): ", end="")
    for i, a in enumerate(a_coefs):
        print(f"{a}*x^{i} + ", end="")
    print(f"x^{num_nodes}")

    print("\nУзлы (корни) w(x):")
    roots = solve_orth_pol_eq_0(a_coefs, seg_beg, seg_end)
    for i, root in enumerate(roots):
        print(f"{i}) {root}")

    print("\nКоэффициенты A_1, ... A_n:")
    coefs = find_coefs(num_nodes, roots, moments_wf)
    for i, coef in enumerate(coefs):
        print(f"{i}) {coef}")

    node_coef_pares = [(root, coef) for root, coef in zip(roots, coefs)]
    approx_integral = find_gauss_integral(node_coef_pares, func)
    print(f"\nApproximate integral: {approx_integral:.12f}")

    test_pol = lambda x: x ** (2 * num_nodes - 1)
    node_coef_pares = [(root, coef) for root, coef in zip(roots, coefs)]
    approx_integral = find_gauss_integral(node_coef_pares, test_pol)

    print(f"\nНа одночлене степени N - 1 = {num_nodes - 1}:")
    print(f"2N - 1 момент функции: {moments_wf[2*num_nodes - 1]:.12f}")
    print(f"Approximate integral: {approx_integral:.12f}")


if __name__ == "__main__":
    print(
        """
    Задание 6. Приближенное вычисление интегралов при помощи КФ НАСТ.
    """
    )
    SEG_DEFAULT_BEG: Final = -1
    SEG_DEFAULT_END: Final = 1

    # Variant 6
    func = lambda x: math.sin(x)
    weight = lambda x: -x * math.log(x) if x > 0 else 0
    func_x_weight = lambda x: weight(x) * func(x)

    EXAMPLE_SEG_BEG = 0
    EXAMPLE_SEG_END = 1

    while True:
        print(f"func  : sin(x)")
        print(f"weight: - x * ln(x)")
        print()

        do_task_1()

        do_task_2()

        if input("\nПродолжить с новыми узлами, a, b? (y, n): ") == "y":
            continue
        break
