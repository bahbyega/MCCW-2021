"""
This module contains solutions for sixth laboratory task,
approximate evaluation of integrals with composite quadratic formulas
for Methods of Computation and Computational Workshop course.
"""

import math
from typing import Callable, Tuple, Final
from hw5.solution5 import (
    print_node_coef_pares,
    compute_gauss_node_coef_pares,
)


def find_gauss_integral_with_composite_quad_formula(
    node_coef_pares: Tuple[float, float],
    num_nodes: int,
    func: Callable,
    len_div_interval: float,
    seg_beg: float,
) -> float:
    """
    Calculates integral using composite quadratic formala for gauss's
    node-coefficients pares
    """
    integral_sum = 0

    for i in range(num_div_intevals):
        new_seg_beg = seg_beg + i * len_div_interval
        new_seg_end = new_seg_beg + len_div_interval
        for j in range(num_nodes):
            root = node_coef_pares[j][0]
            coef = node_coef_pares[j][1]
            new_root = len_div_interval * root + (new_seg_beg + new_seg_end) / 2

            integral_sum += coef * (func(new_root))

    return integral_sum * len_div_interval / 2


def do_task_1() -> None:
    """
    Output for task 1
    """
    node_coef_pares = compute_gauss_node_coef_pares(num_nodes)
    appr_sol = find_gauss_integral_with_composite_quad_formula(
        node_coef_pares, num_nodes, func_x_weight, len_div_interval, seg_beg
    )

    print(f"Checksum: {print_node_coef_pares(node_coef_pares):.12f}")
    print(f"\nApproximate solution: {appr_sol}")


def do_task_2() -> None:
    """pass"""
    pass


if __name__ == "__main__":
    SEG_DEFAULT_BEG: Final = -1
    SEG_DEFAULT_END: Final = 1

    # Variant 6
    func = lambda x: math.sin(x)
    weight = lambda x: -x * math.log(x)
    func_x_weight = lambda x: weight(x) * func(x)

    EXAMPLE_SEG_BEG = 0
    EXAMPLE_SEG_END = 1

    seg_beg = input("Введите начало отрезка (default=0): ")
    seg_end = input("Введите конец отрезка (default=1): ")

    seg_beg = EXAMPLE_SEG_BEG if seg_beg == "" else float(seg_beg)
    seg_end = EXAMPLE_SEG_END if seg_end == "" else float(seg_end)

    num_nodes = int(input("Введите число узлов (N): "))
    num_div_intevals = int(input("Введите число промежутков деления (m): "))
    len_div_interval = (seg_end - seg_beg) / num_div_intevals

    while True:
        do_task_1()

        # do_task_2()

        if input("\nПродолжить с новыми узлами, a, b? (y, n): ") == "y":
            continue
        break
