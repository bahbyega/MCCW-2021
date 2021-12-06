"""
This module contains solutions for fifth laboratory task of
Methods of Computation and Computational Workshop course.
"""

import math
from typing import Callable, Final, Sequence, Tuple

from hw1.solution1 import root_finding, secant


def _legendre(x_val, num_nodes) -> float:
    """
    Calculates legendre function using recurrent formula for homogeneous polynomials
    """
    if num_nodes == 0:
        return 1
    if num_nodes == 1:
        return x_val

    return (2 * num_nodes - 1) / num_nodes * _legendre(x_val, num_nodes - 1) * x_val - (
        num_nodes - 1
    ) / num_nodes * _legendre(x_val, num_nodes - 2)


def print_node_coef_pares(pares: Tuple[float, float]) -> float:
    """
    Prints the pares of nodes and values in "node <-> value" format.
    Returns checksum of coefficients.
    """
    checksum = 0.0

    for node, coef in pares:
        print(f"{node:.12f} <-> {coef:.12f}")
        checksum += coef

    return checksum


def compute_gauss_node_coef_pares(num_nodes: int) -> Sequence[Tuple[float, float]]:
    """
    Calculate pares of nodes and coefficients with gauss' quadratic formula
    """
    node_coef_pares = []

    segments = root_finding(
        lambda x: _legendre(x, num_nodes), SEG_DEFAULT_BEG, SEG_DEFAULT_END
    )

    for seg in segments:
        node = secant(
            lambda x: _legendre(x, num_nodes), seg[0], seg[1], EPS_DEFAULT_VAL
        )
        coef = (2 * (1 - node ** 2)) / (
            num_nodes ** 2 * _legendre(node, num_nodes - 1) ** 2
        )

        node_coef_pares.append((node, coef))

    return node_coef_pares


def compute_meler_node_coef_pares(num_nodes: int) -> Sequence[Tuple[float, float]]:
    """
    Calculate pares of nodes and coefficients with meler's quadratic formula
    """
    node_coef_pares = []

    coef = math.pi / num_nodes
    for num in range(1, num_nodes + 1):
        node = math.cos((2 * num - 1) * math.pi / (2 * num_nodes))
        node_coef_pares.append((node, coef))

    return node_coef_pares


def find_gauss_integral(
    node_coef_pares: Tuple[float, float], func: Callable, seg_a: float, seg_b: float
) -> float:
    """
    Linearly maps node_coef_paresn on [seg_a, seg_b] to [-1, 1],

    Given the pares of nodes and coefficient finds the approximate value of
    an integral on [-1, 1] of function in func.
    """
    similarity_coefficient = (seg_b - seg_a) / (SEG_DEFAULT_END - SEG_DEFAULT_BEG)

    mapped_integral_sum = 0
    for root, coef in node_coef_pares:
        mapped_coef = coef * similarity_coefficient
        mapped_root = seg_a + similarity_coefficient * (root - SEG_DEFAULT_BEG)
        mapped_integral_sum += mapped_coef * func(mapped_root)

    return mapped_integral_sum


def find_meler_integral(node_coef_pares: Tuple[float, float], func: Callable) -> float:
    """
    Given the pares of nodes and coefficient finds the approximate value of
    an integral on [-1, 1] of function in func.
    """
    integral_sum = 0
    num_nodes = len(node_coef_pares)

    for i in range(num_nodes):
        node = node_coef_pares[i][0]
        integral_sum += func(node)

    return math.pi * integral_sum / num_nodes


def do_accuracy_check_for(nodes_list: Sequence[int]) -> None:
    """
    Basically a test (or either 3 tests) of compute compute_gauss_node_coef_pares()
    on three different numbers of nodes for interpolated quadratic formula
    """
    for num_nodes in nodes_list:
        if num_nodes == 3:
            func = lambda x: 6 * x ** 5 + 2 * x + 34

            exact_integral = 68
            node_coef_pares = compute_gauss_node_coef_pares(num_nodes)
            approx_integral = find_gauss_integral(
                node_coef_pares, func, SEG_DEFAULT_BEG, SEG_DEFAULT_END
            )

            print(f"\nPolynom: 6x^5 + 2x + 34, exact integral = {exact_integral}")
            print(f"Approximate integral: {approx_integral:.12f}")

        elif num_nodes == 4:
            func = lambda x: 8 * x ** 7 + 3 * x ** 2

            exact_integral = 2
            node_coef_pares = compute_gauss_node_coef_pares(num_nodes)
            approx_integral = find_gauss_integral(
                node_coef_pares, func, SEG_DEFAULT_BEG, SEG_DEFAULT_END
            )

            print(f"\nPolynom: 8x^7 + 3x^2, exact integral = {exact_integral}")
            print(f"Approximate integral: {approx_integral:.12f}")

        else:
            func = lambda x: 10 * x ** 9 + 5 * x ** 4 + 1

            exact_integral = 4
            node_coef_pares = compute_gauss_node_coef_pares(num_nodes)
            approx_integral = find_gauss_integral(
                node_coef_pares, func, SEG_DEFAULT_BEG, SEG_DEFAULT_END
            )

            print(f"\nPolynom: 10x^9 + 5x^4 + 1, exact integral = {exact_integral}")
            print(f"Approximate integral: {approx_integral:.12f}")


def do_task_1() -> None:
    """
    Runs task 1 output
    """
    print("\n-------------------------------\nЗадание 1.")
    num_nodes_list = list(range(1, 9))

    for num_nodes in num_nodes_list:
        node_coef_pares = compute_gauss_node_coef_pares(num_nodes)
        print(f"\nЧисло узлов: {num_nodes}")
        print(f"Checksum: {print_node_coef_pares(node_coef_pares):.12f}")


def do_task_2() -> None:
    """
    Runs task 2 output
    """
    print("\n-------------------------------\nЗадание 2.")

    num_nodes_list = [3, 4, 5]
    print(f"Узлы: {num_nodes_list}")
    do_accuracy_check_for(num_nodes_list)


def do_task_3() -> None:
    """
    Runs task 3 output
    """
    print("\n-------------------------------\nЗадание 3.")

    num_nodes_list = list(
        map(int, input("Введите список числа узлов (до 4 узлов): ").strip().split())
    )[:4]

    seg_a = input("Введите a (default = 0):")
    seg_b = input("Введите b (default = 1):")

    seg_a = EXAMPLE_GAUSS_SEG_BEG if seg_a == "" else float(seg_a)
    seg_b = EXAMPLE_GAUSS_SEG_END if seg_b == "" else float(seg_b)

    for num_nodes in num_nodes_list:
        node_coef_pares = compute_gauss_node_coef_pares(num_nodes)

        approx_integral = find_gauss_integral(
            node_coef_pares,
            example_gauss,
            seg_a,
            seg_b,
        )

        print(f"\nN = {num_nodes}")
        print(f"Checksum: {print_node_coef_pares(node_coef_pares):.12f}")
        print(f"\nExact integral = {EXAMPLE_GAUSS_EXACT_INTEGRAL}")
        print(f"Approximate integral = {approx_integral}")
        print(f"Error: {abs(approx_integral - EXAMPLE_GAUSS_EXACT_INTEGRAL)}")


def do_task_4() -> None:
    """
    Runs task 4 output
    """
    print("\n-------------------------------\nЗадание 4.")

    num_nodes_list = list(
        map(int, input("Введите список числа узлов (до 3 узлов): ").strip().split())
    )[:3]

    for num_nodes in num_nodes_list:
        print(f"\nN = {num_nodes}")
        node_coef_pares = compute_meler_node_coef_pares(num_nodes)
        approx_integral = find_meler_integral(node_coef_pares, example_meler)

        print(f"Checksum: {print_node_coef_pares(node_coef_pares):.12f}")
        print(f"\nApproximate integral = {approx_integral}")


if __name__ == "__main__":

    print(
        """
    Задание 5. КФ Гаусса, ее узлы и коэффициенты.
    Вычисление интегралов при помощи КФ Гаусса.

    КФ Мелера, ее узлы и коэффициенты.
    Вычисление интегралов при помощи КФ Мелера.
    """
    )

    SEG_DEFAULT_BEG: Final = -1
    SEG_DEFAULT_END: Final = 1
    EPS_DEFAULT_VAL: Final = 10 ** (-20)

    # Variant 6
    example_gauss = lambda x: x * math.log(1 + x)
    example_meler = lambda x: math.cos(x) ** 2

    EXAMPLE_GAUSS_SEG_BEG = 0
    EXAMPLE_GAUSS_SEG_END = 1
    EXAMPLE_GAUSS_EXACT_INTEGRAL = 0.25

    while True:
        do_task_1()

        do_task_2()

        do_task_3()

        do_task_4()

        if input("\nПродолжить с новыми узлами, a, b? (y, n): ") == "y":
            continue
        break
