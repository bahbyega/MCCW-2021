import math

from typing import Final, Sequence, Callable, Tuple
from scipy.special import legendre

from hw1.solution1 import root_finding, secant

SEG_DEFAULT_BEG: Final = -1
SEG_DEFAULT_END: Final = 1


def _legendre(x, n) -> float:
    if n == 0:
        return 1
    if n == 1:
        return x
    return (2 * n - 1) / n * _legendre(x, n - 1) * x - (n - 1) / n * _legendre(x, n - 2)


def print_root_coef_pares(nodes_list: Sequence[int]) -> None:
    epsilon = 10 ** (-21)

    for node in nodes_list:
        print(f"\nN = {node}")

        root_coef_pares = find_gauss_integral(node, _legendre, epsilon)

        for root, coef in root_coef_pares:
            print(f"{root:.12f} <-> {coef:.12f}")


def find_gauss_integral(
    node: int, func: Callable, epsilon: float = 10 ** (-12)
) -> Sequence[Tuple[float, float]]:
    root_coef_pares = []

    segments = root_finding(lambda x: func(x, node), SEG_DEFAULT_BEG, SEG_DEFAULT_END)

    for seg in segments:
        root = secant(lambda x: func(x, node), seg[0], seg[1], epsilon)
        coef = (2 * (1 - root ** 2)) / (node ** 2 * func(root, node - 1) ** 2)

        root_coef_pares.append((root, coef))

    return root_coef_pares


def find_meler_integral(node: int, func: Callable) -> float:
    integral_sum = 0

    for k in range(1, node + 1):
        root = math.cos((2 * k - 1) * math.pi / (2 * node))
        print(f"{root} <-> {math.pi / node}")
        integral_sum += func(root)

    return math.pi * integral_sum / node


def calc_integral_sum_from_roots_and_coeffs(
    root_coef_pares: Tuple[float, float], func: Callable, seg_a: float, seg_b: float
) -> float:
    """
    Linearly maps root_coef_paresn on [seg_a, seg_b] to [-1, 1],
    Calculates mapped integral sum
    """
    similarity_coefficient = (seg_b - seg_a) / (SEG_DEFAULT_END - SEG_DEFAULT_BEG)

    mapped_integral_sum = 0
    for root, coef in root_coef_pares:
        mapped_coef = coef * similarity_coefficient
        mapped_root = seg_a + similarity_coefficient * (root - SEG_DEFAULT_BEG)
        mapped_integral_sum += mapped_coef * func(mapped_root)

    return mapped_integral_sum


def do_accuracy_check_for(nodes_list: Sequence[int]) -> None:
    for node in nodes_list:
        if node == 3:
            f = lambda x: 6 * x ** 5 + 2 * x + 34

            exact_integral = 68
            root_coef_pares = find_gauss_integral(node, _legendre)
            approx_integral = calc_integral_sum_from_roots_and_coeffs(
                root_coef_pares, f, SEG_DEFAULT_BEG, SEG_DEFAULT_END
            )

            print(f"\nPolynom: 6x^5 + 2x + 34, exact integral = {exact_integral}")
            print(f"Approximate integral: {approx_integral:.12f}")

        elif node == 4:
            f = lambda x: 8 * x ** 7 + 3 * x ** 2

            exact_integral = 2
            root_coef_pares = find_gauss_integral(node, _legendre)
            approx_integral = calc_integral_sum_from_roots_and_coeffs(
                root_coef_pares, f, SEG_DEFAULT_BEG, SEG_DEFAULT_END
            )

            print(f"\nPolynom: 8x^7 + 3x^2, exact integral = {exact_integral}")
            print(f"Approximate integral: {approx_integral:.12f}")

        else:
            f = lambda x: 10 * x ** 9 + 5 * x ** 4 + 1

            exact_integral = 4
            root_coef_pares = find_gauss_integral(node, _legendre)
            approx_integral = calc_integral_sum_from_roots_and_coeffs(
                root_coef_pares, f, SEG_DEFAULT_BEG, SEG_DEFAULT_END
            )

            print(f"\nPolynom: 10x^9 + 5x^4 + 1, exact integral = {exact_integral}")
            print(f"Approximate integral: {approx_integral:.12f}")


if __name__ == "__main__":
    print(
        """
    Задание 5. КФ Гаусса, ее узлы и коэффициенты.
    Вычисление интегралов при помощи КФ Гаусса.

    КФ Мелера, ее узлы и коэффициенты.
    Вычисление интегралов при помощи КФ Мелера.
    """
    )

    # Variant 6
    example_gauss = lambda x: x * math.log(1 + x)
    example_gauss_seg_beg = 0
    example_gauss_seg_end = 1
    example_gauss_exact_integral = 0.25

    epsilon = 10 ** (-12)

    example_meler = lambda x: math.cos(x) ** 2

    while True:
        print("\n-------------------------------\nЗадание 1.")
        print_root_coef_pares(list(range(1, 9)))

        print("\n-------------------------------\nЗадание 2.")
        nodes_list = [3, 4, 5]
        print(f"Узлы: {nodes_list}")
        do_accuracy_check_for(nodes_list)

        print("\n-------------------------------\nЗадание 3.")

        nodes_list = list(
            map(int, input("Введите список числа узлов (до 4 узлов): ").strip().split())
        )[:4]

        a = input("Введите a (default = 0):")
        b = input("Введите b (default = 1):")

        a = example_gauss_seg_beg if a == "" else float(a)
        b = example_gauss_seg_end if b == "" else float(b)

        for node in nodes_list:
            root_coef_pares = find_gauss_integral(node, _legendre, epsilon)

            approx_integral = calc_integral_sum_from_roots_and_coeffs(
                root_coef_pares,
                example_gauss,
                a,
                b,
            )

            print(f"\nN = {node}")
            for root, coef in root_coef_pares:
                print(f"{root:.12f} <-> {coef:.12f}")

            print(f"\nExact integral = {example_gauss_exact_integral}")
            print(f"Approximate integral = {approx_integral}")
            print(f"Error: {abs(approx_integral - example_gauss_exact_integral)}")

        print("\n-------------------------------\nЗадание 4.")

        nodes_list = list(
            map(int, input("Введите список числа узлов (до 3 узлов): ").strip().split())
        )[:3]

        for node in nodes_list:
            print(f"\nN = {node}")
            approx_integral = find_meler_integral(node, example_meler)
            print(f"Approximate integral = {approx_integral}")

        if input("\nПродолжить с новыми nodes, a, b? (y, n): ") == "y":
            continue
        break
