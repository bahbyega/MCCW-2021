import math


def print_result(method_name, start_value, step_count,
                 approx_sol, error, residual):
    print(f'''      {method_name}
           Начальное приближение: {start_value:.8f}
           Количество шагов: {step_count}
           Приближенное решение: {approx_sol:.8f}
           |x_m - x_m-1|: {error:.8f}
           Абс. величина невязки: {residual}\n''')


def root_finding(f, A, B, N=1000):
    segments = list()

    H = (B - A) / N
    counter = 0
    X1 = A
    X2 = X1 + H
    Y1 = f(X1)

    while X2 <= B:
        Y2 = f(X2)

        if Y1 * Y2 <= 0:
            counter += 1
            print(f' [{X1:.8f}, {X2:.8f}]')
            segments.append((X1, X2))

        X1 = X2
        X2 = X1 + H
        Y1 = Y2

    print(f'Число корней на промежутке: {counter} ')

    return segments


def bisection(f, A, B, epsilon):
    method_name = "Метод бисекции"
    step_counter = 0
    start_value = A

    last_midpoint = 0
    midpoint = (A + B) / 2

    while B - A > 2 * epsilon:
        last_midpoint = midpoint
        midpoint = (A + B) / 2

        if f(A) * f(midpoint) <= 0:
            B = midpoint
        else:
            A = midpoint

        step_counter += 1

    print_result(method_name,
                 start_value,
                 step_counter,
                 midpoint,
                 abs(midpoint - last_midpoint),
                 abs(f(midpoint)))


def newton(f, Df, A, B, epsilon, max_iter=200, modified=False):
    x_prev = A
    x_sol = x_prev

    if not modified:
        method_name = "Метод Ньютона"
    else:
        method_name = "Модифицированный метод Ньютона"

    for step in range(max_iter):
        if abs(f(x_sol)) < epsilon:
            print_result(method_name,
                         A,
                         step,
                         x_sol,
                         abs(x_sol - x_prev),
                         abs(f(x_sol)))
            return x_sol

        if Df(x_sol) == 0:
            print("No solution")
            return None

        x_prev = x_sol

        if not modified:
            x_sol = x_sol - f(x_sol) / Df(x_sol)
        else:
            x_sol = x_sol - f(x_sol) / Df((A + B) / 2)

    print(f'{method_name} exceeded iterations limit')


def secant(f, A, B, epsilon, max_iter=200):
    method_name = "Метод секущих"

    a_init = A
    b_init = B

    fa = f(a_init)
    fb = f(b_init)

    for step in range(max_iter):
        x_sol = (A * fb - B * fa) / (fb - fa)

        if abs(fb) < epsilon:
            print_result(method_name,
                         a_init,
                         step,
                         A,
                         abs(B - x_sol),
                         abs(f(B)))
            return A

        A, B = B, x_sol
        fa, fb = fb, f(x_sol)

    print(f'{method_name} exceeded iterations limit')


if __name__ == '__main__':
    A = -9
    B = 1
    epsilon = 10**(-7)
    N = 1000

    def f(x): return 8 * math.cos(x) - x - 6
    def Df(x): return -8 * math.sin(x) - 1

    print(f'''Задание 1. ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ
          ----------------------------------------------
          Исходные параметры: A = {A}, B = {B},
          Исходное уравнение: f(x) = 8 * cos(x) - x - 6,
                              epsilon = {epsilon:.7f}''')

    print(f'\n------------- Отделение корней -------------')
    segments = root_finding(f, A, B, N)

    print(f'\n------------- Уточнение корней -------------')
    for segment in segments:
        a_k, b_k = segment[0], segment[1]

        print(f'На отрезке [{a_k:.8f}, {b_k:.8f}]:')

        bisection(f, a_k, b_k, epsilon)
        newton(f, Df, a_k, b_k, epsilon, 20)
        newton(f, Df, a_k, b_k, epsilon, 20, True)
        secant(f, a_k, b_k, epsilon)

        print(f'---')
