from algo import *


def var_1(f):
    a = 50
    b = 2
    f0 = 10
    return f(a, b, f0)


def var_10(f):
    a = 350
    b = 2
    f0 = 110
    return f(a, b, f0)


@var_1
def fun_Rosenbrock(a, b, f0):
    def _fun_Rosenbrock(xs):
        s = 0
        for i in range(len(xs) - 1):
            s += a * (xs[i] ** 2 - xs[i + 1]) ** 2 + b * (xs[i] - 1) ** 2
        return s + f0

    return _fun_Rosenbrock


@var_1
def grad_Rosenbrock(a, b, f0):
    def _grad_Rosenbrock(xs):
        x, y = xs
        return np.array([
            4 * a * x * (x ** 2 - y) + 2 * b * (x - 1),
            -2 * a * (x ** 2 - y),
        ])

    return _grad_Rosenbrock


cond_funcs = [
    lambda xs: xs[0] ** 2 + xs[1] ** 2 - 1,
    lambda xs: -xs[0],
    lambda xs: -xs[1],
]

deriv_cond_funcs = [
    lambda xs: [2 * xs[0], 2 * xs[1]],
    lambda xs: [-1, 0],
    lambda xs: [0, -1],
]

if __name__ == '__main__':
    x_start = [0.279, 0.326]
    epsilon = 0.0001
    lambda_ = 1.0
    betta = 1.5

    f = fun_Rosenbrock
    print(
        'Penalty Method\t',
        *penalty(f, cond_funcs, x_start, 1.0, 0.1, epsilon),
    )
    print(
        'Barier Method\t',
        *barier(f, cond_funcs, x_start, epsilon, 1.0, 0.1),
    )
    print(
        'Lagrange Method\t',
        *lagrange_functions(f, cond_funcs, x_start, epsilon, 0.1, 1.0, [0.01, 0.01, 0.01], [0.01, 0.01, 0.01]),
    )
    print(
        'Gradient Projections Method\t',
        *projection_gradient(f, grad_Rosenbrock, x_start, cond_funcs, deriv_cond_funcs, epsilon),
    )
