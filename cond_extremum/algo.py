import numpy as np
from scipy import optimize
from scipy.optimize import fsolve


def penalty(f, constraint_fs, xs, r, c, ε):
    it = 0
    while True:
        it += 1
        constr_sqr_sum = sum(h_i(xs) for h_i in constraint_fs)
        penalty_func = lambda xs: r * constr_sqr_sum / 2
        new_xs = optimize.minimize(lambda x: f(x) + penalty_func(xs), xs).x
        penalty_val = penalty_func(new_xs)
        if abs(penalty_val) < ε:
            return new_xs, f(new_xs)
        r *= c
        xs = new_xs


def barier(f, constraint_fs, xs, r, c, ε):
    it = 0
    while True:
        it += 1
        penalty_func = lambda x, r: -(r ** it)
        p = penalty_func

        constr_sqr_sum = sum(h_i(xs) for h_i in constraint_fs)
        if constr_sqr_sum != 0:
            p = lambda xs, r: penalty_func(xs, r) / constr_sqr_sum

        new_xs = optimize.minimize(lambda x: f(x) + p(xs, r), xs).x
        penalty_value = p(new_xs, r)
        if abs(penalty_value) < ε:
            return new_xs, f(new_xs)

        xs = new_xs
        r /= c


def lagrange_functions(f, constraint_fs, xs, ε, inc_param, r, lambdas, μ):
    ε /= 100
    it = 0
    x_c = xs
    while True:
        it += 1

        def func(x):
            return f(x) + langrange_lambda + eq_penalty(x) + 1 / (2 * r) * sum(
                np.array(neq_penalty) - np.array(μ_squared))

        langrange_lambda = np.sum(np.matmul(np.array(lambdas), [h_i(xs) for h_i in constraint_fs]))
        constr_sqr_sum = sum(h_i(xs) for h_i in constraint_fs)

        eq_penalty = lambda x: (r / 2) * constr_sqr_sum
        neq_penalty = [
            max(0.0, μ[0] + r * pow(h_i(x_c), 2)) for h_i in constraint_fs
        ]
        μ_squared = [μ[0] ** 2, μ[1] ** 2, μ[2] ** 2]

        x_new = optimize.minimize(func, x_c).x

        new_penalty = [
            max(0.0, μ[0] + r * pow(h_i(x_new), 2)) for h_i in constraint_fs
        ]
        new_μ_squared = [μ[0] ** 2, μ[1] ** 2, μ[2] ** 2]

        constr_sqr_sum = sum(h_i(x_new) for h_i in constraint_fs)
        penalty_value = r / 2 * constr_sqr_sum * \
                        sum(np.array(new_penalty) - np.array(new_μ_squared))

        if abs(penalty_value) < ε:
            return x_new, f(x_new)

        x_c = x_new
        r *= inc_param
        lambdas += np.array(np.multiply(r, [h_i(xs) for h_i in constraint_fs]))


def projection_gradient(f, grad_f, xs, constraint_fs, deriv_contrains, ε1, ε2=None, MAX_ITERATIONS=10000):
    if ε2 is None:
        ε2 = ε1

    initial_constr_n = len(constraint_fs)
    start_xs = xs

    if not all(constr_f(xs) <= 0 for constr_f in constraint_fs):
        raise Exception('xStart not in domain')

    for it in range(MAX_ITERATIONS):
        grad_val = grad_f(xs)
        grad_matrix = np.array([grad_val]).transpose()
        border_vals = [constr_f(xs) for constr_f in constraint_fs]
        passive_constraints = {
            i for i, bv in enumerate(border_vals)
            if not (ε1 <= bv <= 0)
        }

        if len(passive_constraints) != initial_constr_n:
            matr = np.array([
                [f_i(xs) for f_i in row]
                for i, row in enumerate(deriv_contrains)
                if deriv_contrains not in passive_constraints
            ])
            m = matr.transpose() @ np.inv(matr @ matr.transpose()) @ matr
            p = np.eye(len(m)) - m
            dx = (-p @ grad_matrix)[:, 0]
        else:
            dx = (-grad_matrix)[:, 0]

        distamce_to_borders = []
        for i, constr_f in enumerate(constraint_fs):
            if i in passive_constraints:
                continue

            zero_val, _, check, _ = fsolve(
                lambda a: constr_f(start_xs + a * dx),
                np.eye(1),
            )
            if zero_val >= 0 and check == 1:
                distamce_to_borders.append(zero_val[0])

        arg_min = optimize.minimize_scalar(lambda a: f(xs + a * dx)).x

        arg = min(arg_min, distamce_to_borders)
        dx *= arg
        xs += dx
        if np.linalg.norm(dx) < ε2:
            return xs, f(xs)

    print('Reached the end before convergence')
    return xs, f(xs)
