'''
 שירה רחל בורוכוב : 345887046
אביגיל מוסי : 322227711
רותם דינו : 209168442
'''
import math
import numpy as np


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    PINK = '\033[38;5;206m'
    PURPLE = '\033[38;5;135m'




def max_steps(a, b, err):
    s = int(np.floor(- np.log2(err / (b - a)) / np.log2(2) - 1))
    return s

def bisection_method(f, a, b, tol=1e-6):
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a root")
    c, k = 0, 0
    steps = max_steps(a, b, tol)

    print("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}".format("Iteration", "a", "b", "f(a)", "f(b)", "c", "f(c)"))

    while abs(b - a) > tol and k < steps:
        c = a + (b - a) / 2

        if f(c) == 0:
            return c

        if f(c) * f(a) < 0:
            b = c
        else:
            a = c

        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f}".format(k, a, b, f(a), f(b), c, f(c)))
        k += 1

    return c


def find_all_roots(f, intervals, tol=1e-6):
    roots = []
    for interval in intervals:
        try:
            a, b = interval
            root = bisection_method(f, a, b, tol)
            roots.append(root)
        except Exception as e:
            print(bcolors.WARNING, f"Skipping interval {interval}: {e}", bcolors.ENDC)
    return roots


def secant_method(f, x0, x1, TOL, N=50):
    print("{:<10} {:<15} {:<15} {:<15}".format(bcolors.PURPLE, "Iteration: ", "xo", "x1", "p"), bcolors.ENDC)
    for i in range(N):
        if f(x1) - f(x0) == 0:
            print("method cannot continue.")
            return

        p = x0 - f(x0) * ((x1 - x0) / (f(x1) - f(x0)))

        if abs(p - x1) < TOL:
            return p  # Procedure completed successfully
        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f}".format(i, x0, x1, p))
        x0 = x1
        x1 = p
    return p


def find_all_roots1(f, intervals, tol=1e-6):
    roots = []
    for interval in intervals:
        try:
            zip(*intervals)
            x0, x1 = interval
            root = secant_method(f, x0, x1, tol)
            roots.append(root)
        except Exception as e:
            print(bcolors.PURPLE, f"Skipping interval {interval}: {e}", bcolors.ENDC)
    return roots

def newton_raphson(f, df, p0, tol, max_iter=50):
    print(bcolors.OKBLUE,"{:<10} {:<15} {:<15}".format("Iteration", "p0", "p1"),bcolors.ENDC)
    for i in range(max_iter):
        df_p0 = df(p0)
        if df_p0 == 0:
            print("Derivative is zero at p0, method cannot continue.")
            return None

        p = p0 - f(p0) / df_p0

        if abs(p - p0) < tol:
            return p
        print("{:<10} {:<15.9f} {:<15.9f}".format(i, p0, p))
        p0 = p
    return p

def find_all_roots2(f, df, initial_guesses, tol, max_iter):
    roots = []
    for guess in initial_guesses:
        root = newton_raphson(f, df, guess, tol, max_iter)
        if root is not None and all(abs(root - r) > tol for r in roots):
            roots.append(root)
    return roots


if __name__ == '__main__':
    if __name__ == '__main__':
        function = lambda x: x ** 2 - 4 * math.sin(x)
        derivative = lambda x: 2 * x - 4 * math.cos(x)
        initial_guesses = [20, -20, 0, 10, -10]
        tolerance = 1e-6
        max_iterations = 1000

        print(bcolors.PINK, "Using Newton-Raphson method to find roots:", bcolors.ENDC)
        print(bcolors.PINK, "Final result: ",
              find_all_roots2(function, derivative, initial_guesses, tolerance, max_iterations))

        f = lambda x: x ** 2 - 4 * math.sin(x)
        intervals = [(1, 3), (-1, 1), (-3, -1), (3, 5), (-5, -3)]

        print(bcolors.OKGREEN, "Using Bisection method to find roots:", bcolors.ENDC)
        roots = find_all_roots(f, intervals)
        print(bcolors.OKGREEN, f"\nThe equation f(x) has approximate roots at x = {roots}", bcolors.ENDC)

        f = lambda x: x ** 2 - 4 * math.sin(x)
        x0 = 80
        x1 = 100
        TOL = 1e-5
        N = 20
        intervals = [(0, 2), (1, 2)]

        print(bcolors.PURPLE, "Using Secant method to find roots:", bcolors.ENDC)
        roots = find_all_roots1(f, intervals)
        print(bcolors.PINK, f"\nThe equation f(x) has approximate roots at x = {roots}", bcolors.ENDC)
