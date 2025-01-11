from math import e
import sympy as sp
import math


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    GREYBG = '\033[100m'
    REDBG = '\033[101m'
    GREENBG = '\033[102m'
    YELLOWBG = '\033[103m'
    BLUEBG = '\033[104m'
    PINKBG = '\033[105m'
    CYANBG = '\033[106m'


def EvaluateError(startPoint, endPoint):
    exp = pow(10, -10)
    if endPoint - startPoint == 0:
        return 100
    return ((-1) * math.log(exp / (endPoint - startPoint), e)) / math.log(2, e)


def SecantMethod(polynomial, firstGuess, secondGuess, epsilon, iterCounter):
    if iterCounter > 100:
        return None

    if abs(secondGuess - firstGuess) < epsilon:
        print("After ", iterCounter, "iterations, the root found is: ", bcolors.OKBLUE, round(secondGuess, 6), bcolors.ENDC)
        return round(secondGuess, 6)

    next_guess = (firstGuess * polynomial(secondGuess) - secondGuess * polynomial(firstGuess)) / \
                 (polynomial(secondGuess) - polynomial(firstGuess))

    return SecantMethod(polynomial, secondGuess, next_guess, epsilon, iterCounter + 1)


def SecantMethodInRangeIterations(f, check_range, epsilon=0.0001):
    roots = []
    iterCounter = 0
    for i in check_range:
        if i == check_range[-1]:
            break
        for sep in range(0, 10):
            startPoint = round(i + (sep * 0.1), 2)
            endPoint = round(i + ((sep + 1) * 0.1), 2)
            print(bcolors.HEADER, "Checked range:", startPoint, "-", endPoint, bcolors.ENDC)
            local_root = SecantMethod(f, startPoint, endPoint, epsilon, iterCounter)
            if local_root is not None and round(local_root, 6) not in roots:
                roots.append(round(local_root, 6))
            else:
                print(bcolors.FAIL, "Already found that root.", bcolors.ENDC)
    return roots


def NewtonsMethod(f, x0, tries=100, epsilon=0.0001, symbol=sp.symbols('x')):
    if f.subs(symbol, x0) == 0:
        return 0
    for i in range(tries):
        print(bcolors.OKBLUE, "Attempt #", i + 1, ":", bcolors.ENDC)
        print("f({0}) = {1} = {2}".format(x0, f, round(f.subs(symbol, x0), 2)))
        print("f'({0}) = {1} = {2}".format(x0, sp.diff(f, symbol), round(sp.diff(f, symbol).subs(symbol, x0), 2)))
        if sp.diff(f, symbol).subs(symbol, x0) == 0.0:
            continue
        next_x = (x0 - f.subs(symbol, x0) / sp.diff(f, symbol).subs(symbol, x0))

        print("next_X = ", round(next_x, 2))
        t = abs(next_x - x0)
        if t < epsilon:
            print(bcolors.OKGREEN, "Found a Root Solution; X =", round(next_x, 8), bcolors.ENDC)
            return next_x
        x0 = next_x
    print(bcolors.FAIL, "Haven't Found a Root Solution; (returning None)", bcolors.ENDC)
    return None


def NewtonsMethodInRangeIterations(f, check_range, tries=10, epsilon=0.0001, symbol=sp.symbols('x')):
    roots = []
    for i in check_range:
        if i == check_range[-1]:
            break
        for sep in range(0, 10):
            check_number = round(i + (sep * 0.1), 2)
            print(bcolors.HEADER, "First guess:", check_number, bcolors.ENDC)
            local_root = NewtonsMethod(f, check_number, tries, epsilon, symbol)
            if local_root is not None and round(local_root, 6) not in roots:
                roots += [round(local_root, 6)]
            else:
                print(bcolors.FAIL, "Already found that root.", bcolors.ENDC)
    return roots


def BisectionMethod(polynomial, startPoint, endPoint, epsilon, iterCounter):
    roots = []
    middle = (startPoint + endPoint) / 2

    if iterCounter > EvaluateError(startPoint, endPoint):
        print(bcolors.FAIL, "The Method isn't convergent.", bcolors.ENDC)
        return roots

    if (abs(endPoint - startPoint)) < epsilon:
        print("After ", iterCounter, "iterations, the root found is: ", bcolors.OKBLUE, round(middle, 6), bcolors.ENDC)
        roots.append(round(middle, 6))
        return roots

    if polynomial(startPoint) * polynomial(middle) > 0:
        roots += BisectionMethod(polynomial, middle, endPoint, epsilon, iterCounter + 1)
        return roots
    else:
        roots += BisectionMethod(polynomial, startPoint, middle, epsilon, iterCounter + 1)
        return roots


def BisectionMethodSections(polynomial, check_range, epsilon):
    iterCounter = 0
    result = []
    for i in check_range:
        for sep in range(1, 10):
            separate = round(i + (sep * 0.1), 2)
            next_separate = round(separate + 0.1, 2)
            if polynomial(next_separate) == 0:
                print(bcolors.OKBLUE, "Root in ", next_separate, bcolors.ENDC)
                result.append(next_separate)
            if (polynomial(separate) * polynomial(next_separate)) < 0:
                result += BisectionMethod(polynomial, separate, next_separate, epsilon, iterCounter)
    return result


def MainFunction():
    roots = []
    x = sp.symbols('x')

    my_f = 3*x ** 4 + 4*x ** 3 + 2*x

    def Polynomial(X):
        return my_f.subs(x, X)

    my_f_diff = lambda a: sp.diff(my_f, x).subs(x, a)

    checkRange = range(-2, 1)
    epsilon = 0.0001

    print("Finding roots of the equation f(X) = 3X^4 + 4X^3 + 2X\n")
    choice = int(input(
        "Which method do you want? \n\t1. Bisection Method \n\t2. Newton Raphson\n\t3. Secant Method\n"))
    if choice == 1:
        roots += BisectionMethodSections(Polynomial, checkRange, epsilon)
    elif choice == 2:
        roots += NewtonsMethodInRangeIterations(my_f, checkRange, 10, 0.000001)
    elif choice == 3:
        roots += SecantMethodInRangeIterations(Polynomial, checkRange, 0.0000001)
    else:
        print(bcolors.FAIL, "Invalid input", bcolors.ENDC)
        return

    print("\nThere are ", bcolors.OKBLUE, len(roots), "roots: ", roots, bcolors.ENDC)


MainFunction()