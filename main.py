
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
    WHITE = '\033[38;5;15m'

def neville(x_data, y_data, x_interpolate):
    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have the same length.")
    if len(set(x_data)) != len(x_data):
        raise ValueError("x_data must contain unique values.")

    n = len(x_data)

    tableau = [[0.0] * n for _ in range(n)]

    for i in range(n):
        tableau[i][0] = y_data[i]

    for j in range(1, n):
        for i in range(n - j):
            denominator = x_data[i] - x_data[i + j]
            if denominator == 0:
                raise ZeroDivisionError(f"Duplicate x_data values found: x_data[{i}] = x_data[{i + j}]")

            tableau[i][j] = ((x_interpolate - x_data[i + j]) * tableau[i][j - 1] -
                             (x_interpolate - x_data[i]) * tableau[i + 1][j - 1]) / denominator

    return tableau[0][n - 1]

def spline_interpolation(x_data, y_data, x_interpolate):
    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have the same length.")
    if len(set(x_data)) != len(x_data):
        raise ValueError("x_data must contain unique values.")

    n = len(x_data)

    a = list(y_data)
    b = [0.0] + [(a[i] - a[i - 1]) / (x_data[i] - x_data[i - 1]) for i in range(1, n)]
    c = [0.0] * n
    d = [0.0] * (n - 1)

    for i in range(1, n - 1):
        d[i] = 6 * ((a[i + 1] - a[i]) / (x_data[i + 1] - x_data[i]) -
                    (a[i] - a[i - 1]) / (x_data[i] - x_data[i - 1]))

    l, mu, z = [0.0] * n, [0.0] * n, [0.0] * n
    l[0] = 1
    mu[0] = 0
    z[0] = 0

    for i in range(1, n - 1):
        l[i] = 2 * (x_data[i + 1] - x_data[i - 1]) - (x_data[i] - x_data[i - 1]) * mu[i - 1]
        mu[i] = (x_data[i + 1] - x_data[i]) / l[i]
        z[i] = (d[i] - (x_data[i] - x_data[i - 1]) * z[i - 1]) / l[i]

    l[n - 1] = 1
    z[n - 1] = 0
    c[n - 1] = 0

    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / (x_data[j + 1] - x_data[j]) - (x_data[j + 1] - x_data[j]) * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * (x_data[j + 1] - x_data[j]))

    # Step 3: Evaluate the spline at x_interpolate
    for i in range(n - 1):
        if x_data[i] <= x_interpolate <= x_data[i + 1]:
            dx = x_interpolate - x_data[i]
            return a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3

    raise ValueError("x_interpolate is out of the bounds of x_data.")

def lagrange_interpolation(x_data, y_data, x):
    n = len(x_data)
    result = 0.0

    x_data = np.array(x_data)
    y_data = np.array(y_data)

    for i in range(n):
        term = y_data[i]
        term *= np.prod([(x - x_data[j]) / (x_data[i] - x_data[j]) for j in range(n) if j != i])
        result += term

    return result

def linearInterpolation(table_points, point):
    result = 0
    found = False

    for i in range(len(table_points) - 1):
        x1, y1 = table_points[i]
        x2, y2 = table_points[i + 1]

        if x1 <= point <= x2:
            result = (((y2 - y1) / (x2 - x1)) * (point - x1)) + y1
            print(bcolors.WHITE + f" Linear Interpolation Result: " + bcolors.ENDC + f"{round(result, 4)}")
            found = True
            break

    if not found:
        if point < table_points[0][0]:
            x1, y1 = table_points[0]
            x2, y2 = table_points[1]
        else:
            x1, y1 = table_points[-2]
            x2, y2 = table_points[-1]

        m = (y2 - y1) / (x2 - x1)
        result = y1 + m * (point - x1)
        print(bcolors.WHITE + f" Linear Interpolation Result: " + bcolors.ENDC + f"{round(result, 4)}")

    return round(result, 4)

def PolynomialInterpolation(table_points, x):
    n = len(table_points)
    matrix = [[point[0] ** i for i in range(n)] for point in table_points]
    b = [point[1] for point in table_points]

    print(bcolors.OKBLUE, "Matrix from the points: " + bcolors.ENDC, '\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: "+ bcolors.ENDC, b, '\n')

    matrix_solution = SolveMatrix(matrix, b)

    result = sum([matrix_solution[i] * (x ** i) for i in range(n)])

    print(bcolors.OKBLUE, "\nThe Polynomial:"+ bcolors.ENDC)
    print("P(x) = " + ' + '.join([f"({matrix_solution[i]}) * x^{i}" for i in range(n)]))
    print(bcolors.WHITE, f"\nThe Result of P(x={x}) is:", bcolors.ENDC, result)


    return result

def RowXchange(matrix, vector):
    n = len(matrix)
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(matrix[r][i]))
        if i != max_row:
            matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
            vector[i], vector[max_row] = vector[max_row], vector[i]
    return matrix, vector

def InverseMatrix(matrix, vector=None):
    return np.linalg.inv(matrix)

def MulMatrixVector(matrix, vector):
    return np.dot(matrix, vector)

def MakeIMatrix(rows, cols):
    return np.eye(rows, cols)

def MultiplyMatrix(A, B):
    return np.dot(A, B)

def Determinant(matrix, _=None):
    return np.linalg.det(matrix)

def Cond(matrix, inverse_matrix=None):
    if inverse_matrix is None:
        inverse_matrix = np.linalg.inv(matrix)
    return np.linalg.norm(matrix) * np.linalg.norm(inverse_matrix)

def GaussJordanElimination(matrix, vector):
    matrix, vector = RowXchange(matrix, vector)
    invert = np.linalg.inv(matrix)
    return MulMatrixVector(invert, vector)

def LUDecomposition(matrix):
    n = len(matrix)
    U = np.copy(matrix)
    L = MakeIMatrix(n, n)

    for i in range(n):
        for j in range(i + 1, n):
            multiplier = U[j][i] / U[i][i]
            U[j] -= multiplier * U[i]
            L[j][i] = multiplier

    return L, U

def SolveLU(matrix, vector):
    L, U = LUDecomposition(matrix)
    y = np.linalg.solve(L, vector)
    x = np.linalg.solve(U, y)
    return x

def SolveMatrix(matrix, vector):
    detA = Determinant(matrix, 1)
    print(bcolors.WHITE, "DET(A) = ", detA, bcolors.ENDC)

    if detA != 0:
        print(bcolors.WHITE, "Cond(A) = ", Cond(matrix, InverseMatrix(matrix, vector)), bcolors.ENDC)
        result = GaussJordanElimination(matrix, vector)
        return result
    else:
        print("Singular Matrix - Perform LU Decomposition\n")
        L, U = LUDecomposition(matrix)
        print("Matrix L: \n", np.array(L))
        print("\nMatrix U: \n", np.array(U))
        result = MultiplyMatrix(L, U)
        print("\nMatrix A = L * U: \n", np.array(result))
        return result

if __name__ == '__main__':
    table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    x = 1.28

    print(bcolors.PINK, "----------------- Interpolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.PURPLE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.PURPLE, "Finding an approximation for the point: ", bcolors.ENDC, x, '\n')
    print(bcolors.PURPLE, "Using Linear Interpolation:", bcolors.ENDC)
    result_linear = linearInterpolation(table_points, x)

    print(bcolors.PURPLE, "\n Using Lagrange Interpolation: ", bcolors.ENDC)
    x_data, y_data = zip(*table_points)
    result_lagrange = lagrange_interpolation(x_data, y_data, x)
    print(bcolors.WHITE, f"Lagrange Interpolation Result:", bcolors.ENDC, result_lagrange)
    print(bcolors.PURPLE, "\nUsing Polynomial Interpolation:", bcolors.ENDC)
    result_poly = PolynomialInterpolation(table_points, x)

    print(bcolors.PURPLE, "\nUsing Nevil Interpolation: ", bcolors.ENDC)
    interpolated_value = neville(x_data, y_data, x)
    print(bcolors.WHITE, f"Nevil Interpolation Result: ", bcolors.ENDC, interpolated_value)
    print(bcolors.PURPLE, "\nUsing Spline Interpolation: ", bcolors.ENDC)
    interpolated_value = spline_interpolation(x_data, y_data, x)
    print(bcolors.WHITE, f"Spline Interpolation Result: ", bcolors.ENDC, interpolated_value)

    print(bcolors.PINK, "\n------------------------------------------------------", bcolors.ENDC)


