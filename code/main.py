import math
import timeit
from scipy import integrate
from scipy.optimize import root_scalar
import math

import scipy.special
from scipy import integrate

# E ~ Bin(N,p)
# Описание модели

# Pb(k | N,p)
def Bin(N, p):
    # k - успешное число испытаний
    # N - количество испытаний
    # p - вероятность успеха

    startTime = timeit.default_timer()
    result = 0

    leftBound = math.trunc(N*p - 1)
    rightBound = math.trunc(N*p + 1)

    for k in range(leftBound, rightBound + 1):
        result += math.factorial(N) / (math.factorial(k) * math.factorial(N-k)) \
                  * math.pow(p, k) * math.pow(1-p, N-k)

    endTime = timeit.default_timer()
    workTime = (endTime - startTime) * 1000

    return result, f"Время работы программы:{workTime:.4f} (мс)"

def TheoryPuasson(N,p):

    startTime = timeit.default_timer()
    result = 0

    leftBound = math.trunc(N * p - 1)
    rightBound = math.trunc(N * p + 1)

    #Pb(m | N,p) -> P(m|Lamba) = (Lamba**m / m!) * e^(-Lamba)
    Lamba = N*p

    for m in range(leftBound, rightBound + 1):
        result += (math.pow(Lamba, m) / math.factorial(m)) * math.exp(-Lamba)

    endTime = timeit.default_timer()
    workTime = (endTime - startTime) * 1000

    return result, f"Время работы программы:{workTime:.4f} (мс)"

# нормальная функция распределения
# Ф(х) ~ F(x) = Intergral(phi(t), -inf, x) , где phi(t) = (1/sqrt(2*pi)) * exp( (-t^2/2) )
def Phi(t):
    return (1 / math.sqrt(2*math.pi)) * math.exp( -math.pow(t,2) / 2 )

def F(x):
    return integrate.quad(Phi, -math.inf, x)

def TheoryML_LocalType(N,p):
    startTime = timeit.default_timer()
    result = 0

    leftBound = math.trunc(N * p - 1)
    rightBound = math.trunc(N * p + 1)

    # Pb(m | N,p) = 1 / delta * Phi( (m-mu) / delta )

    mu = N*p
    delta = math.sqrt(N*p*(1-p))

    for m in range(leftBound, rightBound + 1):
        result += 1 / delta * Phi( (m-mu) / delta )

    endTime = timeit.default_timer()
    workTime = (endTime - startTime) * 1000

    return result, f"Время работы программы:{workTime:.4f} (мс)"

def TheoryML_IntegralType(N,p):
    startTime = timeit.default_timer()
    result = 0

    leftBound = math.trunc(N * p - 1)
    rightBound = math.trunc(N * p + 1)

    # P (a <= (Ksi - mu) / delta < b) = F(b) - F(a)

    mu = N * p
    delta = math.sqrt(N * p * (1 - p))

    # Используя Z5:
    result = F((rightBound + 0.5 - mu) / delta)[0] - F((leftBound - mu) / delta)[0]

    endTime = timeit.default_timer()
    workTime = (endTime - startTime) * 1000

    return result, f"Время работы программы:{workTime:.4f} (мс)"

DataSet_of_p = [
        0.5, 0.5, 0.5, 0.5,
        0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.1, 0.05,
        0.05, 0.05, 0.05, 0.05,
        0.05, 0.01, 0.01, 0.01,
        0.01, 0.01, 0.01, 0.01
    ]

DataSet_of_N = [
        16, 24, 50, 100,
        10, 30, 50, 80,
        120, 250, 500, 20,
        60, 100, 240, 500,
        1000, 50, 100, 300,
        500, 800, 1200, 2500
    ]

DataSet_Bin_Result = []
DataSet_ThPuas_Result = []
DataSet_ThML_Int_Result = []
DataSet_ThML_Local_Result = []

DataSet_ThPuas_Inaccuracy = []
DataSet_ThML_Int_Inaccuracy = []
DataSet_ThML_Local_Inaccuracy = []

if __name__ == '__main__':
    "Сравнить точность аппроксимаций в предельных теоремах" \
    "при различных сочетаниях вероятности успеха p" \
    "и числа испытаний n" \
    "Вычислить вероятность попадания числа успехов в" \
    "интервал A = [Np-1, Np+1] по формулам биномиального распределения" \
    "и приближ. фформулам Пуассона и Муавра-Лапласа (как инт, так и лок)"

    for index in range(0, len(DataSet_of_p)):

        N = DataSet_of_N[index]
        p = DataSet_of_p[index]

        DataSet_Bin_Result.append(Bin(N, p))
        DataSet_ThPuas_Result.append(TheoryPuasson(N, p))
        DataSet_ThML_Int_Result.append((TheoryML_IntegralType(N, p)))
        DataSet_ThML_Local_Result.append(TheoryML_LocalType(N, p))

        DataSet_ThPuas_Inaccuracy.append( abs(Bin(N,p)[0] - TheoryPuasson(N, p)[0]) )
        DataSet_ThML_Int_Inaccuracy.append( abs(Bin(N,p)[0] - TheoryML_IntegralType(N, p)[0]) )
        DataSet_ThML_Local_Inaccuracy.append( abs(Bin(N,p)[0] - TheoryML_LocalType(N, p)[0]) )

    print("Биномиальная формула:")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_Bin_Result[_][0]:.4f}, {DataSet_Bin_Result[_][1]}", end = "\t")

    print(), print()

    print("Теорема Пуассона:")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_ThPuas_Result[_][0]:.4f}, {DataSet_ThPuas_Result[_][1]}", end = "\t")

    print()

    print("Погрешность Пуассона:")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_ThPuas_Inaccuracy[_]:.4f}", end = "\t")

    print(), print()

    print("Теорема Муавра-Лапласа (локальная):")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_ThML_Local_Result[_][0]:.4f}, {DataSet_ThML_Local_Result[_][1]}", end = "\t")

    print()

    print("Погрешность Муавра-Лапласа (локальная):")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_ThML_Local_Inaccuracy[_]:.4f}", end = "\t")

    print(), print()

    print("Теорема Муавра-Лапласа (интегральная):")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_ThML_Int_Result[_][0]:.4f}, {DataSet_ThML_Int_Result[_][1]}", end = "\t")

    print()

    print("Погрешность Муавра-Лапласа (интегральная):")
    for _ in range(0, len(DataSet_of_p)):
        if _%4 == 0 and _ != 0:
            print()
        print(f"{_+1}: {DataSet_ThML_Int_Inaccuracy[_]:.4f}", end = "\t")

    # ЗАДАНИЕ НОМЕР 88
    print(), print(), print()
    print("Задание номер 88.")
    print()

    p = 0.515
    n = 7280
    # Функция, которую мы интегрируем
    def func(x):
        return (1 / math.sqrt(2 * math.pi)) * math.exp(-x ** 2 / 2)
    # Нижний предел интегрирования
    a = float('-inf')

    # Значение интеграла
    integral_value = 0.98

    # Функция для вычисления значения интеграла
    def calculate_t(t):
        result, error = integrate.quad(func, a, t)
        return result - integral_value


    # Поиск значения t с использованием функции root_scalar
    solution = root_scalar(calculate_t, method='brentq', bracket=[-10, 10])

    t = solution.root

    print("Значение t:", t)

    print(
        f"[{(p + (math.pow(t,2) / n) - t * math.sqrt( (p*(1-p))/(n) + (math.pow(t,2))/(4*math.pow(n,2)) )) / (1 + (math.pow(t,2))/(n))};"
        f"{(p + (math.pow(t,2) / n) + t * math.sqrt( (p*(1-p))/(n) + (math.pow(t,2))/(4*math.pow(n,2)) )) / (1 + (math.pow(t,2))/(n))}]"
    )
