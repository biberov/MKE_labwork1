#pragma once
#pragma once
#include <cmath>
#include <functional>
/// Здесь определяются:
/// Полиномы степеней 0 1 2 3 4 5 6
/// Функции для вычисления точных интегралов(первообразные).
/// Аналитические оценки главного члена погрешности.


// Константы для интервала интегрирования (толщина оболочки)
const double H = 2.0;      // толщина
const double A = -H / 2.0; // -1
const double B = H / 2.0;  // +1

//-----------------------------------------------------------------------------
// ПОЛИНОМЫ РАЗНЫХ СТЕПЕНЕЙ (для тестирования)
//-----------------------------------------------------------------------------

// Полином 0 степени: f(x) = 1
double poly0(double x) { return 1.0; }

// Полином 1 степени: f(x) = x
double poly1(double x) { return x; }

// Полином 2 степени: f(x) = x^2
double poly2(double x) { return x * x; }

// Полином 3 степени: f(x) = x^3
double poly3(double x) { return x * x * x; }

// Полином 4 степени: f(x) = x^4
double poly4(double x) { double x2 = x * x; return x2 * x2; }

// Полином 5 степени: f(x) = x^5
double poly5(double x) { double x2 = x * x; double x3 = x2 * x; return x3 * x2; }

// Полином 6 степени: f(x) = x^6
double poly6(double x) { double x2 = x * x; double x4 = x2 * x2; return x4 * x2; }

// Тригонометрическая функция (для демонстрации порядка точности)
double trig_func(double x) { return cos(x); }

//-----------------------------------------------------------------------------
// ТОЧНЫЕ ЗНАЧЕНИЯ ИНТЕГРАЛОВ (первообразные)
//-----------------------------------------------------------------------------

// ∫ x^n dx от A до B
double exact_integral_power(int n) {
    if (n == -1) return log(B) - log(A); // для 1/x, но нам не нужно
    return (pow(B, n + 1) - pow(A, n + 1)) / (n + 1);
}

// Удобные функции для точных значений интегралов от наших полиномов
double exact_poly0() { return exact_integral_power(0); } // = H
double exact_poly1() { return exact_integral_power(1); } // = 0
double exact_poly2() { return exact_integral_power(2); } // = H^3/12
double exact_poly3() { return exact_integral_power(3); } // = 0
double exact_poly4() { return exact_integral_power(4); } // = H^5/80
double exact_poly5() { return exact_integral_power(5); } // = 0
double exact_poly6() { return exact_integral_power(6); } // = H^7/448

// Точное значение для cos(x)
double exact_trig() { return sin(B) - sin(A); }

//-----------------------------------------------------------------------------
// ГЛАВНЫЙ ЧЛЕН ПОГРЕШНОСТИ (аналитические оценки)
//-----------------------------------------------------------------------------

// Для метода трапеций на одном интервале: -(b-a)^3/12 * f''(ξ)
// Для составной формулы: -h^2/12 * (b-a) * f''(ξ)
double trapezoid_error_estimate(double (*f)(double), double a, double b, int n, double f_double_prime) {
    double h = (b - a) / n;
    return -h * h * (b - a) * f_double_prime / 12.0;
}

// Для метода Симпсона: -h^4/180 * (b-a) * f^{(4)}(ξ)
double simpson_error_estimate(double (*f)(double), double a, double b, int n, double f_fourth_prime) {
    double h = (b - a) / n;
    return -h * h * h * h * (b - a) * f_fourth_prime / 180.0;
}

// Для метода Гаусса с n узлами: (b-a)^{2n+1} * (n!)^4 / ((2n+1) * ((2n)!)^3) * f^{(2n)}(ξ)
// Но мы для тестов будем просто сравнивать с полиномами степени 2n