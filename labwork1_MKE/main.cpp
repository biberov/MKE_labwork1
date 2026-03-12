#include <iostream>
#include <iomanip>
#include <vector>
#include "Methods.h"
#include "TestFunctions.h"
#include "ConvergenceTest.h"
#include <fstream>
using namespace std;

// Вспомогательная функция для форматированного вывода в файл
void print_header(const std::string& title, ofstream& File) {
    
    File << "\n" << std::string(80, '=') << "\n";
    File << "  " << title << "\n";
    File << std::string(80, '=') << "\n";
    
}
// Вспомогательная функция для форматированного вывода в файл
void print_subheader(const std::string& subtitle, ofstream& File) {
    
    File << "\n--- " << subtitle << " ---\n";
    
}

// Тестирование методов Гаусса на полиномах различной степени
void test_gauss_on_polynomials() {
/// Проверяет теоретическое свойство: метод Гаусса с n узлами точен для полиномов степени до 2n-1.
/// Для Гаусса-3 (n=3) точность ожидается до степени 5 включительно.
/// Для Гаусса-5 (n=5) точность ожидается до степени 9 включительно.
/// Для Гаусса-7 (n=7) точность ожидается до степени 13 включительно (все наши полиномы до x^6).
/// Результаты сохраняются в файл output.txt
    ofstream File("output.txt", ios::app);
    print_header("ТЕСТИРОВАНИЕ МЕТОДОВ ГАУССА", File);

    File << std::fixed << std::setprecision(12);
    File << "Интервал интегрирования: [" << A << ", " << B << "], h = " << H << "\n\n";

    // Точные значения для разных степеней
    double exact[] = { exact_poly0(), exact_poly1(), exact_poly2(),
                      exact_poly3(), exact_poly4(), exact_poly5(), exact_poly6() };

    // Массивы функций и их описаний
    double (*polys[])(double) = { poly0, poly1, poly2, poly3, poly4, poly5, poly6 };
    const char* poly_names[] = { "x^0", "x^1", "x^2", "x^3", "x^4", "x^5", "x^6" };


    File << std::left << std::setw(20) << "Полином"
         << std::setw(35) << "Точное значение"
        << std::setw(25) << "Гаусс-3"
         << std::setw(20) << "Гаусс-5"
        << std::setw(20) << "Гаусс-7" << "\n";
    File << std::string(90, '-') << "\n";

    for (int i = 0; i < 7; i++) {
        double g3 = gauss_3(polys[i], A, B);
        double g5 = gauss_5(polys[i], A, B);
        double g7 = gauss_7(polys[i], A, B);

        File << std::left << std::setw(10) << poly_names[i]
            << std::scientific << std::setprecision(8)
            << std::setw(20) << exact[i]
            << std::setw(20) << g3
            << std::setw(20) << g5
            << std::setw(20) << g7 << "\n";

        // Проверка точности для полиномов степени <= 2n-1
        if (i <= 5) { // до x^5 включительно - Гаусс-3 должен быть точен
            double eps = std::abs(exact[i] - g3);
            if (eps > 1e-12)
                File << "  !!! Гаусс-3 неточен для " << poly_names[i] << " (погр = " << eps << ")\n";
        }
        if (i <= 9) { // до x^9 - Гаусс-5 должен быть точен (но у нас только до x^6)
            double eps = std::abs(exact[i] - g5);
            if (eps > 1e-12)
                File << "  !!! Гаусс-5 неточен для " << poly_names[i] << " (погр = " << eps << ")\n";
        }
        // Гаусс-7 точен для всех наших полиномов (до x^6 < 13)
    }
    File.close();
}

// Тестирование составной формулы трапеций согласно пункту 3 лабораторной работы
void test_trapezoid() {
/// Выполняет три вида тестов:
/// 1. Проверка точности на полиномах 0 и 1 степени (где f'' = 0, метод должен быть точен)
/// 2. Проверка оценки погрешности на полиноме x^2 (где f'' = const)
/// 3. Исследование порядка сходимости на тригонометрической функции cos(x)
/// Результаты сохраняются в файл output.txt
    ofstream File("output.txt", ios::app);
    print_header("ТЕСТИРОВАНИЕ СОСТАВНОЙ ФОРМУЛЫ ТРАПЕЦИЙ", File);

    print_subheader("Полиномы, для которых метод точен (f'' = 0)", File);

    std::vector<int> n_values = { 1, 2, 4, 8, 16 };

    // Полином 0 степени (f(x)=1) - должен быть точен при любом n
    File << "\nf(x) = 1 (полином 0 степени):\n";
    for (int n : n_values) {
        double approx = composite_trapezoid(poly0, A, B, n);
        double error = std::abs(exact_poly0() - approx);
        File << "  n=" << n << ": approx=" << approx << ", error=" << error << "\n";
    }

    // Полином 1 степени (f(x)=x) - должен быть точен при любом n
    File << "\nf(x) = x (полином 1 степени):\n";
    for (int n : n_values) {
        double approx = composite_trapezoid(poly1, A, B, n);
        double error = std::abs(exact_poly1() - approx);
        File << "  n=" << n << ": approx=" << approx << ", error=" << error << "\n";
    }

    print_subheader("Полином, для которого оценка погрешности точна (f'' = const)", File);
    // f(x) = x^2, f''(x) = 2
    File << "\nf(x) = x^2 (полином 2 степени):\n";
    File << std::string(90, '-') << "\n";
    File << std::left << std::setw(10) << "n"
        << std::fixed << std::setprecision(8) << std::setw(30) << "Приближение"
        << std::setw(30) << "Точное"
         << std::setw(30) << "Погрешность"
        << std::setw(30) << "Оценка погрешности" << "\n";
    File << std::string(90, '-') << "\n";
    double h_prev = 0, error_prev = 0;
    for (int n : n_values) {
        double approx = composite_trapezoid(poly2, A, B, n);
        double error = exact_poly2() - approx; // знак учтен
        double error_abs = std::abs(error);

        // Оценка главного члена: -h^2*(b-a)*f''/12, f''=2
        double estimate = trapezoid_error_estimate(poly2, A, B, n, 2.0);

        File << std::left << std::setw(10) << n
            << std::fixed << std::setprecision(8) << std::setw(20) << approx
            << std::setw(20) << exact_poly2()
            << std::scientific << std::setprecision(8) << std::setw(20) << error
            << std::setw(30) << estimate << "\n";

        if (n > 1) {
            double order = log2(error_prev / error_abs);
            File << "    Порядок точности: " << std::fixed << std::setprecision(2) << order << "\n";
        }
        error_prev = error_abs;
    }

    print_subheader("Исследование порядка точности на тригонометрической функции", File);
    File << "\nf(x) = cos(x)\n";
    auto result = test_trapezoid_convergence(trig_func, exact_trig(), A, B, n_values);
    print_convergence("Составная трапеция", result, File);
    File.close();
}


// Тестирование составной формулы Симпсона согласно пункту 4 лабораторной работы
void test_simpson() {
/// Выполняет три вида тестов:
/// 1. Проверка точности на полиномах 0-3 степени (где f^(4) = 0, метод должен быть точен)
/// 2. Проверка оценки погрешности на полиноме x^4 (где f^(4) = const)
/// 3. Исследование порядка сходимости на тригонометрической функции sin(x)
/// Результаты сохраняются в файл output.txt
    ofstream File("output.txt", ios::app);
    print_header("ТЕСТИРОВАНИЕ СОСТАВНОЙ ФОРМУЛЫ СИМПСОНА", File);

    std::vector<int> n_values = { 2, 4, 8, 16, 32 }; // n должно быть четным

    print_subheader("Полиномы, для которых метод точен (до 3 степени)", File);

    // Полиномы 0-3 степени
    double (*polys[])(double) = { poly0, poly1, poly2, poly3 };
    const char* names[] = { "x^0", "x^1", "x^2", "x^3" };
    double exact_vals[] = { exact_poly0(), exact_poly1(), exact_poly2(), exact_poly3() };

    for (int idx = 0; idx < 4; idx++) {
        File << "\nf(x) = " << names[idx] << ":\n";
        for (int n : n_values) {
            double approx = composite_simpson(polys[idx], A, B, n);
            double error = std::abs(exact_vals[idx] - approx);
            File << "  n=" << n << ": error=" << std::scientific << error << "\n";
            if (error > 1e-12)
                File << "    !!! Ошибка не равна нулю!\n";
        }
    }

    print_subheader("Полином, для которого оценка погрешности точна (f^{(4)} = const)", File);
    // f(x) = x^4, f^{(4)}(x) = 24
    File << "\nf(x) = x^4 (полином 4 степени):\n";
    File << std::string(90, '-') << "\n";
    File << std::left << std::setw(10) << "n"
        << std::fixed << std::setprecision(8) << std::setw(30) << "Приближение"
        << std::setw(30) << "Точное"
        << std::scientific << std::setprecision(8) << std::setw(30) << "Погрешность"
        << std::setw(30) << "Оценка погрешности" << "\n";
    File << std::string(90, '-') << "\n";
    double error_prev = 0;
    for (int n : n_values) {
        double h = (B - A) / n;
        double approx = composite_simpson(poly4, A, B, n);
        double error = exact_poly4() - approx;
        double error_abs = std::abs(error);

        // Оценка: -h^4*(b-a)*f^{(4)}/180, f^{(4)}=24
        double estimate = -h * h * h * h * (B - A) * 24.0 / 180.0;

        File << std::left << std::setw(10) << n
            << std::fixed << std::setprecision(8) << std::setw(20) << approx
            << std::setw(20) << exact_poly4()
            << std::scientific << std::setprecision(8) << std::setw(20) << error
            << std::setw(30) << estimate << "\n";

        if (n > n_values[0]) {
            double order = log2(error_prev / error_abs);
            File << "    Порядок точности: " << std::fixed << std::setprecision(2) << order << "\n";
        }
        error_prev = error_abs;
    }

    print_subheader("Исследование порядка точности на тригонометрической функции", File);
    File << "\nf(x) = cos(x)\n";
    auto result = test_simpson_convergence(trig_func, exact_trig(), A, B, n_values);
    print_convergence("Составной Симпсон", result, File);
    File.close();
}

// Тестирование точности методов Гаусса согласно пункту 5 лабораторной работы.
void test_gauss_precision() {
/// Проверяет границы применимости методов Гаусса:
/// - Гаусс-3 должен быть точен для x^5 (степень 5 = 2*3-1), но неточен для x^6; 
/// - Гаусс-5 должен быть точен для x^5 и x^6 (т.к. 6 < 2*5-1 = 9). 
/// Результаты сохраняются в файл output.txt
    ofstream File("output.txt", ios::app);
    print_header("ТЕСТИРОВАНИЕ ТОЧНОСТИ МЕТОДОВ ГАУССА", File);

    File << std::fixed << std::setprecision(12);

    // Для Гаусса с 3 узлами
    print_subheader("Гаусс с 3 узлами (точен для полиномов до 5 степени)", File);
    File << "f(x) = x^5 (степень 5 - должна быть точность):\n";
    double g3_5 = gauss_3(poly5, A, B);
    double error_5 = std::abs(exact_poly5() - g3_5);
    File << "  Точное: " << exact_poly5() << ", приближение: " << g3_5
        << ", погрешность: " << error_5 << "\n";

    File << "\nf(x) = x^6 (степень 6 - погрешность должна быть ненулевой):\n";
    double g3_6 = gauss_3(poly6, A, B);
    double error_6 = std::abs(exact_poly6() - g3_6);
    File << "  Точное: " << exact_poly6() << ", приближение: " << g3_6
        << ", погрешность: " << error_6 << "\n";

    // Для Гаусса с 5 узлами
    print_subheader("Гаусс с 5 узлами (точен для полиномов до 9 степени)", File);
    File << "f(x) = x^5 (степень 5):\n";
    double g5_5 = gauss_5(poly5, A, B);
    error_5 = std::abs(exact_poly5() - g5_5);
    File << "  Погрешность: " << error_5 << "\n";

    File << "\nf(x) = x^6 (степень 6):\n";
    double g5_6 = gauss_5(poly6, A, B);
    error_6 = std::abs(exact_poly6() - g5_6);
    File << "  Погрешность: " << error_6 << "\n";
    File.close();
}


// Вычисление усилия N и момента M в оболочечном элементе (пункт 6).
void compute_shell_forces() {
/// Для заданного закона напряжения sigma(xi) = xi^3 - 2*xi^2 + 3 вычисляет:
/// - Усилие N = ∫ sigma(xi) dxi;
/// - Момент M = ∫ sigma(xi) * xi dxi.
/// Сравнивает точные аналитические значения с приближениями, полученными различными методами:
/// составная трапеция, составной Симпсон, Гаусс с 3, 5 и 7 узлами. 
/// Результаты сохраняются в файл output.txt
    ofstream File("output.txt", ios::app);
    print_header("ВЫЧИСЛЕНИЕ УСИЛИЙ И МОМЕНТОВ В ОБОЛОЧЕЧНОМ ЭЛЕМЕНТЕ", File);

    // Закон напряжения: sigma(xi) = xi^3 - 2*xi^2 + 3 (произвольный пример)
    auto sigma = [](double xi) -> double {
        return xi * xi * xi - 2.0 * xi * xi + 3.0;
        };

    // Для момента нужно sigma * xi
    auto sigma_times_xi = [](double xi) -> double {
        return xi * (xi * xi * xi - 2.0 * xi * xi + 3.0); // xi^4 - 2xi^3 + 3xi
        };

    File << "Закон напряжения: sigma(xi) = xi^3 - 2*xi^2 + 3\n";
    File << "Толщина оболочки h = " << H << "\n\n";

    // Точные значения (аналитическое интегрирование)
    double exact_N = (pow(B, 4) / 4.0 - 2.0 * pow(B, 3) / 3.0 + 3.0 * B) -
        (pow(A, 4) / 4.0 - 2.0 * pow(A, 3) / 3.0 + 3.0 * A);

    double exact_M = (pow(B, 5) / 5.0 - 2.0 * pow(B, 4) / 4.0 + 3.0 * pow(B, 2) / 2.0) -
        (pow(A, 5) / 5.0 - 2.0 * pow(A, 4) / 4.0 + 3.0 * pow(A, 2) / 2.0);

    File << "Точные значения:\n";
    File << "  N = " << exact_N << "\n";
    File << "  M = " << exact_M << "\n\n";
    File << std::string(90, '-') << "\n";
    File << std::left << std::setw(15) << "Метод"
        << std::setw(20) << "N"
        << std::setw(20) << "M"
        << std::setw(15) << "Погр N"
        << std::setw(15) << "Погр M" << "\n";
    File << std::string(85, '-') << "\n";

    // Составная трапеция (с большим числом разбиений)
    double N_trap = composite_trapezoid(sigma, A, B, 1000);
    double M_trap = composite_trapezoid(sigma_times_xi, A, B, 1000);
    File << std::left << std::setw(18) << "Трапеция"
        << std::fixed << std::setprecision(8) << std::setw(20) << N_trap
        << std::setw(20) << M_trap
        << std::scientific << std::setprecision(4)
        << std::setw(15) << std::abs(exact_N - N_trap)
        << std::setw(15) << std::abs(exact_M - M_trap) << "\n";

    // Составной Симпсон
    double N_simp = composite_simpson(sigma, A, B, 1000);
    double M_simp = composite_simpson(sigma_times_xi, A, B, 1000);
    File << std::left << std::setw(17) << "Симпсон"
        << std::fixed << std::setprecision(8) << std::setw(20) << N_simp
        << std::setw(20) << M_simp
        << std::scientific << std::setprecision(4)
        << std::setw(15) << std::abs(exact_N - N_simp)
        << std::setw(15) << std::abs(exact_M - M_simp) << "\n";

    // Гаусс с разным числом узлов
    double N_g3 = gauss_3(sigma, A, B);
    double M_g3 = gauss_3(sigma_times_xi, A, B);
    File << std::left << std::setw(15) << "Гаусс-3"
        << std::fixed << std::setprecision(8) << std::setw(20) << N_g3
        << std::setw(20) << M_g3
        << std::scientific << std::setprecision(4)
        << std::setw(15) << std::abs(exact_N - N_g3)
        << std::setw(15) << std::abs(exact_M - M_g3) << "\n";

    double N_g5 = gauss_5(sigma, A, B);
    double M_g5 = gauss_5(sigma_times_xi, A, B);
    File << std::left << std::setw(15) << "Гаусс-5"
        << std::fixed << std::setprecision(8) << std::setw(20) << N_g5
        << std::setw(20) << M_g5
        << std::scientific << std::setprecision(4)
        << std::setw(15) << std::abs(exact_N - N_g5)
        << std::setw(15) << std::abs(exact_M - M_g5) << "\n";

    double N_g7 = gauss_7(sigma, A, B);
    double M_g7 = gauss_7(sigma_times_xi, A, B);
    File << std::left << std::setw(15) << "Гаусс-7"
        << std::fixed << std::setprecision(8) << std::setw(20) << N_g7
        << std::setw(20) << M_g7
        << std::scientific << std::setprecision(4)
        << std::setw(15) << std::abs(exact_N - N_g7)
        << std::setw(15) << std::abs(exact_M - M_g7) << "\n";
    File.close();
}


int main() {


    ofstream File("output.txt");
    File << "ЛАБОРАТОРНАЯ РАБОТА №1\n";
    File << "Анализ погрешности вычисления усилий и моментов в оболочечном КЭ\n\n";

    // Пункт 2 (теоретический) - пояснения в начале программы
    File << "ТЕОРЕТИЧЕСКИЙ АНАЛИЗ:\n";
    File << "- Метод трапеций: 2-й порядок точности, точен для полиномов 0-1 степени\n";
    File << "- Метод Симпсона: 4-й порядок точности, точен для полиномов 0-3 степени\n";
    File << "- Гаусс с 3 узлами: точен для полиномов до 5 степени\n";
    File << "- Гаусс с 5 узлами: точен для полиномов до 9 степени\n";
    File << "- Гаусс с 7 узлами: точен для полиномов до 13 степени\n\n";
    File.close();

    // Пункт 3: Тестирование трапеции
    test_trapezoid();

    // Пункт 4: Тестирование Симпсона
    test_simpson();

    // Пункт 5: Тестирование Гаусса
    test_gauss_on_polynomials();
    test_gauss_precision();

    // Пункт 6: Расчет для оболочечного элемента
    compute_shell_forces();

    return 0;
}