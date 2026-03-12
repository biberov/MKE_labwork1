#pragma once
#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Methods.h"
#include "TestFunctions.h"
#include <fstream>
using namespace std;

///Здесь реализованы функции для исследования сходимости:
///Вычисление погрешности для разных количеств разбиений
///Вычисление порядка точности(как меняется погрешность при измельчении сетки) 




// Структура для хранения результатов исследования сходимости
struct ConvergenceResult {
    std::vector<int> n_values;           // количество интервалов
    std::vector<double> errors;           // абсолютные погрешности
    std::vector<double> orders;           // наблюдаемый порядок сходимости на каждом шаге.
};


/// <summary>
/// Функция исследования сходимости составного метода трапеций
/// </summary>
/// <param name="f">Указатель на функцию(интегрируемая функция)</param>
/// <param name="exact_value">Точное значение интеграла функции (аналитическое)</param>
/// <param name="a">Нижний предел интегрирования</param>
/// <param name="b">Верхний предел интегрирования</param>
/// <param name="n_list">Список количества разбиений</param>
/// <returns>Возвращается объект структуры, который содержит три вектора с данными о процессе сходимости метода (n_values, errors, orders)</returns>
ConvergenceResult test_trapezoid_convergence(double (*f)(double), double exact_value,
    double a, double b, const std::vector<int>& n_list) {
    ConvergenceResult result;
    result.n_values = n_list;

    double prev_error = 0.0;
    for (int i = 0; i < n_list.size(); i++) {
        int n = n_list[i];
        double approx = composite_trapezoid(f, a, b, n);
        double error = std::abs(exact_value - approx);
        result.errors.push_back(error);

        // Вычисляем порядок точности (кроме первого шага)
        if (i == 0) {
            result.orders.push_back(0.0);
        }
        else {
            double ratio = prev_error / error;
            double order = log2(ratio); // так как шаг уменьшается в 2 раза
            result.orders.push_back(order);
        }
        prev_error = error;
    }

    return result;
/// Метод трапеций имеет второй порядок точности O(h²).
/// При удвоении количества разбиений (уменьшении шага в 2 раза) погрешность
/// должна уменьшаться примерно в 4 раза, что соответствует порядку точности 2.
}

/// <summary>
/// Функция исследования сходимости составного метода Симпсона
/// </summary>
/// <param name="f">Указатель на интегрируемую функцию</param>
/// <param name="exact_value">Точное значение интеграла (аналитическое)</param>
/// <param name="a">Нижний предел интегрирования</param>
/// <param name="b">Верхний предел интегрирования</param>
/// <param name="n_list">Список количества разбиений (должны быть четными)</param>
/// <returns>Возвращается объект структуры, который содержит три вектора с данными о процессе сходимости метода (n_values, errors, orders)</returns>
ConvergenceResult test_simpson_convergence(double (*f)(double), double exact_value,
    double a, double b, const std::vector<int>& n_list) {
    ConvergenceResult result;
    result.n_values = n_list;

    double prev_error = 0.0;
    for (int i = 0; i < n_list.size(); i++) {
        int n = n_list[i];
        double approx = composite_simpson(f, a, b, n);
        double error = std::abs(exact_value - approx);
        result.errors.push_back(error);

        if (i == 0) {
            result.orders.push_back(0.0);
        }
        else {
            double ratio = prev_error / error;
            double order = log2(ratio);
            result.orders.push_back(order);
        }
        prev_error = error;
    }

    return result;
/// Метод Симпсона имеет четвертый порядок точности O(h^4).
/// При удвоении количества разбиений (уменьшении шага в 2 раза) погрешность
/// должна уменьшаться примерно в 16 раз, что соответствует порядку точности 4.
}

/// <summary>
/// Вывод результатов исследования сходимости численного метода в файл
/// </summary>
/// <param name="method_name">Название метода (например, "Составная трапеция" или "Составной Симпсон")</param>
/// <param name="result">Структура с результатами исследования сходимости, содержащая векторы:
///        n_values - количество интервалов разбиения,
///        errors - значения погрешности,
///        orders - вычисленные порядки точности</param>
/// <param name="File">Ссылка на открытый файловый поток для записи результатов</param>
void print_convergence(const std::string& method_name, const ConvergenceResult& result, ofstream& File) {
   
    File << "\n" << method_name << " - исследование сходимости:\n";
    File << std::string(60, '-') << "\n";
    File << std::left << std::setw(30) << "n (интервалов)"
        << std::setw(30) << "Погрешность" << std::setw(30)
        << "Порядок точности\n";
    File << std::string(60, '-') << "\n";

    for (int i = 0; i < result.n_values.size(); i++) {
        File << std::left << std::setw(20) << result.n_values[i]
            << std::scientific << std::setprecision(8) << std::setw(30) << result.errors[i]
            << std::fixed << std::setprecision(2) << result.orders[i] << "\n";
    }
/// Функция форматирует вывод в виде таблицы с тремя колонками:
/// 1. Количество интервалов разбиения n
/// 2. Погрешность метода (в научном формате с точностью 8 знаков)
/// 3. Порядок точности (в фиксированном формате с точностью 2 знака)
/// 
/// Порядок точности вычисляется как log2(error_prev/error_curr) и показывает,
/// как быстро уменьшается погрешность при измельчении сетки.
/// Для метода трапеций ожидается порядок ~2, для метода Симпсона ~4.
}