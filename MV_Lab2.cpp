#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> 

using namespace std;
const double EXACT_VALUE = 10.5844484649508098; //точное значение 

// Функция для вычисления частичной суммы ряда
double partialSum(int n, double alpha) {
    double sum = 0.0;
    for (int i = 1; i <= n; ++i) {
        sum += 1.0 / pow(i, alpha);
    }
    return sum;
}

double richardsonExtrapolation(double zn, double z2n, double Q, double k) {
    double Qk = pow(Q, k);
    return (Qk * z2n - zn) / (Qk - 1.0);
}

int main() {
    setlocale(LC_ALL, "ru");

    double alpha = 1.1;               // параметр α ряда
    double Q = 2.0;                   // коэффициент увеличения n
    double k1 = alpha - 1.0;          // первый показатель k1 = α-1
    int startN = 1;                   // начальное n (2^i)
    int numLevels = 17;

    vector<int> nValues(numLevels); //n
    vector<double> z(numLevels); // Z_n вычисленная сумма для n

    for (int i = 0; i < numLevels; ++i) {
        nValues[i] = startN * pow(Q, i);
        z[i] = partialSum(nValues[i], alpha);
    }

    vector<double> z1(numLevels, 0.0);    // Z_n^(1) первое экстраполированное значение
    vector<double> z2(numLevels, 0.0);    // Z_n^(2) второе экстраполированное значение

    for (int i = 1; i < numLevels; ++i) {
        // Z_n = z[i-1], Z_2n = z[i]
        z1[i] = richardsonExtrapolation(z[i - 1], z[i], Q, k1);
    }

    double k2 = k1 + 1.0;
    for (int i = 2; i < numLevels; ++i) {
        // Z_n^(1) = z1[i-1], Z_2n^(1) = z1[i]
        z2[i] = richardsonExtrapolation(z1[i - 1], z1[i], Q, k2);
    }


    vector<double> delta_n(numLevels, 0.0); //оценка погрешности после первой экстраполяции
    vector<double> delta_n1(numLevels, 0.0);// оценка погрешности после 2 экстраполяции

    cout << "Результаты экстраполяции для ряда ∑1/i^" << alpha << " (Q=" << Q << ", k=" << k1 << ")\n";
    cout << "Точное значение (S_inf): Z_TOЧH = " << scientific << EXACT_VALUE << "\n\n";

    cout << left << setw(10) << "n"
        << setw(12) << "Z_n-Z_TOЧH"
        << setw(14) << "Z_n-Z_n(1)"
        << setw(15) << "Z_n(1)-Z_TOЧH"
        << setw(12) << "delta"
        << setw(15) << "Z_n(2)-Z_TOЧH"
        << setw(12) << "delta(1)" << endl;

    cout << string(88, '-') << endl;

    cout << scientific << setprecision(3);

    for (int i = 0; i < numLevels; ++i) {
        
        double error_z = z[i] - EXACT_VALUE; // Погрешность (Z_n - Z_TOЧH) погрешность вычесленной суммы
        double error_z1 = z1[i] - EXACT_VALUE; // Z_n(1) - Z_TOЧH погрешность после первой экстраполяции
        double error_z2 = (i >= 2) ? z2[i] - EXACT_VALUE : 0.0; // Z_n(2)-Z_TOЧH – погрешность после 2 экстраполяции

        double diff_z_z1 = z[i] - z1[i]; // Z_n-Z_n(1) Разность между исходной и экстраполированной оценкой.


        delta_n[i] = error_z1 / error_z;
        delta_n1[i] = error_z2 / error_z1;



        cout << left << setw(8) << nValues[i]
            << setw(12) << error_z
            << setw(14) << diff_z_z1
            << setw(15) << error_z1
            << setw(12) << delta_n[i]
            << setw(15) << error_z2
            << setw(12) << delta_n1[i] << endl;
    }

    return 0;
}