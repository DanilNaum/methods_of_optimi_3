#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
using namespace std;

int coun = 0;

double Func(complex<double> x, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    coun++;
    complex<double>A = complex<double>(6.0, 12.0);
    complex<double>B = complex<double>(12.0, 6.0);
    complex<double>C = complex<double>(6.0, 7.0);
    complex<double>D = first_root;
    complex<double>E = second_root;
    double f;
    if (number_roots == 0)
        f = abs((x - A) * (x - B) * (x - C));
    else if (number_roots == 1)
        f = abs(x * x + (D - A - B - C) * x + (D * D - A * D - B * D - C * D + A * B + A * C + B * C));
    else
        f = abs(x + E + D - A - B - C);
    return (f);
}

double DFunc(complex<double> x, bool per1, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    coun++;

    complex<double>A = complex<double>(6.0, 12.0);
    complex<double>B = complex<double>(12.0, 6.0);
    complex<double>C = complex<double>(6.0, 7.0);
    complex<double>D = first_root;
    double n_r1cof_b_real = (D - A - B - C).real();
    double n_r1cof_b_imag = (D - A - B - C).imag();
    double n_r1cof_c_real = (D * D - A * D - B * D - C * D + A * B + A * C + B * C).real();
    double n_r1cof_c_imag = (D * D - A * D - B * D - C * D + A * B + A * C + B * C).imag();
    complex<double>E = second_root;
    double n_r2cof_b_real = (E + D - A - B - C).real();
    double n_r2cof_b_imag = (E + D - A - B - C).imag();
    double f;
    double x1 = x.real();
    double y1 = x.imag();
    if (!per1) {
        if (number_roots == 0)
            f = 2 * (-120 * pow(x1, 4) + 3 * pow(x1, 5) + 2 * pow(x1, 3) * (1165 - 50 * y1 + 3 * pow(y1, 2)) - 18 * pow(x1, 2) * (1443 - 138 * y1 + 8 * pow(y1, 2)) + x1 * (165240 - 27252 * y1 + 2402 * pow(y1, 2) - 100 * pow(y1, 3) + 3 * pow(y1, 4)) - 6 * (78300 - 19140 * y1 + 2283 * pow(y1, 2) - 138 * pow(y1, 3) + 4 * pow(y1, 4)));
        else if (number_roots == 1)
            f = 2 * (pow(n_r1cof_b_real, 2) * x1 + pow(n_r1cof_b_imag, 2) * x1 + 3 * n_r1cof_b_real * pow(x1, 2) + 2 * pow(x1, 3) + n_r1cof_c_real * (n_r1cof_b_real + 2 * x1) + 2 * n_r1cof_b_imag * x1 * y1 + (n_r1cof_b_real + 2 * x1) * pow(y1, 2) + n_r1cof_c_imag * (n_r1cof_b_imag + 2 * y1));
        else
            f = 2 * (n_r2cof_b_real + x1);
        return (f);
    }
    else {
        if (number_roots == 0)
            f = 2 * (-502200 + pow(x1, 3) * (828 - 96 * y1) + 178200 * y1 - 27918 * pow(y1, 2) + 2474 * pow(y1, 3) - 125 * pow(y1, 4) + 3 * pow(y1, 5) + pow(x1, 4) * (-25 + 3 * y1) + 2 * pow(x1, 2) * (-6813 + 1201 * y1 - 75 * pow(y1, 2) + 3 * pow(y1, 3)) - 12 * x1 * (-9570 + 2283 * y1 - 207 * pow(y1, 2) + 8 * pow(y1, 3)));
        else if (number_roots == 1)
            f = 2 * (n_r1cof_b_imag * pow(x1, 2) + n_r1cof_c_imag * (n_r1cof_b_real + 2 * x1) + pow(n_r1cof_b_real, 2) * y1 + pow(n_r1cof_b_imag, 2) * y1 + 2 * n_r1cof_b_real * x1 * y1 + 2 * pow(x1, 2) * y1 + 3 * n_r1cof_b_imag * pow(y1, 2) + 2 * pow(y1, 3) - n_r1cof_c_real * (n_r1cof_b_imag + 2 * y1));
        else
            f = 2 * (n_r2cof_b_imag + y1);

        return (f);
    }
}

double g(complex<double> z) {
    const double R = sqrt(84);
    double x = z.real(), y = z.imag();
    return x * x + y * y - R * R;
}

//Градиент функции ограничения
complex<double> g_grad(complex<double> z) {
    return z * 2.;
}
complex<double>Grad(complex<double> x, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    complex<double>gr = complex<double>(DFunc(x, bool(0), number_roots, first_root, second_root), DFunc(x, bool(1), number_roots, first_root, second_root));
    return gr;
}

void surf_border(complex<double> z) {
    double eps = 0.;
    const double EPS = 1e-5;
    const double L = 0.5;
    const double APPROACH_STEP = 1.;
    int k = 1;

    complex<double> fg = Grad(z);
    //быстро и грубо шагами длиной APPROACH_STEP подходим к границе
    while (abs(g(z)) > EPS) {
        fg /= abs(fg);
        //если шаг выбивает нас за границу допустимой области,
        //линейно аппроксимируем g(x, y) и делаем шаг с такой длиной,
        //чтобы попасть ровненько на границу (g(x,y) = 0)
        //подробнее - в методичке
        if (g(z - APPROACH_STEP * fg) > EPS) {
            z -= g(z) / (g_grad(z).real()* fg.real()+ g_grad(z).imag() * fg.imag()) * fg;
        }
        else {
            z -= APPROACH_STEP * fg;
        }
        fg = Grad(z);
        printf("%d; %d; (%f, %f);  %f; (%f, %f); -; -; -; -\n",
            k, coun, z.real(), z.imag(),  g(z),
            fg.real(), fg.imag());
        ++k;
    }
    
    complex<double> gg = g_grad(z), prev_tan = 0.;
    double cos = (-fg.real()*gg.real() + (-fg).imag() * gg.imag())/(abs(-fg)*abs(gg));
    const double TARGET_COS = 0.9994;
    double alpha = 0.1;
    //начинаем скользить по границе допустимой области
    //останавливаемся, когда косинус между антиградиентом f и градиентом g
    //станет больше либо равен TARGET_COS (т.е. меньше двух градусов)
    do {
        //Проекция антиградиента на касательную плоскость
        complex<double> tan = -fg - gg * (-fg.real()* gg.real() + -fg.imag() * gg.imag()) / (abs(gg)* abs(gg));
        //если проекция антиградиента резко изменила направление
        //т.е. угол между прыдудущей проекцией и текущей тупой,
        //то это означает, что мы сделали слишком большой шаг и проскочили минимум
        //раз шаг слишком большой, то делим его пополам
        if ((tan.real()* prev_tan.real()+ tan.imag() * prev_tan.imag()) < 0)
            alpha /= 2;
        prev_tan = tan;

       tan /= abs(tan);
        z += tan * alpha;
 
        //Если слишком сильно вышли за границу, то возвращаемся обратно.
        //возвращаемся по направлению антиградиента с такой длиной, 
        //что линейная аппроксимация g в новой точке равна 0
        if (g(z) > EPS)
            z -= g(z) / (abs(gg)* abs(gg)* abs(gg) * abs(gg)) * gg;

        fg = Grad(z);
        gg = g_grad(z);
        cos = (-fg.real()*gg.real() + (-fg).imag() * gg.imag())/(abs(-fg)*abs(gg));
        ++k;

        printf("%d; %d; (%f, %f);  %f; (%f, %f); (%f, %f); (%f, %f); %f; %f\n",
            k, coun, z.real(), z.imag(), g(z),
            fg.real(), fg.imag(), gg.real(), gg.imag(), tan.real(), tan.imag(), cos,
            acos(cos) * 180 / 3.1415);
    } while (cos < TARGET_COS);
}


int main() {
    surf_border(complex<double>(0., 0.));
    return 0;
}