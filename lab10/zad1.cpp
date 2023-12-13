#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

const int nx = 150;
const int nt = 1000;
const double delta = 0.1;
const double delta_t = 0.05;
const double tmax = nt * delta_t;
double xA = 7.5;
const double sigma = 0.5;
double b = 0.0; // damping coefficient
double alfa = 0.0; // spring constant
double aF = 0.0; // forcing term

std::vector<double> u(nx, 0.0);
std::vector<double> v(nx, 0.0);
std::vector<double> a(nx, 0.0);
std::vector<double> u0(nx, 0.0);
std::vector<double> vp(nx, 0.0);

double aFf(double x, double t, double tmax, double alfa, double xF) {
    return cos(50 * t / tmax) * (x==xF ? 1 : 0);
}
void initialize() {
    for(int i = 0; i < nx; i++) {
        double x = i * delta;
        u[i] = exp(-(x - xA)*(x - xA) / (2 * sigma * sigma));
        v[i] = 0.0;
    }
    u0 = u;
    for(int i = 1; i < nx - 1; i++) {
        a[i] = (u[i+1] - 2*u[i] + u[i-1]) / (delta * delta) - b * (u[i] - u0[i]) / delta_t + alfa*aFf(i*delta, 0, tmax, alfa, xA);
    }
}
void initialize_zero() {
    for(int i = 0; i < nx; i++) {
        u[i] = 0.0;
        v[i] = 0.0;
    }
    u0 = u;
    for(int i = 0; i < nx ; i++) {
        a[i] = (u[i+1] - 2*u[i] + u[i-1]) / (delta * delta) - b * (u[i] - u0[i]) / delta_t + alfa*aFf(i*delta, 0, tmax, alfa, xA);
    }
}
void update(double t) {
    for(int i = 0; i < nx ; i++) {
        vp[i] = v[i] + delta_t * a[i] / 2;
    }
    u0 = u;
    for(int i = 0; i < nx ; i++) {
        u[i] = u[i] + delta_t * vp[i];
    }
    for(int i = 0; i < nx ; i++) {
        a[i] = (u[i+1] - 2*u[i] + u[i-1]) / (delta * delta) - b * (u[i] - u0[i]) / delta_t + alfa*aFf(i*delta, t, tmax, alfa, xA);
        v[i] = vp[i] + delta_t * a[i] / 2;
    }
    // Warunki brzegowe
    u[0] = 0;
    u[nx-1] = 0;
    v[0] = 0;
    v[nx-1] = 0;
}
double calculate_energy() {
    double E = 0.0;
    double energy_1 = pow((u[1] - u[0]) / delta, 2);
    double energy_2 = pow((u[nx] - u[nx-1]) / delta, 2);
    E += delta / 4 * (energy_1 + energy_2);
    // Obliczanie energii dla środka
    for(int i = 0; i < nx ; i++) {
        energy_1 = v[i] * v[i];
        energy_2 = pow((u[i+1] - u[i-1]) / (2*delta), 2);
        E += (delta / 2)*(energy_1 + energy_2);
    }
    return E;
}
int main() {
    std::vector<double> b_values = {0.0, 0.1, 1.0, 1.0};  // wartości beta
    std::vector<double> a_values = {0.0, 0.0, 0.0, 1.0};  // wartości alfa
    for (int i = 0; i < b_values.size(); i++) {
        b = b_values[i];
        alfa = a_values[i];
        std::ostringstream stream;
        stream << std::fixed << std::setprecision(1) << b;
        std::string a_str = stream.str();
        stream.str("");
        stream << std::fixed << std::setprecision(1) << alfa;
        std::string b_str = stream.str();
        stream.str("");
        printf("b = %f, alfa = %f\n", b, alfa);
        std::ofstream file("out/output_alfa_" + b_str + "_beta_" + a_str + ".txt");
        std::ofstream energy_file("out/energy_alfa_" + b_str + "_beta_" + a_str + ".txt");
        if (alfa == 1.0) {
            xA = 2.5;
            initialize_zero();
        }
        else {
            initialize();
        }
        for(double t = 0; t < tmax; t+=delta_t) {
            update(t);
            double E = calculate_energy();
            energy_file << E << "\n";
            for(int i = nx; i > 0; i--) {
                file << u[i] << " ";
            }
            file << "\n";
        }
        file.close();
        energy_file.close();
    }
    return 0;
}