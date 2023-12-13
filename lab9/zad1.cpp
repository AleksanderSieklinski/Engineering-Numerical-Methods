#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

const int nx = 40;
const int ny = 40;
const int N = (nx + 1) * (ny + 1);
const double delta = 1.0;
const double delta_t = 1.0;
const double TA = 40.0;
const double TB = 0.0;
const double TC = 30.0;
const double TD = 0.0;
const double kB = 0.1;
const double kD = 0.6;
const int IT_MAX = 2000;

int main() {
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    gsl_matrix *B = gsl_matrix_alloc(N, N);
    gsl_vector *c = gsl_vector_alloc(N);
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            int l = i * (ny + 1) + j;
            if (i == 0 || i == nx) {  // warunek Dirichleta na lewej i prawej krawędzi
                gsl_matrix_set(A, l, l, 1.0);
                gsl_matrix_set(B, l, l, 1.0);
                gsl_vector_set(c, l, 0.0);
            } else if (j == 0) {  // warunek von Neumanna na dole
                gsl_matrix_set(A, l, l, 1.0 + 1.0 / (kD * delta));
                gsl_matrix_set(A, l, l + nx + 1, -1.0 / (kD * delta));
                gsl_vector_set(c, l, TD);
            } else if (j == ny) {  // warunek von Neumanna na górze
                gsl_matrix_set(A, l, l - nx - 1, -1.0 / (kB * delta));
                gsl_matrix_set(A, l, l, 1.0 + 1.0 / (kB * delta));
                gsl_vector_set(c, l, TB);
            } else {  // środek
                double a_ll = -2.0 * delta_t / (delta * delta) - 1.0;
                double a_others = delta_t / (2.0 * delta * delta);
                gsl_matrix_set(A, l, l - nx - 1, a_others);
                gsl_matrix_set(A, l, l - 1, a_others);
                gsl_matrix_set(A, l, l, a_ll);
                gsl_matrix_set(A, l, l + 1, a_others);
                gsl_matrix_set(A, l, l + nx + 1, a_others);
                double b_ll = 2.0 * delta_t / (delta * delta) - 1.0;
                double b_others = -delta_t / (2.0 * delta * delta);
                gsl_matrix_set(B, l, l - nx - 1, b_others);
                gsl_matrix_set(B, l, l - 1, b_others);
                gsl_matrix_set(B, l, l, b_ll);
                gsl_matrix_set(B, l, l + 1, b_others);
                gsl_matrix_set(B, l, l + nx + 1, b_others);
            }
        }
    }
    gsl_vector *T = gsl_vector_alloc(N);
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            int l = i * (ny + 1) + j;
            if (i == 0) {
                gsl_vector_set(T, l, TA);  // lewa krawędź
            } else if (i == nx) {
                gsl_vector_set(T, l, TC);  // prawa krawędź
            } else {
                gsl_vector_set(T, l, 0.0);  // środek
            }
        }
    }
    // Dekompozycja LU
    gsl_permutation *p = gsl_permutation_alloc(N);
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);
    // A*x = c
    gsl_vector *x = gsl_vector_alloc(N);
    gsl_linalg_LU_solve(A, p, c, x);
    // Implementacja algorytmu Cranka-Nicolsona
    gsl_vector *d = gsl_vector_alloc(N);
    for (int it = 0; it <= IT_MAX; ++it) {
        // d = B * T + c
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0.0, d);
        gsl_vector_add(d, c);
        // A*T = d
        gsl_linalg_LU_solve(A, p, d, T);
        if (it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000) {
            std::string it_str = it < 1000 ? "0" + std::to_string(it) : std::to_string(it);
            std::ofstream file("out/temperature_" + it_str + ".dat");
            for (int i = 0; i <= nx; ++i) {
                for (int j = 0; j <= ny; ++j) {
                    int l = i * (ny + 1) + j;
                    double T_ij = gsl_vector_get(T, l);
                    file << i << " " << j << " " << T_ij << "\n";
                }
                file << "\n";
            }
            file.close();
            std::ofstream file_laplacian("out/laplacian_" + it_str + ".dat");
            for (int i = 1; i < nx; ++i) {
                for (int j = 1; j < ny; ++j) {
                    int l = i * (ny + 1) + j;
                    double T_ij = gsl_vector_get(T, l);
                    double T_im1j = gsl_vector_get(T, l - ny - 1);
                    double T_ip1j = gsl_vector_get(T, l + ny + 1);
                    double T_ijm1 = gsl_vector_get(T, l - 1);
                    double T_ijp1 = gsl_vector_get(T, l + 1);
                    double laplacian = (T_im1j + T_ip1j + T_ijm1 + T_ijp1 - 4 * T_ij) / (delta * delta);
                    file_laplacian << i << " " << j << " " << laplacian << "\n";
                }
                file_laplacian << "\n";
            }
            file_laplacian.close();
        }
    }
    gsl_vector_free(d);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    return 0;
}