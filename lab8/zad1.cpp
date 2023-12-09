#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;
const int nx = 400;
const int ny = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const float delta = 0.01;
const float sigma = 10.0 * delta;
const float x_a = 0.45;
const float y_a = 0.45;
const int it_max = 10000;

float x(int i){
    return delta*i;
}
float y(int j){
    return delta*j;
}
float t(int k, float dt){
    return dt*k;
}
float u(float x, float y){
    return exp(-(pow((x - x_a),2) + pow((y - y_a),2)) / (2.0 * pow(sigma,2))) / (2.0 * M_PI * pow(sigma,2));
}
bool is_in(int i, int j){
    return (i >= i_1 && i <= i_2 && j >= 0 && j <= j_1);
}
void save_to_file(fstream& file, vector<vector<float>> u0){
    file << fixed;
    for(int j=0; j<=ny; ++j){
        for(int i=0; i<=nx; ++i){
            file << i << " " << j << " "<< u0[i][j] << "\n";
        }
        file << "\n";
    }
}
void calc_loop(float dd, float vx[nx+1][ny+1], float vy[nx+1][ny+1], float dt){
    vector<vector<float>> u0(nx+1, vector<float>(ny+1));
    vector<vector<float>> u1(nx+1, vector<float>(ny+1));
    // wypełnienie zerami i warunki brzegowe
    for(int j=0; j<=ny; ++j){
        for(int i=0; i<=nx; ++i){
            u0[i][j] = u(x(i), y(j));
            u1[i][j] = 0.0;
        }
    }
    // Stworzenie wektora do przechowywania danych
    vector<vector<vector<float>>> stored_data(50);
    // Zmienne pomocnicze w celu przyspieszenia obliczeń
    float power = pow(delta,2);
    float dtdd2 = dt * dd / 2.0;
    float delta2 = delta * 2.0;
    double divider = 1.0 + (2.0 * dd * dt) / (power);
    auto start = chrono::high_resolution_clock::now();
    int k=0;
    stringstream ss;
    ss << "out/zad5_c" << dd << ".txt";
    fstream file_c(ss.str(), ios::out);
    ss.str("");
    ss << "out/zad5_xsr" << dd << ".txt";
    fstream file_xsr(ss.str(), ios::out);
    file_c << fixed;
    file_xsr << fixed;
    for(int it=0; it<it_max; ++it){
        u1 = u0;
        for(int kit=1; kit<=20; ++kit){
            for(int i=0; i<=nx; ++i){
                for(int j=1; j<ny; ++j){
                    if (is_in(i,j)){
                        continue;
                    }
                    else if(i==0){
                        u1[i][j] = (u0[i][j] - dt * vx[i][j] / 2.0 * (((u0[i+1][j] - u0[nx][j]) / (delta2)) + ((u1[i+1][j] - u1[nx][j]) / (delta2))) - dt * vy[i][j] / 2.0 * (((u0[i][j+1] - u0[i][j-1]) / (delta2)) + ((u1[i][j+1] - u1[i][j-1]) / (delta2))) + dtdd2
                                   * ((u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4.0 * u0[i][j]) / (power) + ((u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1]) / (power)))) / divider;
                    }
                    else if(i==nx){
                        u1[i][j] = (u0[i][j] - dt * vx[i][j] / 2.0 * (((u0[0][j] - u0[i-1][j]) / (delta2)) + ((u1[0][j] - u1[i-1][j]) / (delta2))) - dt * vy[i][j] / 2.0 * (((u0[i][j+1] - u0[i][j-1]) / (delta2)) + ((u1[i][j+1] - u1[i][j-1]) / (delta2))) + dtdd2
                                   * ((u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0 * u0[i][j]) / (power) + ((u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1]) / (power)))) / divider;
                    }
                    else{
                        u1[i][j] = (u0[i][j] - dt * vx[i][j] / 2.0 * (((u0[i+1][j] - u0[i-1][j]) / (delta2)) + ((u1[i+1][j] - u1[i-1][j]) / (delta2))) - dt * vy[i][j] / 2.0 * (((u0[i][j+1] - u0[i][j-1]) / (delta2)) + ((u1[i][j+1] - u1[i][j-1]) / (delta2))) + dtdd2
                                   * ((u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0 * u0[i][j]) / (power) + ((u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1]) / (power)))) / divider;
                    }
                }
            }
        }
        u0 = u1;
        if(it%10==0){
            auto end = chrono::high_resolution_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            cout << it << endl;
            cout << "Time elapsed: " << elapsed << " ms" << endl;
            start = chrono::high_resolution_clock::now();
        }
        // zapis do pliku
        float xsr_sum = 0.0;
        float c_sum = 0.0;
        for(int i=0; i<=nx; ++i){
            for(int j=0; j<=ny; ++j){
                c_sum += u0[i][j] * power;
                xsr_sum += u0[i][j] * x(i) * power;
            }
        }
        file_xsr << xsr_sum << endl;
        file_c << c_sum << endl;
        if(it == (it_max*k)/50){
            stored_data[k] = u0;
            k++;
        }
    }
    // create 50 files with u, save data to them and close them
    for(int i=0; i<50; ++i){
        stringstream ss;
        ss << "out/zad5_it=" << i << "_" << dd << ".txt";
        fstream file(ss.str(), ios::out);
        file << fixed;
        save_to_file(file, stored_data[i]);
        file.close();
    }
}
void zad(){
    // wczytanie psi
    ifstream file("psi.dat");
    string line;
    float psi[nx+1][ny+1];
    while (getline(file, line)){
        istringstream iss(line);
        size_t i, j;
        float temp_psi;
        if (iss >> i >> j >> temp_psi) {
            psi[i][j] = temp_psi;
        } else {
            cerr << "Error parsing line: " << line << endl;
        }
    }
    // obliczenie prędkości
    float vx[nx+1][ny+1];
    float vy[nx+1][ny+1];
    // wypełnienie zerami
    for (int i=0; i<=nx; ++i) {
        for (int j=0; j<=ny; ++j) {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    // obliczenie prędkości
    for (int i=1; i<nx; ++i) {
        for (int j=1; j<ny; ++j) {
            vx[i][j] = (psi[i][j+1] - psi[i][j-1]) / (delta * 2.0);
            vy[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (delta * 2.0);
        }
    }
    // warunki brzegowe zastawka
    for(int i=i_1; i<=i_2; ++i){
        for(int j=0; j<=j_1; ++j){
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    // warunki brzegowe górna i dolna ściana
    for(int i=1; i<nx; ++i){
        vx[i][0] = 0.0;
        vy[i][ny] = 0.0;
    }
    // warunki brzegowe lewa i prawa ściana
    for(int j=0; j<=ny; ++j){
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx-1][j];
    }
    for (int k=0;k<2;k++){
        stringstream ss;
        if (k==0){
        ss << "out/zad3_vx.txt";
        }
        else{
        ss << "out/zad3_vy.txt";
        }
        fstream file(ss.str(), ios::out);
        file << fixed;
        for(int j=0; j<=ny; ++j){
            for(int i=0; i<=nx; ++i){
                if (k==0){
                file << vx[i][j] << " ";
                }
                else{
                file << vy[i][j] << " ";
                }
            }
            file << endl;
        }
        file.close();
    }
    // obliczenie dt na podstawie prędkości
    float v_max = 0.0;
    float v_temp = 0.0;
    for (int i=0; i<=nx; ++i) {
        for (int j=0; j<=ny; ++j) {
            v_temp = sqrt(pow(vx[i][j],2) + pow(vy[i][j],2));
            if (v_temp > v_max) {
                v_max = v_temp;
            }
        }
    }
    float dt = delta / (4.0 * v_max);
    stringstream ss;
    ss << "out/zad5_tt.txt";
    fstream file_tt(ss.str(), ios::out);
    file_tt << fixed;
    for(int it=0;it<it_max;it++){
        file_tt << t(it, dt) << endl;
    }
    file_tt.close();
    ss.str("");

    calc_loop(0.0, vx, vy, dt);
    calc_loop(0.1, vx, vy, dt);
}
int main(){
    zad();
    return 0;
}