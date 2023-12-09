#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

const double delta=0.01;
const double rho=1.0;
const double mi=1.0;
const int nx=200;
const int ny=90;
const int i_1=50;
const int j_1=55;
const int IT_MAX=20000; 

double x(int i){
    return delta * i;
}
double y(int j){
    return delta * j;
}
double q_wy(double q_we){
    return q_we * (pow(y(ny), 3) - pow(y(j_1), 3) - 3.0 * y(j_1) * pow(y(ny), 2) + 3.0 * pow(y(j_1), 2) * y(ny)) / pow(y(ny), 3);
}
bool is_edge(int i, int j){
    if ((i == 0 && j >= j_1 && j <= ny) || 
        (i == 0 && j >= 0 && j <= i_1) || 
        (i >= 0 && i <= i_1 && j == j_1) ||
        (i == nx && j >= 0 && j <= ny) || 
        (i >= i_1 && i <= nx && j == 0) || 
        (i == i_1 && j >= 0 && j <= j_1)){
        return true;
    }
    return false;
}
void warunek_brzegowy_psi(double psi[nx+1][ny+1], int q){
    // brzeg A
    for(int j=j_1; j<=ny; ++j){
		psi[0][j] = q/(2.0*mi) * (pow(y(j), 3)/3.0 - pow(y(j),2)*(y(j_1) + y(ny))/2.0 + y(j)*y(j_1)*y(ny));
	}
    // brzeg B
	for(int j=0; j<=ny; ++j){
		psi[nx][j] = q_wy(q)/(2.0*mi) * (pow(y(j), 3)/3.0 - pow(y(j),2)*y(ny)/2.0) + (q*pow(y(j_1),2)*(-y(j_1) + 3.0*y(ny))) / (12.0*mi);
	}
    // brzeg C
	for(int i=1; i<nx; ++i){
		psi[i][ny] = psi[0][ny];
	}
    // brzeg D
	for(int i=i_1; i<nx; ++i){
		psi[i][0] = psi[0][j_1];
	}
    // brzeg E
	for(int j=1; j<=j_1; ++j){
		psi[i_1][j] = psi[0][j_1];
	}
    // brzeg F
	for(int i=1; i<=i_1; ++i){
		psi[i][j_1] = psi[0][j_1];
	}
}
void warunek_brzegowy_zeta(double zeta[nx+1][ny+1], double psi[nx+1][ny+1], int q){
    // brzeg A
    for(int j=j_1; j<=ny; ++j){
        zeta[0][j] = q/(2.0*mi) * (2.0*y(j) - y(j_1) - y(ny));
    }
    // brzeg B
    for(int j=0; j<=ny; ++j){
        zeta[nx][j] = q_wy(q)/(2.0*mi) * (2.0*y(j) - y(ny));
    }
    // brzeg C
    for(int i=1; i<nx; ++i){
        zeta[i][ny] = 2.0/pow(delta,2) * (psi[i][ny-1] - psi[i][ny]);
    }
    // brzeg D
    for(int i=i_1+1; i<nx; ++i){
        zeta[i][0] = 2.0/pow(delta,2) * (psi[i][1] - psi[i][0]);
    }
    // brzeg E
    for(int j=1; j<j_1; ++j){
        zeta[i_1][j] = 2.0/pow(delta,2) * (psi[i_1+1][j] - psi[i_1][j]);
    }
    // brzeg F
    for(int i=1; i<=i_1; ++i){
        zeta[i][j_1] = 2.0/pow(delta,2) * (psi[i][j_1+1] - psi[i][j_1]);
    }
    // wierzcholek E/F
    zeta[i_1][j_1] = (zeta[i_1-1][j_1] + zeta[i_1][j_1-1]) / 2.0;
}
double control_error(double psi[nx+1][ny+1], double zeta[nx+1][ny+1]){
    double gamma = 0.0;
    int j_2 = j_1 + 1;
    for(int i=1; i<nx; ++i){
        gamma += psi[i+1][j_2] + psi[i-1][j_2] + psi[i][j_2+1] + psi[i][j_2-1] - 4.0 * psi[i][j_2] - pow(delta,2) * zeta[i][j_2];
    }
    return gamma;
}
void border_wipe(double psi[nx+1][ny+1], double zeta[nx+1][ny+1]){
    for(int i=0; i<=nx; ++i){
        psi[i][0] = NAN;
        psi[i][ny] = NAN;
        zeta[i][0] = NAN;
        zeta[i][ny] = NAN;
    }
    for(int j=0; j<=ny; ++j){
        psi[0][j] = NAN;
        psi[nx][j] = NAN;
        zeta[0][j] = NAN;
        zeta[nx][j] = NAN;
    }
}
void save_to_file(fstream &file, double tab[nx+1][ny+1]){
    for(int j=0; j<=ny; ++j){
        for(int i=0; i<=nx; ++i){
            file << tab[i][j] << " ";
        }
        file << endl;
    }
}
void zad(int q){
    double omega;
    double gamma;
    double psi[nx+1][ny+1];
    double zeta[nx+1][ny+1];
    double u[nx+1][ny+1];
    double v[nx+1][ny+1];
    printf("q = %d\n", q);
    for(int j=0; j<=ny; ++j){
        for(int i=0; i<=nx; ++i){
            psi[i][j]=0.0;
            zeta[i][j]=0.0;
            u[i][j]=0.0;
            v[i][j]=0.0;
        }
    }
    // ustalamy warunki brzegowe
    warunek_brzegowy_psi(psi, q);
    for (int it = 1; it<=IT_MAX; ++it) {
        if(it<2000){
            omega=0.0;
        }
        else{
            omega=1.0;
        }
        for (int j=1; j<ny; ++j){
            for (int i=1; i<nx; ++i){
                if (!is_edge(i, j)){
                    psi[i][j] = (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - pow(delta,2) * zeta[i][j]) / 4.0;
                    if(omega){
                        zeta[i][j] = (zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]) / 4.0 - omega * rho / (16.0 * mi) * ((psi[i][j+1] - psi[i][j-1]) * (zeta[i+1][j] - zeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j]) * (zeta[i][j+1] - zeta[i][j-1]));
                    }
                    u[i][j] = (psi[i][j+1] - psi[i][j-1]) / (2.0 * delta);
                    v[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2.0 * delta);
                }
            }
        }
        // modyfikujemy wartosci na brzegach
        warunek_brzegowy_zeta(zeta, psi, q);
        // kontrola bledu
        gamma = control_error(psi, zeta);
    }
    // usuwamy brzegi
    border_wipe(psi, zeta);
    // zapisujemy do plikow
    stringstream ss;
    ss << "out/u_" << q << ".txt";
    fstream file(ss.str(), ios::out);
    save_to_file(file, u);
    file.close();
    ss.str("");
    ss << "out/v_" << q << ".txt";
    file.open(ss.str(), ios::out);
    save_to_file(file, v);
    file.close();
    ss.str("");
    ss << "out/psi_" << q << ".txt";
    file.open(ss.str(), ios::out);
    save_to_file(file, psi);
    file.close();
    ss.str("");
    ss << "out/zeta_" << q << ".txt";
    file.open(ss.str(), ios::out);
    save_to_file(file, zeta);
    file.close();
}

int main() {
    int qvalues[] = {-4000, -1000, 4000};
    for(int i = 0; i < 3; ++i) {
        zad(qvalues[i]);
    }
    return 0;
}