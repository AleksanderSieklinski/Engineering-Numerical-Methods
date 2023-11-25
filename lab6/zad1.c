#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mgmres.h"
const double delta = 0.1;
int nx = 4;
int ny = 4;
double V1 = 10;
double V2 = -10;
double V3 = 10;
double V4 = -10;
double ep1 = 1;
double ep2 = 1;
double rho(int i, int j,double xmax, double ymax){
	double sigma = 0.1 * xmax;
	double rho1 = exp(-pow((delta * i - 0.25 * xmax) / sigma, 2) - pow((delta * j - 0.5 * ymax) / sigma, 2));
	double rho2 = -exp(-pow((delta * i - 0.75 * xmax) / sigma, 2) - pow((delta * j - 0.5 * ymax) / sigma, 2));
	return rho1 + rho2;
}
int j(int l, int nx){return l / (nx + 1);}
int i(int l, int nx){return l - j(l, nx) * (nx + 1);}
void metodaPoissona(char useRho, FILE *f, char writeFile){
	int N =(nx + 1) * (ny + 1);
	double xmax = delta * nx;
	double ymax = delta * ny;
	double *a = calloc(5 * N, sizeof(double));
	double *b = calloc(N, sizeof(double));
	double *V = calloc(N, sizeof(double));
	int *ja = calloc(5 * N, sizeof(int));
	int *ia = calloc(N + 1, sizeof(int));
	printf("V1 = %f, V2 = %f, V3 = %f, V4 = %f\n", V1, V2, V3, V4);
	for(int i = 0; i < N+1; i++){
		ia[i] = -1;
	}
	int k = -1;
	for(int l = 0; l < N; l++){
		int brzeg = 0;
		double vb = 0;
		if(i(l, nx) == 0){
			brzeg = 1;
			vb = V1;
		}
		if(j(l, nx) == ny){
			brzeg = 1;
			vb = V2;
		}
		if(i(l, nx) == nx){
			brzeg = 1;
			vb = V3;
		}
		if(j(l, nx) == 0){
			brzeg = 1;
			vb = V4;
		}
		if(brzeg == 1){
			b[l] = vb;
		}
		else{
			if(useRho == 'y'){
				b[l] = -rho(i(l, nx), j(l, nx), xmax, ymax);
			}
			else
				b[l] = 0;
		}
		ia[l] = -1;
		//lewa skrajna przekatna
		if(l - nx - 1 >= 0 && brzeg == 0){
			k++;
			if(ia[l] < 0){
				ia[l] = k;
			}
			double epsilon = i(l, nx) <= nx / 2 ? ep1 : ep2;
			a[k] = epsilon / pow(delta, 2);
			ja[k] = l - nx - 1;
		}
		//poddiagonala
		if(l - 1 >= 0 && brzeg == 0){
			k++;
			if(ia[l] < 0){
				ia[l] = k;
			}
			double epsilon = i(l, nx) <= nx / 2 ? ep1 : ep2;
			a[k] = epsilon / pow(delta, 2);
			ja[k] = l - 1;
		}
		//diagonala
		k++;
		if(ia[l] < 0){
			ia[l] = k;
		}
		if(brzeg == 0){
			double epsilon0 = i(l, nx) <= nx / 2 ? ep1 : ep2;
			double epsilon1 = i(l + 1, nx) <= nx / 2 ? ep1 : ep2;
			double epsilon2 = i(l + nx + 1, nx) <= ny / 2 ? ep1 : ep2;
			a[k] = -(2 * epsilon0 + epsilon1 + epsilon2) / pow(delta, 2);

		}
		else
			a[k] = 1;
		ja[k] = l;
		//nad diagonala
		if(l < N && brzeg == 0){
			k++;
			double epsilon = i(l + 1, nx) <= nx / 2 ? ep1 : ep2;
			a[k] = epsilon / pow(delta, 2);
			ja[k] = l + 1;
		}
		//prawa skrajna przekatna
		if(l < N - nx - 1 && brzeg == 0){
			k++;
			double epsilon = i(l + nx + 1, nx) <= nx / 2 ? ep1 : ep2;
			a[k] = epsilon / pow(delta, 2);
			ja[k] = l + nx + 1;
		}
	}
	int nz_num = k + 1;
	ia[N] = nz_num;
	pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, 500, 500, 1e-8, 1e-8);

	if(writeFile == 'm'){
		fprintf(f, "# l \t i_l \t j_l \t b[l]\n");
		for(int l = 0; l < N; l++){
			fprintf(f, "%d \t\t%d \t\t%d \t\t%lf\n", l, i(l, nx), j(l, nx),b[l]);
		}
		fprintf(f, "# k \t a[l]\n");
		for(int l = 0; l < 5*N; l++){
			fprintf(f, "%d \t %lf\n", l,a[l]);
		}
	}
	else if(writeFile == 'w'){
		for(int i = 0; i < N-ny-1; i++){
			fprintf(f, "%f ",V[i]);
			if(i % (nx + 1) == nx)
				fprintf(f, "\n");
		}
	}
	free(a);
	free(b);
	free(V);
	free(ja);
	free(ia);
}
int main(){
	// zad 3
	FILE *f;
	f = fopen("matrixA_vectorB.txt", "w");
	metodaPoissona('n', f, 'm');
	fclose(f);
	// zad 5 podpunkt a
	f = fopen("z5_nx_ny_50.txt", "w");
	nx = 50;
	ny = 50;
	metodaPoissona('n', f, 'w');
	fclose(f);
	// zad 5 podpunkt b
	f = fopen("z5_nx_ny_100.txt", "w");
	nx = 100;
	ny = 100;
	metodaPoissona('n', f, 'w');
	fclose(f);
	// zad 5 podpunkt c
	f = fopen("z5_nx_ny_200.txt", "w");
	nx = 200;
	ny = 200;
	metodaPoissona('n', f, 'w');
	fclose(f);
	// zad 6 podpunkt a
	f = fopen("z6_e1_1_e2_1.txt", "w");
	nx = 100;
	ny = 100;
	V1=0;
	V2=0;
	V3=0;
	V4=0;
	metodaPoissona('y', f, 'w');
	fclose(f);
	// zad 6 podpunkt b
	f = fopen("z6_e1_1_e2_2.txt", "w");
	ep2=2;
	metodaPoissona('y', f, 'w');
	fclose(f);
	// zad 6 podpunkt c
	f = fopen("z6_e1_1_e2_10.txt", "w");
	ep2=10;
	metodaPoissona('y', f, 'w');
	fclose(f);
	return 0;
}