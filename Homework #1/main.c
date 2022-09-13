#include<stdio.h>
#include<math.h>

void machar(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, float *eps, float *epsneg,
	float *xmin, float *xmax);

void machar_double(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, double *eps, double *epsneg,
	double *xmin, double *xmax);

float get_eps_float() {
	int count = 1;
	float n = 1;
	while (1 + n / 2 != 1) {
		n /= 2;
		count++;
	}
	printf("The result of float count : %d\n", count);
	return n;
}

double get_eps_double() {
	int count = 1;
	double n = 1;
	while (1 + n / 2 != 1) {
		n /= 2;
		count++;
	}
	printf("The result of double count : %d\n", count);
	return n;
}

int main() {
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	float f_eps, f_epsneg, f_xmin, f_xmax;
	double d_eps, d_epsneg, d_xmin, d_xmax;

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &f_eps, &f_epsneg, &f_xmin, &f_xmax);
	machar_double(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &d_eps, &d_epsneg, &d_xmin, &d_xmax);

	printf("Method 1\n");
	printf("The machine accuracy of float : %12.6g\n", f_eps);
	printf("The machine accuracy of double : %12.6g\n", d_eps);
	printf("\n");
	printf("Method 2\n");
	printf("The machine accuracy of get_eps_float : %12.6g\n", get_eps_float());
	printf("The machine accuracy of get_eps_double : %12.6g\n", get_eps_double());
}
