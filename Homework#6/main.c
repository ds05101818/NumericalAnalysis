#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*

int* ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int* v;

	v = (int*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
	//if (!v) nrerror("allocation failure in ivector()");
	return v - nl + NR_END;
}

void free_ivector(int* v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	/* allocate pointers to rows */
	m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float*)));
	//if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (float *)malloc((size_t)((nrow*ncol + NR_END) * sizeof(float)));
	//if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void gaussj(float **a, int n, float **b, int m)
{
	int* indxc, *indxr, *ipiv;
	int i, icol, irow, j, k, l, ll;
	float big, dum, pivinv, temp;

	indxc = ivector(1, n);
	indxr = ivector(1, n);
	ipiv = ivector(1, n);
	for (j = 1; j <= n; j++) ipiv[j] = 0;
	for (i = 1; i <= n; i++) {
		big = 0.0;
		for (j = 1; j <= n; j++)
			if (ipiv[j] != 1)
				for (k = 1; k <= n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
				for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) {
			//nrerror("gaussj: Singular Matrix");
			printf("nerror,,, Gaussj: Singular Matrix,,,\n");
			return;
		}
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
		for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
		for (ll = 1; ll <= n; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
				for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
			}
	}
	for (l = n; l >= 1; l--) {
		if (indxr[l] != indxc[l])
			for (k = 1; k <= n; k++)
				SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	}
	free_ivector(ipiv, 1, n);
	free_ivector(indxr, 1, n);
	free_ivector(indxc, 1, n);
}

void read_data(FILE *file, float *x, float *y, float *xp, float *yp) {
	for (int i = 0; i < 77; i++) {
		fscanf(file, "%f %f %f %f", &x[i], &y[i], &xp[i], &yp[i]);
	}
}

void calculation(float **Jy, float *x, float *y, float *xp, float *yp, int n, int m, int k) {
	float ** JJ = matrix(1, m, 1, m);
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= m; j++) {
			JJ[i][j] = 0;
			if (j <= 2) {
				Jy[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		JJ[1][1] += x[i] * x[i];
		JJ[1][2] += x[i] * y[i];
		JJ[1][3] += x[i];
		JJ[2][2] += y[i] * y[i];
		JJ[2][3] += y[i];
		Jy[1][1] += x[i] * xp[i];
		Jy[1][2] += x[i] * yp[i];
		Jy[2][1] += xp[i] * y[i];
		Jy[2][2] += y[i] * yp[i];
		Jy[3][1] += xp[i];
		Jy[3][2] += yp[i];
	}
	JJ[2][1] = JJ[1][2];
	JJ[3][1] = JJ[1][3];
	JJ[3][2] = JJ[2][3];
	JJ[3][3] = n;

	gaussj(JJ, m, Jy, k);
}

int main() {
	FILE *fp;
	float *x = (float*)malloc(sizeof(float)*77);
	float *y = (float*)malloc(sizeof(float)*77);
	float *xp = (float*)malloc(sizeof(float)*77);
	float *yp = (float*)malloc(sizeof(float)*77);
	float **Jy = matrix(1, 3, 1, 2);
	int n = 1;

	fp = fopen("fitdata1.dat", "r"); 
	printf("The Answer of fitdata1\n"); 
	read_data(fp, x, y, xp, yp);
	calculation(Jy, x, y, xp, yp, 77, 3, 2);
	for (int i = 1; i <= 2; i++) {
		for (int j = 1; j <= 3; j++) {
			printf("a%d : %f\n", n++, Jy[j][i]);
		}
	}
	fclose(fp);
	printf("\n");

	fp = fopen("fitdata2.dat", "r");
	printf("The Answer of fitdata2\n");
	read_data(fp, x, y, xp, yp);
	calculation(Jy, x, y, xp, yp, 77, 3, 2);
	n = 1;
	for (int i = 1; i <= 2; i++) {
		for (int j = 1; j <= 3; j++) {
			printf("a%d : %f\n", n++, Jy[j][i]);
		}
	}
	fclose(fp);
	printf("\n");

	fp = fopen("fitdata3.dat", "r");
	printf("The Answer of fitdata3\n");
	read_data(fp, x, y, xp, yp);
	calculation(Jy, x, y, xp, yp, 77, 3, 2);
	n = 1;
	for (int i = 1; i <= 2; i++) {
		for (int j = 1; j <= 3; j++) {
			printf("a%d : %f\n", n++, Jy[j][i]);
		}
	}
	fclose(fp);
	printf("\n");
}