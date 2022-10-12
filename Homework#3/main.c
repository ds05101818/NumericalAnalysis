#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#define NR_END 1
#define FREE_ARG char*

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static float maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}

int* ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int* v;

	v = (int*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v - nl + NR_END;
}

void free_ivector(int* v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

float* vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float* v;

	v = (float*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v - nl + NR_END;
}

void free_vector(float* v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

//pythag.c
float pythag(float a, float b)
{
	float absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

//gaussj.c
void gaussj(a, n, b, m)
float** a, ** b;
int m, n;
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

//ludcmp.c
void ludcmp(a, n, indx, d)
float** a, *d;
int* indx, n;
{
	int i, imax, j, k;
	float big, dum, sum, temp;
	float* vv;

	vv = vector(1, n);
	*d = 1.0;
	for (i = 1; i <= n; i++) {
		big = 0.0;
		for (j = 1; j <= n; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i] = 1.0 / big;
	}
	for (j = 1; j <= n; j++) {
		for (i = 1; i < j; i++) {
			sum = a[i][j];
			for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i <= n; i++) {
			sum = a[i][j];
			for (k = 1; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 1; k <= n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i <= n; i++) a[i][j] *= dum;
		}
	}
	free_vector(vv, 1, n);
}

//svdcmp.c
void svdcmp(a, m, n, w, v)
float** a, ** v, w[];
int m, n;
{
	float pythag(float a, float b);
	int flag, i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z, *rv1;

	rv1 = vector(1, n);
	g = scale = anorm = 0.0;
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
				}
				for (k = i; k <= m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
				}
				for (k = l; k <= n; k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--) {
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++)
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n); i >= 1; i--) {
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++) a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j <= n; j++) {
				for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
				f = (s / a[i][i]) * g;
				for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
			}
			for (j = i; j <= m; j++) a[j][i] *= g;
		}
		else for (j = i; j <= m; j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k = n; k >= 1; k--) {
		for (its = 1; its <= 30; its++) {
			flag = 1;
			for (l = k; l >= 1; l--) {
				nm = l - 1;
				if ((float)(fabs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((float)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((float)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_vector(rv1, 1, n);
}

//mprove.c
void mprove(a, alud, n, indx, b, x)
float** a, ** alud, b[], x[];
int indx[], n;
{
	void lubksb();
	int j, i;
	double sdp;
	float* r;

	r = vector(1, n);
	for (i = 1; i <= n; i++) {
		sdp = -b[i];
		for (j = 1; j <= n; j++) sdp += a[i][j] * x[j];
		r[i] = sdp;
	}
	lubksb(alud, n, indx, r);
	for (i = 1; i <= n; i++) x[i] -= r[i];
	free_vector(r, 1, n);
}

//lubksb.c
void lubksb(a, n, indx, b)
float** a, b[];
int* indx, n;
{
	int i, ii = 0, ip, j;
	float sum;

	for (i = 1; i <= n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
		else if (sum) ii = i;
		b[i] = sum;
	}
	for (i = n; i >= 1; i--) {
		sum = b[i];
		for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}

//svbksb.c
void svbksb(u, w, v, m, n, b, x)
float** u, ** v, b[], w[], x[];
int m, n;
{
	int jj, j, i;
	float s, *tmp;

	tmp = vector(1, n);
	for (j = 1; j <= n; j++) {
		s = 0.0;
		if (w[j]) {
			for (i = 1; i <= m; i++) s += u[i][j] * b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j = 1; j <= n; j++) {
		s = 0.0;
		for (jj = 1; jj <= n; jj++) s += v[j][jj] * tmp[jj];
		x[j] = s;
	}
	free_vector(tmp, 1, n);
}

void gaussj_solution(FILE* fp, float** a, float** b, int* n) {
	int m;
	fscanf(fp, "%d %d", n, &m);

	a = (float**)malloc(sizeof(float*) * (m + 1));
	b = (float**)malloc(sizeof(float*) * (m + 1));
	for (int i = 1; i <= m; i++) {
		a[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	for (int i = 1; i <= m; i++) {
		b[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= m; j++) {
			fscanf(fp, "%f", &a[i][j]);
		}
	}
	for (int i = 1; i <= m; i++) {
		fscanf(fp, "%f", &b[i][1]);
	}
	printf("The answer by using gaussj\n");
	gaussj(a, *n, b, 1);
	for (int i = 1; i <= *n; i++) {
		printf("%f ", b[i][1]);
	}
	printf("\n");
	for (int i = 1; i <= *n; i++) {
		free(a[i]);
		free(b[i]);
	}
	free(a);
	free(b);
}

void ludcmp_solution(FILE* fp, float** a, float* b, int* n) {
	int m;
	fscanf(fp, "%d %d", n, &m);
	a = (float**)malloc(sizeof(float*) * (m + 1));
	b = (float*)malloc(sizeof(float) * (m + 1));
	for (int i = 1; i <= m; i++) {
		a[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= m; j++) {
			fscanf(fp, "%f", &a[i][j]);
		}
	}
	for (int i = 1; i <= m; i++) {
		fscanf(fp, "%f", &b[i]);
	}
	printf("The answer by using ludcmp\n");
	int* indx = (int*)malloc(sizeof(int) * (*n + 1));
	float d;
	ludcmp(a, *n, indx, &d);
	lubksb(a, *n, indx, b);
	for (int i = 1; i <= *n; i++) {
		printf("%f ", b[i]);
	}
	printf("\n");
	for (int i = 1; i <= *n; i++) {
		free(a[i]);
	}
	free(a);
	free(b);
	free(indx);
}

void svdcmp_solution(FILE* fp, float** a, float* b, int* n) {
	int m;
	fscanf(fp, "%d %d", n, &m);
	a = (float**)malloc(sizeof(float*) * (m + 1));
	b = (float*)malloc(sizeof(float) * (m + 1));
	for (int i = 1; i <= m; i++) {
		a[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= m; j++) {
			fscanf(fp, "%f", &a[i][j]);
		}
	}
	for (int i = 1; i <= m; i++) {
		fscanf(fp, "%f", &b[i]);
	}
	printf("The answer by using svdcmp\n");
	float* w = (float*)malloc(sizeof(float) * (m + 1));
	float* x = (float*)malloc(sizeof(float) * (m + 1));
	float** v = (float**)malloc(sizeof(float*) * (m + 1));
	for (int i = 1; i <= *n; i++) {
		v[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	svdcmp(a, *n, *n, w, v);
	svbksb(a, w, v, *n, *n, b, x);
	for (int i = 1; i <= *n; i++) {
		printf("%f ", x[i]);
	}
	printf("\n");
	for (int i = 1; i <= *n; i++) {
		free(a[i]);
		free(v[i]);
	}
	free(a);
	free(b);
	free(w);
	free(x);
	free(v);
}

void itmprove(FILE *fp, float **a, float *b, int *n) {
	int m;
	fscanf(fp, "%d %d", n, &m);
	a = (float**)malloc(sizeof(float*) * (m + 1));
	b = (float*)malloc(sizeof(float) * (m + 1));
	for (int i = 1; i <= m; i++) {
		a[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= m; j++) {
			fscanf(fp, "%f", &a[i][j]);
		}
	}
	for (int i = 1; i <= m; i++) {
		fscanf(fp, "%f", &b[i]);
	}
	printf("The answer by using mprove,,,\n");

	float** alud = (float**)malloc(sizeof(float*) * (*n + 1));
	float* x = (float*)malloc(sizeof(float) * (*n + 1));
	int* indx = (int*)malloc(sizeof(int) * (*n + 1));
	float d;
	for (int i = 1; i <= *n; i++) {
		alud[i] = (float*)malloc(sizeof(float) * (*n + 1));
	}
	for (int i = 1; i <= *n; i++) {
		for (int j = 1; j <= *n; j++) {
			alud[i][j] = a[i][j];
		}
	}
	ludcmp(alud, *n, indx, &d);
	mprove(a, alud, *n, indx, b, x);
	for (int i = 1; i <= *n; i++) {
		printf("%f ", x[i]);
	}
	printf("\n");
	for (int i = 1; i <= *n; i++) {
		free(alud[i]);
		free(a[i]);
	}
	free(alud);
	free(a);
	free(b);
	free(indx);
	free(x);
}

void inverse(FILE* fp, float** a, int* n) {
	printf("The inverse is,,,\n");
	int m;
	fscanf(fp, "%d %d", n, &m);
	a = (float**)malloc(sizeof(float*) * (m + 1));
	float** y = (float**)malloc(sizeof(float*) * (m + 1));
	for (int i = 1; i <= m; i++) {
		a[i] = (float*)malloc(sizeof(float) * (m + 1));
		y[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	for (int i = 1; i <= *n; i++) {
		for (int j = 1; j <= *n; j++) {
			fscanf(fp, "%f", &a[i][j]);
		}
	}
	float* col = (float*)malloc(sizeof(float) * (m + 1));
	int* indx = (int*)malloc(sizeof(int) * (m + 1));
	float d;

	ludcmp(a, *n, indx, &d);
	for (int j = 1; j <= *n; j++) {
		for (int i = 1; i <= *n; i++) {
			col[i] = 0.0;
		}
		col[j] = 1.0;
		lubksb(a, *n, indx, col);
		for (int i = 1; i <= *n; i++) {
			y[i][j] = col[i];
		}
	}
	for (int i = 1; i <= *n; i++) {
		for (int j = 1; j <= *n; j++) {
			printf("%f ", y[i][j]);
		}
		printf("\n");
	}
	for (int i = 1; i <= *n; i++) {
		free(a[i]);
		free(y[i]);
	}
	free(a);
	free(y);
	free(col);
	free(indx);
}

void determinant(FILE *fp, float **a, int *n) {
	int m;
	fscanf(fp, "%d %d", n, &m);
	a = (float**)malloc(sizeof(float*) * (m + 1));
	for (int i = 1; i <= m; i++) {
		a[i] = (float*)malloc(sizeof(float) * (m + 1));
	}
	int* indx = (int*)malloc(sizeof(int) * (m + 1));
	float det;
	for (int i = 1; i <= *n; i++) {
		for (int j = 1; j <= *n; j++) {
			fscanf(fp, "%f", &a[i][j]);
		}
	}
	ludcmp(a, *n, indx, &det);
	for (int j = 1; j <= *n; j++) {
		det *= a[j][j];
	}
	printf("The result of determinant\n");
	printf("%f\n", det);
	for (int i = 1; i <= *n; i++) {
		free(a[i]);
	}
	free(a);
	free(indx);
}

int main() {
	FILE* fp;
	int n;
	float** a = NULL;
	float* b = NULL;
	float** gauss_b = NULL;

	// lineq1.dat
	printf("\nThe Answer of Linear Equation1\n\n");

	//gaussj_solution
	fp = fopen("lineq1.dat", "r");
	gaussj_solution(fp, a, gauss_b, &n);
	fclose(fp);
	printf("\n");

	//ludcmp_solution
	fp = fopen("lineq1.dat", "r");
	ludcmp_solution(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//svdcmp_solution
	fp = fopen("lineq1.dat", "r");
	svdcmp_solution(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//itmprove
	fp = fopen("lineq1.dat", "r");
	itmprove(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//no inverse

	//determinant
	fp = fopen("lineq1.dat", "r");
	determinant(fp, a, &n);
	printf("The determinant of linear equation1 is 0.\n So, the inverse of linear equation1 does not exist.");
	fclose(fp);
	printf("\n");

	// lineq2.dat
	printf("\n\nThe Answer of Linear Equation2\n\n");

	//gaussj_solution
	fp = fopen("lineq2.dat", "r");
	gaussj_solution(fp, a, gauss_b, &n);
	fclose(fp);
	printf("\n");

	//ludcmp_solution
	fp = fopen("lineq2.dat", "r");
	ludcmp_solution(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//svdcmp_solution
	fp = fopen("lineq2.dat", "r");
	svdcmp_solution(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//itmprove
	fp = fopen("lineq2.dat", "r");
	itmprove(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//inverse
	fp = fopen("lineq2.dat", "r");
	inverse(fp, a, &n);
	fclose(fp);
	printf("\n");

	//determinant
	fp = fopen("lineq2.dat", "r");
	determinant(fp, a, &n);
	fclose(fp);
	printf("\n");

	// lineq3.dat
	printf("\n\nThe Answer of Linear Equation3\n\n");

	//gaussj_solution
	fp = fopen("lineq3.dat", "r");
	gaussj_solution(fp, a, gauss_b, &n);
	fclose(fp);
	printf("\n");

	//ludcmp_solution
	fp = fopen("lineq3.dat", "r");
	ludcmp_solution(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//svdcmp_solution
	fp = fopen("lineq3.dat", "r");
	svdcmp_solution(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//itmprove
	fp = fopen("lineq3.dat", "r");
	itmprove(fp, a, b, &n);
	fclose(fp);
	printf("\n");

	//inverse
	fp = fopen("lineq3.dat", "r");
	inverse(fp, a, &n);
	fclose(fp);
	printf("\n");

	//determinant
	fp = fopen("lineq3.dat", "r");
	determinant(fp, a, &n);
	fclose(fp);
	printf("\n");
}
#undef SWAP
#undef TINY
