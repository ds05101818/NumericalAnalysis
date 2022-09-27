#include<stdio.h>
#include<math.h>

//rtbis.c
#define JMAX 40
float rtbis(float(*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j, count = 0;
	float dx, f, fmid, xmid, rtb;

	f = (*func)(x1);
	fmid = (*func)(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		count++;
		fmid = (*func)(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0) {
			printf("The number of iterations : %d. ", count);
			return rtb;
		}
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
#undef JMAX

//rtflsp.c
#define MAXIT 30
float rtflsp(float(*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j, count = 0;
	float fl, fh, xl, xh, swap, dx, del, f, rtf;

	fl = (*func)(x1);
	fh = (*func)(x2);
	if (fl*fh > 0.0) nrerror("Root must be bracketed in rtflsp");
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xl = x2;
		xh = x1;
		swap = fl;
		fl = fh;
		fh = swap;
	}
	dx = xh - xl;
	for (j = 1; j <= MAXIT; j++) {
		count++;
		rtf = xl + dx * fl / (fl - fh);
		f = (*func)(rtf);
		if (f < 0.0) {
			del = xl - rtf;
			xl = rtf;
			fl = f;
		}
		else {
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		dx = xh - xl;
		if (fabs(del) < xacc || f == 0.0) {
			printf("The number of iterations : %d. ", count);
			return rtf;
		}
	}
	nrerror("Maximum number of iterations exceeded in rtflsp");
	return 0.0;
}
#undef MAXIT

//rtsec.c
#define MAXIT 30
float rtsec(float(*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j, count = 0;
	float fl, f, dx, swap, xl, rts;

	fl = (*func)(x1);
	f = (*func)(x2);
	if (fabs(fl) < fabs(f)) {
		rts = x1;
		xl = x2;
		swap = fl;
		fl = f;
		f = swap;
	}
	else {
		xl = x1;
		rts = x2;
	}
	for (j = 1; j <= MAXIT; j++) {
		count++;
		dx = (xl - rts)*f / (f - fl);
		xl = rts;
		fl = f;
		rts += dx;
		f = (*func)(rts);
		if (fabs(dx) < xacc || f == 0.0) {
			printf("The number of iterations : %d. ", count);
			return rts;
		}
	}
	nrerror("Maximum number of iterations exceeded in rtsec");
	return 0.0;
}
#undef MAXIT

//rtnewt.c
#define JMAX 20
float rtnewt(void(*funcd)(float, float *, float *), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j, count = 0;
	float df, dx, f, rtn;

	rtn = 0.5*(x1 + x2);
	for (j = 1; j <= JMAX; j++) {
		count++;
		(*funcd)(rtn, &f, &df);
		dx = f / df;
		rtn -= dx;
		if ((x1 - rtn)*(rtn - x2) < 0.0)
			nrerror("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) {
			printf("The number of iterations : %d. ", count);
			return rtn;
		}
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}
#undef JMAX

//rtsafe.c
#define MAXIT 100
float rtsafe(void(*funcd)(float, float *, float *), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j, count = 0;
	float df, dx, dxold, f, fh, fl;
	float temp, xh, xl, rts;

	(*funcd)(x1, &fl, &df);
	(*funcd)(x2, &fh, &df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	rts = 0.5*(x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	(*funcd)(rts, &f, &df);
	for (j = 1; j <= MAXIT; j++) {
		count++;
		if ((((rts - xh)*df - f)*((rts - xl)*df - f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold = dx;
			dx = 0.5*(xh - xl);
			rts = xl + dx;
			if (xl == rts) {
				printf("The number of iterations : %d. ", count);
				return rts;
			}
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts) {
				printf("The number of iterations : %d. ", count);
				return rts;
			}
		}
		if (fabs(dx) < xacc) {
			printf("The number of iterations : %d. ", count);
			return rts;
		}
		(*funcd)(rts, &f, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}
#undef MAXIT

//rtsafe.c for my equation
#define MAXIT 100
float rtsafe1(void(*funcd)(float, float *, float *), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j, count = 0;
	float df, dx, dxold, f, fh, fl;
	float temp, xh, xl, rts;

	(*funcd)(x1, &fl, &df);
	(*funcd)(x2, &fh, &df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	rts = 0.5*(x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	(*funcd)(rts, &f, &df);
	for (j = 1; j <= MAXIT; j++) {
		count++;
		if ((((rts - xh)*df - f)*((rts - xl)*df - f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold = dx;
			dx = 0.5*(xh - xl);
			rts = xl + dx;
			if (xl == rts) {
				//printf("The number of iterations : %d. ", count);
				return rts;
			}
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts) {
				//printf("The number of iterations : %d. ", count);
				return rts;
			}
		}
		if (fabs(dx) < xacc) {
			//printf("The number of iterations : %d. ", count);
			return rts;
		}
		(*funcd)(rts, &f, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}
#undef MAXIT

//muller.c
#define MAXIT 100
float muller(float(*func)(float), float x1, float x2, float xacc) {

	void nrerror(char error_text[]);
	int j, count = 0;
	float p0, p1, p2, p3, f0, f1, f2, a, b, c, d0, d1, d2;

	p0 = x1;
	p1 = x2;
	p2 = (x1 + x2) / 2;

	for (j = 1; j <= MAXIT; j++) {
		count++;

		f0 = (*func)(p0);
		f1 = (*func)(p1);
		f2 = (*func)(p2);

		d0 = p0 - p1;
		d1 = p0 - p2;
		d2 = p1 - p2;

		a = (d2*(f0 - f2) - d1 * (f1 - f2)) / (d0 * d1*d2);
		b = (d1*d1*(f1 - f2) - d2 * d2*(f0 - f2)) / (d0 * d1*d2);
		c = f2;

		if (b > 0) {
			p3 = p2 - (2 * c) / (b + sqrt(b*b - 4 * a*c));
		}
		else {
			p3 = p2 - (2 * c) / (b - sqrt(b*b - 4 * a*c));
		}

		if (fabs(p3 - p2) < xacc*p3) {
			printf("The number of iterations : %d. ", count);
			return p3;
		}

		p0 = p1;
		p1 = p2;
		p2 = p3;
	}
	nrerror("Maximum number of iterations exceeded in muller");
	return 0.0;
}
#undef MAXIT

//equation.c
float equation(float x) {
	return x * x*x - 12 * x*x + 39 * x - 28;
}

float derivative_equation(float x) {
	return 3 * x*x - 24 * x + 39;
}

void solve(float x, float *x1, float *x2) {
	*x1 = equation(x);
	*x2 = derivative_equation(x);
}

//zbrak.c
void zbrak(float(*fx)(float), float x1, float x2, int n, float xb1[], float xb2[], int *nb)
{
	int nbb, i;
	float x, fp, fc, dx;

	nbb = 0;
	dx = (x2 - x1) / n;
	fp = (*fx)(x = x1);
	for (i = 1; i <= n; i++) {
		fc = (*fx)(x += dx);
		if (fc*fp <= 0.0) {
			xb1[++nbb] = x - dx;
			xb2[nbb] = x;
			if (*nb == nbb) return;

		}
		fp = fc;
	}
	*nb = nbb;
}

//bessj0.c
float bessj0(float x)
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
			+ y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
			+ y * (59272.64853 + y * (267.8532712 + y * 1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
			+ y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
			+ y * (-0.6911147651e-5 + y * (0.7621095161e-6
				- y * 0.934945152e-7)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z * sin(xx)*ans2);
	}
	return ans;
}

//bessj1.c
float bessj1(float x)
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
			+ y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
		ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
			+ y * (99447.43394 + y * (376.9991397 + y * 1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
			+ y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		ans2 = 0.04687499995 + y * (-0.2002690873e-3
			+ y * (0.8449199096e-5 + y * (-0.88228987e-6
				+ y * 0.105787412e-6)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z * sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

void func(float x, float *x1, float *x2) {
	*x1 = bessj0(x);
	*x2 = -bessj1(x);
}

//nrutil.c
/* CAUTION: This is the traditional K&R C (only) version of the Numerical
Recipes utility file nrutil.c.  Do not confuse this file with the
same-named file nrutil.c that is supplied in the same subdirectory or
archive as the header file nrutil.h.  *That* file contains both ANSI and
traditional K&R versions, along with #ifdef macros to select the
correct version.  *This* file contains only traditional K&R.           */
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void exit();

	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}

float *vector(nl, nh)
long nh, nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v = (float *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v - nl + NR_END;
}

int *ivector(nl, nh)
long nh, nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v = (int *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v - nl + NR_END;
}

unsigned char *cvector(nl, nh)
long nh, nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v = (unsigned char *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v - nl + NR_END;
}

unsigned long *lvector(nl, nh)
long nh, nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v = (unsigned long *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v - nl + NR_END;
}

double *dvector(nl, nh)
long nh, nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v = (double *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v - nl + NR_END;
}

float **matrix(nrl, nrh, ncl, nch)
long nch, ncl, nrh, nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	/* allocate pointers to rows */
	m = (float **)malloc((unsigned int)((nrow + NR_END) * sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (float *)malloc((unsigned int)((nrow*ncol + NR_END) * sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl, nrh, ncl, nch)
long nch, ncl, nrh, nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	double **m;

	/* allocate pointers to rows */
	m = (double **)malloc((unsigned int)((nrow + NR_END) * sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (double *)malloc((unsigned int)((nrow*ncol + NR_END) * sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl, nrh, ncl, nch)
long nch, ncl, nrh, nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int **m;

	/* allocate pointers to rows */
	m = (int **)malloc((unsigned int)((nrow + NR_END) * sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl] = (int *)malloc((unsigned int)((nrow*ncol + NR_END) * sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a, oldrl, oldrh, oldcl, oldch, newrl, newcl)
float **a;
long newcl, newrl, oldch, oldcl, oldrh, oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
	float **m;

	/* allocate array of pointers to rows */
	m = (float **)malloc((unsigned int)((nrow + NR_END) * sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for (i = oldrl, j = newrl; i <= oldrh; i++, j++) m[j] = a[i] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a, nrl, nrh, ncl, nch)
float *a;
long nch, ncl, nrh, nrl;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;

	/* allocate pointers to rows */
	m = (float **)malloc((unsigned int)((nrow + NR_END) * sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl] = a - ncl;
	for (i = 1, j = nrl + 1; i<nrow; i++, j++) m[j] = m[j - 1] + ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(nrl, nrh, ncl, nch, ndl, ndh)
long nch, ncl, ndh, ndl, nrh, nrl;
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t = (float ***)malloc((unsigned int)((nrow + NR_END) * sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl] = (float **)malloc((unsigned int)((nrow*ncol + NR_END) * sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl] = (float *)malloc((unsigned int)((nrow*ncol*ndep + NR_END) * sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j - 1] + ndep;
	for (i = nrl + 1; i <= nrh; i++) {
		t[i] = t[i - 1] + ncol;
		t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
		for (j = ncl + 1; j <= nch; j++) t[i][j] = t[i][j - 1] + ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v, nl, nh)
float *v;
long nh, nl;
/* free a float vector allocated with vector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_ivector(v, nl, nh)
int *v;
long nh, nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_cvector(v, nl, nh)
long nh, nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_lvector(v, nl, nh)
long nh, nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_dvector(v, nl, nh)
double *v;
long nh, nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_matrix(m, nrl, nrh, ncl, nch)
float **m;
long nch, ncl, nrh, nrl;
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void free_dmatrix(m, nrl, nrh, ncl, nch)
double **m;
long nch, ncl, nrh, nrl;
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void free_imatrix(m, nrl, nrh, ncl, nch)
int **m;
long nch, ncl, nrh, nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void free_submatrix(b, nrl, nrh, ncl, nch)
float **b;
long nch, ncl, nrh, nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG)(b + nrl - NR_END));
}

void free_convert_matrix(b, nrl, nrh, ncl, nch)
float **b;
long nch, ncl, nrh, nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG)(b + nrl - NR_END));
}

void free_f3tensor(t, nrl, nrh, ncl, nch, ndl, ndh)
float ***t;
long nch, ncl, ndh, ndl, nrh, nrl;
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG)(t[nrl] + ncl - NR_END));
	free((FREE_ARG)(t + nrl - NR_END));
}

int main() {
	float x1 = 1.0;
	float x2 = 10.0;
	float xacc = 0.000001;
	float xb1[10], xb2[10], root;
	int nb = 10;

	zbrak(bessj0, x1, x2, 9, xb1, xb2, &nb);
	printf("Bisection Method\n");
	for (int i = 1; i <= 3; i++) {
		root = rtbis(bessj0, xb1[i], xb2[i], xacc);
		printf(" The answer of %dth root : %f\n", i, root);
	}
	
	printf("\nLinear Interpolation Method\n");
	for (int i = 1; i <= 3; i++) {
		root = rtflsp(bessj0, xb1[i], xb2[i], xacc);
		printf(" The answer of %dth root : %f\n", i, root);
	}
	
	printf("\nSecant Method\n");
	for (int i = 1; i <= 3; i++) {
		root = rtsec(bessj0, xb1[i], xb2[i], xacc);
		printf(" The answer of %dth root : %f\n", i, root);
	}
	
	printf("\nNewton-Phapson Method\n");
	for (int i = 1; i <= 3; i++) {
		root = rtnewt(func, xb1[i], xb2[i], xacc);
		printf(" The answer of %dth root : %f\n", i, root);
	}
	
	printf("\nNewton with Bracketing Method\n");
	for (int i = 1; i <= 3; i++) {
		root = rtsafe(func, xb1[i], xb2[i], xacc);
		printf(" The answer of %dth root : %f\n", i, root);
	}
	
	printf("\nMuller Method\n");
	for (int i = 1; i <= 3; i++) {
		root = muller(bessj0, xb1[i], xb2[i], xacc);
		printf(" The answer of %dth root : %f\n", i, root);
	}

	printf("\n");
	zbrak(equation, x1, x2, 70, xb1, xb2, &nb);
	printf("Nonlinear Equation Using Newton with Bracketing Method\n");
	for (int i = 1; i <= 3; i++) {
		root = rtsafe(solve, xb1[i], xb2[i], xacc);
		printf("The answer of %dth root : %f\n", i, root);
	}
}