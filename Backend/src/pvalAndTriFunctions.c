/**
* This file is entirely taken from the GenEx Project
*
*
* GenExR C Functions
*
* This class was taken from the GenEx Project - written by Rachel Xu
**/

#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rdefines.h>


/**
* pval
* Given vector of correlations r & degrees of freedom d, empty vector p
* (same length as r), calculates the p-values
*
* Reference:
* Code for betai, betacf, gammln; modified code for pval from
* Press, et. al., Numerical Recipes in C: The Art of Scientific Computing, 2nd ed., Cambridge University Press (2002).
* 
* R equivalent by
* Bill Venables, 2000
* https://stat.ethz.ch/pipermail/r-help/2000-January/009758.html
* 
**/

//for pval
#define TINY 1.0e-30

//for betacf
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30


SEXP pval(SEXP r, SEXP d, SEXP p) {
	double betai(double a, double b, double x); //helper function protoype

	//digest R objects
	unsigned int df = asInteger(d); //degrees of freedom
	int i, L = LENGTH(r), PL = LENGTH(p); //length of vectors
	double x, *cor = REAL(r), *pval = REAL(p); //point to real part of Rvector

	//check vector size
	if (PL != L)
		return R_NilValue;

	//compute pval
	for (i = 0; i < L; ++i) {
		//subtract TINY to prevent result from being NaN
		x = ((cor[i])*(cor[i]) * (df)) / ((1 - ((cor[i])*(cor[i]))) - TINY);
		
		if (x > df)
			pval[i] = betai(0.5 * (df), 0.5, df / (df + x));
		else
			pval[i] = 1.0 - betai(0.5, 0.5 * (df), x/(x + df));
	}
	return p;
}

/*Returns the incomplete beta function Ix(a,b).*/
double betai(double a, double b, double x) {
	double betacf(double a, double b, double x);
	double gammln(double xx);
	double bt;

	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		//Factors in front of the continued fraction
		bt = exp(gammln(a+b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0-x));
	if (x < (a + 1.0)/(a + b + 2.0)) //Use continued fraction directly.
		return bt * betacf(a, b, x) / a;
	else	
		//Use continued fraction after symmetry transformation
		return 1.0 - bt * betacf(b, a, 1.0 - x)/b;
}

/*Used by betai: Evaluates continued fraction for incomplete beta function by modiﬁed Lentz’s method (§5.2).*/
double betacf(double a, double b, double x) {
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	//These q’s will be used in factors that occur in the coefficients
	qab = a+b;
	qap = a+1.0;
	qam = a-1.0;

	//First step of Lentz’s method.
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN)
		d=FPMIN;
	d=1.0/d;
	h=d;
	for (m = 1; m <= MAXIT;m++) {
		m2 = 2 * m;
		aa = m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d; //One step (the even one) of the recurrence.
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; //Next step of the recurrence (the odd one).
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break; //Are we done?
	}
	return h;
}


/*Returns the value ln[Γ(xx)] for xx > 0.*/
double gammln(double xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
		return -tmp+log(2.5066282746310005*ser/x);
}



/*****************
*
* Triangular vector functions
*
******************/


/**
* retvec
* Given a vector representative of a symmetric square matrix m; the matrix
* dimension n; and an empty vector of length n * (n - 1) / 2; fills the 
* vector with the lower triangle of the matrix
**/

SEXP retvec(SEXP m, SEXP n, SEXP v) {
	//digest R objects
	int N = asInteger(n);
	int V = LENGTH(v);

	if (isMatrix(m))
		m = coerceVector(m, REALSXP); // does it by column-first orientation

	double *mat = REAL(m); 
	double *vec = REAL(v); 

	int i, j, l = 0; //indices
	int k;

	//check vector size
	if (V != N*(N-1)/2)
		return R_NilValue;

	//iterate and copy
	for (i = 0; i < N; ++i) {
		for (j = i+1; j < N; ++j) {
			k = (i * (N) + j);
			vec[l] = mat[k];
			++l;
		}
	}
	return v;
}

/**
* triIndex
* Given a vector of square matrix index numbers index; the length of a 
* larger vector nold, and a new vector newv of length nold ( nold - 1) / 2;
* fills newv with indices of nold that correspond to a square matrix of
* columns in 'index'
**/

SEXP triIndex(SEXP index, SEXP nold, SEXP newv) {
	//digest from R
	int n = asInteger(nold), *ind = INTEGER(index), *v = INTEGER(newv);
	int L = LENGTH(index), V = LENGTH(newv);
	//indices
	int i, j, k, a, b;

	//check vector size
	if (V != L*(L-1)/2)
		return R_NilValue;

	//generate index
	for (k = 0, i = 0; i < L - 1; ++i) {
		for (j = (i+1); j < L; ++j, ++k) {
			if ((a = ind[i] - 1) <= (b = ind[j] - 1))
				v[k] = (a* (n - 1) - (a-1)*(a)/2 + b - a);
			else
				v[k] = (b* (n - 1) - (b-1)*(b)/2 + a - b);
		}
	}
	return newv;
}


/**
* sumTriCols
* Given triagnular matrix vector vR, side length n of the square matrix, and
* empty vector s (whose length = n), use vR to fill s with the sum of
* absolute values of square matrix columns
**/ 
SEXP sumTriCols(SEXP vR, SEXP nR, SEXP sR) {
	//digest R objects
	int n = asInteger(nR);
	double *v = REAL(vR);
	double *s = REAL(sR);

	//indices
	int i, j, k;

	if (n <= 0)
		return R_NilValue;

	//loop through and sum
	for (k = 0, i = 0; i < n - 1; ++i) {
		for (j = (i + 1); j < n; ++j, ++k) {
			if (!ISNA(v[k])) {
				s[i] += fabs(v[k]);
				s[j] += fabs(v[k]);
			}
		}
	}
	return sR;
}


