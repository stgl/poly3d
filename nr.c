/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: nr.c
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* The functions in this file are adapted by permission from the book:
*
*       Numerical Recipes in C: The Art of Scientific Computing, 2nd Ed.
*       Press, W.H., Teukolsky, S.A., Vetterling, W.T. & Flannery, B.P.
*       1992, Cambridge University Press, Cambridge, 994 p.
*****************************************************************************/


/***************************** Includes/Defines *****************************/
#include <math.h>
#include "nr.h"
#include "nrutil.h"


/*************************** Function: d_jacobi ******************************
* Computes all eigenvalues and eigenvectors of a real symmetric matrix
* a[1..n][1..n].  On output, elements of a above the diagonal are destroyed.
* d[1..n] returns the eigenvalues of a.  v[1..n][1..n] is a matrix whose
* columns contain, on output, the normalized eigenvectors of a.  nrot returns
* the number of Jacobi rotations that were required.
*****************************************************************************/
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void d_jacobi(double **a, int n, double d[], double **v, int *nrot)
#else
void d_jacobi(a, n, d, v, nrot)
double	**a;
int		n;
double	d[];
double	**v;
int		*nrot;
#endif

{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=dvector(1,n);
	z=dvector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_dvector(z,1,n);
			free_dvector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine djacobi");
}
#undef ROTATE


/************************** Function: d_eigsrt ******************************
* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output
* from d_jacobi, this routine sorts the eigenvalues into descending order,
* and rearranges the columns of v correspondingly.  The method is straight
* insertion.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void d_eigsrt(double d[], double **v, int n)
#else
void d_eigsrt(d, v, n)
double	d[];
double	**v;
int		n;
#endif

{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}


/****************************** Function: d_ludcmp ***************************
* Given a matrix a[1..n][1..n], this routine replaces it by the LU
* decomposition of a rowwise permutation of itself.  a and n are input. a is
* output, arranged as in equation (2.3.14) of "Numerical Recipes".  indx[1..n]
* is an output vector that records the row permutation effected by the
* partial pivoting.  d is output as +/- 1 depending on whether the number
* of row interchanges was even or odd, respectively. This routine is used in
* combiniation with d_lubksb to solve linear equations or invert a matrix.
*****************************************************************************/
#define TINY 1.0e-20;
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void d_ludcmp(double **a, int n, int *indx, double *d)
#else
void d_ludcmp(a, n, indx, d)
double	**a;
int		n;
int		*indx;
double	*d;
#endif
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY


/************************* Function: d_lubksb ********************************
* Solves the set of n linear equations Ax = b.  Here a[1..n][1..n] is input,
* not as the matrix A but rather as its LU decomposition, determined by the
* routine d_ludcmp.  indx[1..n] is input as the permutation vector returned
* by d_ludcmp.  b[1..n] is input as the right-hand side vector b, and returns
* with the solution vector x.  a, n, and indx are not modified by this routine
* and can be left in place for successive calls with different right-hand
* sides b.  This routine takes into account the possibility that b will begin
* with many zero elements, so it is efficient for use in matrix inversion.
******************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void d_lubksb(double **a, int n, int *indx, double b[])
#else
void d_lubksb(a, n, indx, b)
double	**a;
int		n;
int		*indx;
double	b[];
#endif
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
