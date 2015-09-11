/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: infcoeff.c
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* The following functions calculate displacement and strain influence
* coefficients for the angular dislocation as described by Comninou &
* Dunders (1975).  The contributions from the dislocation and its image
* have been split into separate equations, so that whole- and half-space
* problems may be addressed with a single set of functions.
*****************************************************************************/


/***************************** Includes/Defines *****************************/
#include <math.h>
#include "infcoeff.h"
#include "safetan.h"
#include "pi.h"


/*************************** Function Declarations *************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
double	reduce_angle(double angle);
#else
double	reduce_angle();
#endif


/*********************** Function: comninou_displ_ics ***********************
* Calculates displacement influence coefficients for the angular dislocation.
*
* In:	y				- coords of point at which to calc inf coeffs
*		a				- depth to angular dislocation vertex
*		beta			- angle subtended by dislocation
*		pr				- poisson's ratio
*		half_space		- half/whole-space flag
* Out:	displ_ic[i][j]	- jth component of displ at y due to unit ith Burgers
*						  vector component
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     comninou_displ_ics(double y[3], double a, double beta, double pr,int	 half_space, double displ_ic[3][3])
#else
void comninou_displ_ics(y, a, beta, pr, half_space, displ_ic)
double y[3];
double a;
double beta;
double pr;
int	 half_space;
double displ_ic[3][3];
#endif
{
	/*double trash1[3][3];*/
	/*double trash2[3][3];*/
	double dic2[3][3];
	double term1, term2, term3;

	int		i, j;
	double	displ_ic_ft[3][3];

	double	y1		= y[0];
	double	y2		= y[1];
	double	y3		= y[2];
	double	sinb	= sin(beta);
	double	cosb	= cos(beta);
	double	cotb	= 1/tan(beta);

	double	y3b		= y3 + 2*a;
	double	z1		= y1*cosb - y3*sinb;
	double	z3		= y1*sinb + y3*cosb;
	double	r		= sqrt(y1*y1 + y2*y2 + y3*y3);
	double	z1b		= y1*cosb + y3b*sinb;
	double	z3b		= -y1*sinb + y3b*cosb;
	double	rb		= sqrt(y1*y1 + y2*y2 + y3b*y3b);

	/* For debugging only
	---------------------*/
	/*
	double	f		= reduce_angle(-safe_atan2(y2,y1) + safe_atan2(y2,z1) +
					  safe_atan2(y2*r*sinb,(y1*z1+(y2*y2)*cosb)));
	double	fb		= reduce_angle(-safe_atan2(y2,y1) + safe_atan2(y2,z1b) + 
					  safe_atan2(y2*rb*sinb,(y1*z1b+(y2*y2)*cosb)));
	*/
	double	f		= -safe_atan2(y2,y1) + safe_atan2(y2,z1) +
					  safe_atan2(y2*r*sinb,(y1*z1+(y2*y2)*cosb));
	double	fb		= -safe_atan2(y2,y1) + safe_atan2(y2,z1b) + 
					  safe_atan2(y2*rb*sinb,(y1*z1b+(y2*y2)*cosb));


	/* Dislocation-Induced Displacements
	------------------------------------*/
		
	displ_ic[0][0] = 
		-y1*y2*(1/(r*(r-y3))) - 
		y2*cosb*((r*sinb-y1)/(r*(r-z3)));
		
	displ_ic_ft[0][0] = 2*(1-pr)*(f);
		
	displ_ic[0][1] = 
		(1-2*pr)*(log(r-y3) - cosb * (log(r-z3))) -
		(y2*y2)*(1/(r*(r-y3)) - cosb*(1/(r*(r-z3))));
		
	displ_ic_ft[0][1] = 0.0;
		
	displ_ic[0][2] = 
		y2 * (1/r - cosb * ((r*cosb - y3)/(r*(r-z3))));
		
	displ_ic_ft[0][2] = 0.0;
		
	displ_ic[1][0] = 
		(-(1-2*pr))*(log(r-y3)-cosb*(log(r-z3))) +
		(y1*y1)*(1/(r*(r-y3))) + z1*(r*sinb-y1)/(r*(r-z3));
		
	displ_ic_ft[1][0] = 0.0;
		
	displ_ic[1][1] = 
		y1*y2*(1/(r*(r-y3))) - y2*(z1/(r*(r-z3)));
		
	displ_ic_ft[1][1] = 2*(1-pr)*(f);
		
	displ_ic[1][2] = 
		-(1-2*pr)*sinb*(log(r-z3)) -
		y1*(1/r) + z1*(r*cosb-y3)/(r*(r-z3));
		
	displ_ic_ft[1][2] = 0.0;
		
	displ_ic[2][0] = 
		y2*sinb*((r*sinb-y1)/(r*(r-z3)));
		
	displ_ic_ft[2][0] = 0.0;
		
	displ_ic[2][1] = 
		(1-2*pr)*sinb*(log(r-z3)) -
		(y2*y2)*sinb*(1/(r*(r-z3)));
		
	displ_ic_ft[2][1] = 0.0;
		
	displ_ic[2][2] = 
		y2*sinb*((r*cosb-y3)/(r*(r-z3)));
		
	displ_ic_ft[2][2] = 2*(1-pr)*(f);
		

	if (half_space) {
		
		/* Image-Induced Displacements
		------------------------------*/
	
		displ_ic[0][0] += 
			-y1*y2*(1/(rb*(rb+y3b))) -
			y2*cosb*((rb*sinb-y1)/(rb*(rb+z3b)));
			
		displ_ic_ft[0][0] += 2*(1-pr)*(fb);
			
		displ_ic[0][1] += 
			(1-2*pr)*(log(rb+y3b) - cosb * (log(rb+z3b))) -
			(y2*y2)*(1/(rb*(rb+y3b)) - cosb*(1/(rb*(rb+z3b))));
			
		displ_ic_ft[0][1] += 0.0;
			
		displ_ic[0][2] += 
			y2 * (-1/rb - cosb * (-(rb*cosb + y3b)/(rb*(rb+z3b))));
			
		displ_ic_ft[0][2] += 0.0;
			
		displ_ic[1][0] += 
			(-(1-2*pr))*(log(rb+y3b)-cosb*(log(rb+z3b))) +
			(y1*y1)*(1/(rb*(rb+y3b))) + z1b*(rb*sinb-y1)/(rb*(rb+z3b));
			
		displ_ic_ft[1][0] += 0.0;
			
		displ_ic[1][1] += 
			y1*y2*(1/(rb*(rb+y3b))) - y2*(z1b/(rb*(rb+z3b)));
			
		displ_ic_ft[1][1] += 2*(1-pr)*(fb);
			
		displ_ic[1][2] += 
			-(1-2*pr)*sinb*(-log(rb+z3b)) -
			y1*(-1/rb) - z1b*(rb*cosb+y3b)/(rb*(rb+z3b));
			
		displ_ic_ft[1][2] += 0.0;
			
		displ_ic[2][0] += 
			y2*sinb*((rb*sinb - y1)/(rb*(rb+z3b)));
			
		displ_ic_ft[2][0] += 0.0;
			
		displ_ic[2][1] += 
			(1-2*pr)*sinb*(log(rb+z3b)) -
			(y2*y2)*sinb*(1/(rb*(rb+z3b)));
			
		displ_ic_ft[2][1] += 0.0;
			
		displ_ic[2][2] += 
			y2*sinb*(-(rb*cosb+y3b)/(rb*(rb+z3b)));
			
		displ_ic_ft[2][2] += 2*(1-pr)*(-fb);
			
	
		/* Corrective Displacements
		---------------------------*/
			
		displ_ic[0][0] += 2.0*(
			/* term1 */
			(1-2*pr)*y2/(rb+y3b) *
			((1-2*pr-a/rb)*cotb - y1/(rb+y3b)*(pr+a/rb)) +
			/* term2 */
			(1-2*pr)*y2*cosb*cotb/(rb+z3b)*(cosb+a/rb) +
			/* term3 */
			a*y2*(y3b-a)*cotb/pow(rb,3.0) +
			/* term4 */
			y2*(y3b-a)/(rb*(rb+y3b))*(-(1-2*pr)*cotb + y1/(rb+y3b) *
			(2*pr+a/rb) + a*y1/(rb*rb)) +
			/* term5 */
			y2*(y3b-a)/(rb*(rb+z3b))*(cosb/(rb+z3b)*((rb*cosb+y3b) *
			((1-2*pr)*cosb-a/rb)*cotb + 
			2*(1-pr)*(rb*sinb-y1)*cosb) -
			a*y3b*cosb*cotb/(rb*rb))
			);
			
		displ_ic_ft[0][0] += 2* -2*(1-pr)*(1-2*pr)*fb*(cotb*cotb);
			
		displ_ic[0][1] += 2.0*(
			/* term1 */
			(1-2*pr)*((2*(1-pr)*(cotb*cotb)-pr)*log(rb+y3b) -
			(2*(1-pr)*(cotb*cotb)+1-2*pr)*cosb*log(rb+z3b)) -
			/* term2 */
			(1-2*pr)/(rb+y3b)*(y1*cotb*(1-2*pr-a/rb) + pr*y3b - a +
			(y2*y2)/(rb+y3b)*(pr+a/rb)) -
			/* term3 */
			(1-2*pr)*z1b*cotb/(rb+z3b)*(cosb+a/rb) - 
			a*y1*(y3b-a)*cotb/(pow(rb,3.0)) +
			/* term4 */
			(y3b-a)/(rb+y3b)*(-2*pr + 1/rb*((1-2*pr)*y1*cotb-a) +
			(y2*y2)/(rb*(rb+y3b))*(2*pr+a/rb)+a*(y2*y2)/pow(rb,3.0)) +
			/* term5 */
			(y3b-a)/(rb+z3b)*((cosb*cosb) - 
			1/rb*((1-2*pr)*z1b*cotb+a*cosb) +
			a*y3b*z1b*cotb/pow(rb,3.0) - 1/(rb*(rb+z3b)) *
			((y2*y2)*(cosb*cosb) - a*z1b*cotb/rb*(rb*cosb+y3b)))
			);
			
		displ_ic_ft[0][1] += 0.0;
			
		displ_ic[0][2] += 2.0*(
			/* term1 */
			2*(1-pr)*(y2/(rb+y3b)*(2*pr+a/rb) -
			y2*cosb/(rb+z3b)*(cosb+a/rb)) +
			/* term2 */
			y2*(y3b-a)/rb*(2*pr/(rb+y3b)+a/(rb*rb)) +
			/* term3 */
			y2*(y3b-a)*cosb/(rb*(rb+z3b))*
			(1-2*pr-(rb*cosb+y3b)/(rb+z3b) *
			(cosb + a/rb) - a*y3b/(rb*rb))
			);
			
		displ_ic_ft[0][2] += 2* 2*(1-pr)*((1-2*pr)*fb*cotb);
			
		displ_ic[1][0] += 2.0*(
			/* term1 */
			(1-2*pr)*((2*(1-pr)*(cotb*cotb)+pr)*log(rb+y3b) -
			(2*(1-pr)*(cotb*cotb)+1)*cosb*log(rb+z3b)) +
			/* term2 */
			(1-2*pr)/(rb+y3b)*
			(-(1-2*pr)*y1*cotb+pr*y3b-a+a*y1*cotb/rb + 
			(y1*y1)/(rb+y3b)*(pr+a/rb)) -
			/* term3 */
			(1-2*pr)*cotb/(rb+z3b)*(z1b*cosb - 
			a*(rb*sinb-y1)/(rb*cosb)) -
			/* term4 */
			a*y1*(y3b-a)*cotb/pow(rb,3.0) +
			/* term5 */
			(y3b-a)/(rb+y3b)*(2*pr + 1/rb*((1-2*pr)*y1*cotb+a) -
			(y1*y1)/(rb*(rb+y3b))*(2*pr+a/rb) - a*(y1*y1)/pow(rb,3.0)) +
			/* term6 */
			(y3b-a)*cotb/(rb+z3b)*
			(-cosb*sinb+a*y1*y3b/(pow(rb,3.0)*cosb) +
			(rb*sinb-y1)/rb*(2*(1-pr)*cosb - 
			(rb*cosb+y3b)/(rb+z3b) *
			(1+a/(rb*cosb))))
			);
			
		displ_ic_ft[1][0] += 0.0;
			
		displ_ic[1][1] += 2.0*(
			/* term1 */
			(1-2*pr)*y2/(rb+y3b)*(-(1-2*pr-a/rb)*cotb + 
			y1/(rb+y3b)*(pr+a/rb)) -
			/* term2 */
			(1-2*pr)*y2*cotb/(rb+z3b)*(1+a/(rb*cosb)) -
			/* term3 */
			a*y2*(y3b-a)*cotb/pow(rb,3.0) +
			/* term4 */
			y2*(y3b-a)/(rb*(rb+y3b))*((1-2*pr)*cotb - 2*pr*y1/(rb+y3b) -
			a*y1/rb*(1/rb+1/(rb+y3b))) +
			/* term5 */
			y2*(y3b-a)*cotb/(rb*(rb+z3b))*(-2*(1-pr)*cosb +
			(rb*cosb+y3b)/(rb+z3b)*(1+a/(rb*cosb)) +
			a*y3b/((rb*rb)*cosb))
			);
			
		displ_ic_ft[1][1] += 2* 2*(1-pr)*(1-2*pr)*fb*(cotb*cotb);
			
		displ_ic[1][2] += 2.0*(
			/* term1 */
			-2*(1-pr)*(1-2*pr)*cotb * 
			(log(rb+y3b)-cosb*log(rb+z3b)) -
			/* term2 */
			2*(1-pr)*y1/(rb+y3b)*(2*pr+a/rb) + 
			/* term3 */
			2*(1-pr)*z1b/(rb+z3b)*(cosb+a/rb) +
			/* term4 */
			(y3b-a)/rb*((1-2*pr)*cotb-2*pr*y1/(rb+y3b)-a*y1/(rb*rb)) -
			/* term5 */
			(y3b-a)/(rb+z3b)*(cosb*sinb +
			(rb*cosb+y3b)*cotb/rb*(2*(1-pr)*cosb - 
			(rb*cosb+y3b)/(rb+z3b)) +
			a/rb*(sinb - y3b*z1b/(rb*rb) - 
			z1b*(rb*cosb+y3b)/(rb*(rb+z3b))))
			);
			
		displ_ic_ft[1][2] += 0.0;

		term1 = 
                        (1-2*pr)*(y2/(rb+y3b)*(1+a/rb) -
                        y2*cosb/(rb+z3b)*(cosb+a/rb));

		term2 = y2*(y3b-a)/rb*(a/(rb*rb) + 1/(rb+y3b));

		term3 = y2*(y3b-a)*cosb/(rb*(rb+z3b))*
                        ((rb*cosb+y3b)/(rb+z3b)*(cosb+a/rb) +
                        a*y3b/(rb*rb)); 
			
		displ_ic[2][0] += 2.0*(
			/* term1 */
			(1-2*pr)*(y2/(rb+y3b)*(1+a/rb) - 
			y2*cosb/(rb+z3b)*(cosb+a/rb)) -
			/* term2 */
			y2*(y3b-a)/rb*(a/(rb*rb) + 1/(rb+y3b)) +
			/* term3 */
			y2*(y3b-a)*cosb/(rb*(rb+z3b))*
			((rb*cosb+y3b)/(rb+z3b)*(cosb+a/rb) +
			a*y3b/(rb*rb))
			);
			
		displ_ic_ft[2][0] += 0.0;
			
		displ_ic[2][1] += 2.0*(
			/* term1 */
			(1-2*pr)*(-sinb*log(rb+z3b) - y1/(rb+y3b)*(1+a/rb) +
			z1b/(rb+z3b)*(cosb+a/rb)) +
			/* term2 */
			y1*(y3b-a)/rb*(a/(rb*rb) + 1/(rb+y3b)) -
			/* term3 */
			(y3b-a)/(rb+z3b)*(sinb*(cosb-a/rb) + 
			z1b/rb*(1+a*y3b/(rb*rb)) -
			1/(rb*(rb+z3b))*((y2*y2)*cosb*sinb - 
			a*z1b/rb*(rb*cosb+y3b)))
			);
			
		displ_ic_ft[2][1] += 0.0;
			
		displ_ic[2][2] += 2.0*(
			/* term1 */
			2*(1-pr) * (y2*sinb/(rb+z3b)*(cosb + a/rb)) +
			/* term2 */
			y2*(y3b-a)*sinb/(rb*(rb+z3b))*
			(1 + (rb*cosb+y3b)/(rb+z3b) *
			(cosb+a/rb) + a*y3b/(rb*rb))
			);
			
		displ_ic_ft[2][2] += 2* 2*(1-pr) * (fb);

	}

	/* For debugging only -- Delete me */
	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			dic2[i][j] = 2.0*displ_ic[i][j];
		}
	}


	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {

			/* For debugging only
			----------------------
			displ_ic_ft[i][j] /= 2.*2.*(1-pr);
			trash1[i][j] = displ_ic_ft[i][j];
			displ_ic_ft[i][j] = reduce_angle(displ_ic_ft[i][j]);
			trash2[i][j] = displ_ic_ft[i][j];
			displ_ic_ft[i][j] *= 2.*2.*(1-pr);
			*/

			/* Add _ft terms to inf coeff eqns
			----------------------------------*/
			displ_ic[i][j] += displ_ic_ft[i][j];

			/* Divide through by consts on l.h.s.
			-------------------------------------*/
			displ_ic[i][j] /= 8.0*PI*(1-pr);
		}
	}

}

/***************** Function: comninou_strain_ics **********************
* Calculates the influence coefficients for strain for an angular dislocation.
*
* In:	y					- coords of point at which to calc inf coeffs
*		a					- depth to angular dislocation vertex
*		beta				- angle subtended by dislocation
*		pr					- poisson's ratio
*		half_space			- half/whole-space flag
* Out:	strain_ic[i][j][k]	- jk component of strain at y due to unit ith
*							  Burgers vector component
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     comninou_strain_ics(double y[3], double a, double beta, double pr,int	half_space, double strain_ic[3][3][3])
#else
void comninou_strain_ics(y, a, beta, pr, half_space, strain_ic)
double y[3];
double a;
double beta;
double pr;
int	 half_space;
double strain_ic[3][3][3];
#endif
{

	double	dummy[3][3][3];

	double	dg_ic[3][3][3];
	int		i, j, k;

	double	y1	= y[0];
	double	y2	= y[1];
	double	y3	= y[2];
	double	y3b	= y3 + 2.*a;

	double	sinb = sin(beta);
	double	cosb = cos(beta);
	double	cotb = 1./tan(beta);

	double	z1	= y1*cosb - y3*sinb;
	double	z3	= y1*sinb + y3*cosb;
	double	r	= sqrt(y1*y1 + y2*y2 + y3*y3);

	double	z1b	= y1*cosb + y3b*sinb;
	double	z3b = -y1*sinb + y3b*cosb;
	double	rb	= sqrt(y1*y1 + y2*y2 + y3b*y3b);

	double	rMy3P2		= (r-y3)*(r-y3);
	double	rP2			= r*r;
	double	y1P2		= y1*y1;
	double	y2P2		= y2*y2;
	double	y2P3		= y2P2*y2;
	double	y2P4		= y2P2*y2P2;
	double	rMz3P2		= (r-z3)*(r-z3);
	double	cosbP2		= cosb*cosb;
	double	sinbP2		= sinb*sinb;
	double	cotbP2		= cotb*cotb;
	double	rbPy3bP2	= (rb+y3b)*(rb+y3b);
	double	rbP2		= rb*rb;
	double	rbPz3bP2	= (rb+z3b)*(rb+z3b);
	double	rM1			= 1./r;
	double	rbM1		= 1./rb;
	double	rbP3		= rbP2*rb;
	double	rbP4		= rbP2*rbP2;
	double	rbPy3bM1	= 1./(rb+y3b);
	double	rbPz3b		= rb+z3b;
	double	rbPy3b		= rb+y3b;
	double	rMz3		= r-z3;
	double	rTsinbMy1	= r*sinb - y1;
	double	rTcosbMy3	= r*cosb - y3;
	double	rbTcosbPy3b	= rb*cosb + y3b;
	double	rbTsinbMy1	= rb*sinb - y1;
	double	rMy3		= r-y3;
	double	y3bMa		= y3b-a;
	double	pr2			= 2.*(1.-pr);
	double	pr3			= (1.-2.*pr);
	double	aQrb		= a/rb;
	double	z1P2		= z1*z1;
	double	z1bP2		= z1b*z1b;

	double	Dy3bDy3		= 1.0;
	double	Dz1Dy1		= cosb;
	double	Dz1Dy3		= -sinb;
	double	Dz3Dy1		= sinb;
	double	Dz3Dy3		= cosb;
	double	DrDy1		= y1/r;
	double	DrDy2		= y2/r;
	double	DrDy3		= y3/r;
	double	Dz1bDy1		= cosb;
	double	Dz1bDy3b	= sinb;
	double	Dz3bDy1		= -sinb;
	double	Dz3bDy3b	= cosb;
	double	DrbDy1		= y1/rb;
	double	DrbDy2		= y2/rb;
	double	DrbDy3b		= y3b/rb;


	double	Df_denom = 
		rP2*sinbP2*y2P2 + cosbP2*y2P4 + 2*cosb*y1*y2P2*z1 + y1P2*z1P2;
		
	double	Dfb_denom =
		rbP2*sinbP2*y2P2 + cosbP2*y2P4 + 2*cosb*y1*y2P2*z1b + y1P2*z1bP2;
		
	double	DfDy1 = 
		y2/(y1P2 + y2P2) - (r*sinb*y2*z1)/ Df_denom;
		
	double	DfDy2 = 
		-(y1/(y1P2 + y2P2)) + z1/(y2P2 + z1P2) + 
		(-(cosb*r*sinb*y2P2) + r*sinb*y1*z1)/ Df_denom;
		
	double	DfDz1 = 
		-(y2/(y2P2 + z1P2)) - (r*sinb*y1*y2)/ Df_denom;
		
	double	DfDr = 
		(sinb*y2*(cosb*y2P2 + y1*z1))/ Df_denom;
		
	double	DfbDy1 = 
		y2/(y1P2 + y2P2) - (rb*sinb*y2*z1b)/ Dfb_denom;
		
	double	DfbDy2 = 
		-(y1/(y1P2 + y2P2)) + z1b/(y2P2 + z1bP2) + 
		(-(cosb*rb*sinb*y2P2) + rb*sinb*y1*z1b)/ Dfb_denom;
		
	double	DfbDz1b = 
		-(y2/(y2P2 + z1bP2)) - (rb*sinb*y1*y2)/ Dfb_denom;
		
	double	DfbDrb = 
		(sinb*y2*(cosb*y2P2 + y1*z1b))/ Dfb_denom;
	
	
	/* Dislocation-Indiced Displacement Gradients */
		
	dg_ic[0][0][0] =
		pr2*(DrDy1* DfDr + Dz1Dy1* DfDz1 + DfDy1) - 
		y2/(rMy3*r) - (cosb*y2*(-1. + sinb*DrDy1))/ (rMz3*r) + 
		(cosb*y2*rTsinbMy1* (-Dz3Dy1 + DrDy1))/ (rMz3P2*r) + 
		(y1*y2*DrDy1)/ (rMy3*rP2) + (cosb*y2*rTsinbMy1*DrDy1)/
		(rMz3*rP2) + (y1*y2*DrDy1)/ (rMy3P2*r);
		
	dg_ic[0][0][1] =
		pr2*(DrDy2* DfDr + DfDy2) - y1/(rMy3*r) - 
		(cosb*rTsinbMy1)/ (rMz3*r) + (y1*y2*DrDy2)/ (rMy3*rP2) + 
		(cosb*y2*rTsinbMy1*DrDy2)/ (rMz3*rP2) + (y1*y2*DrDy2)/ (rMy3P2*r) + 
		(cosb*y2*rTsinbMy1*DrDy2)/ (rMz3P2*r) - 
		(cosb*sinb*y2*DrDy2)/ (rMz3*r);
		
	dg_ic[0][0][2] =
		pr2*(DrDy3* DfDr + Dz1Dy3* DfDz1) + 
		(y1*y2*(-1. + DrDy3))/ (rMy3P2*r) + 
		(cosb*y2*rTsinbMy1* (-Dz3Dy3 + DrDy3))/ (rMz3P2*r) + 
		(y1*y2*DrDy3)/ (rMy3*rP2) + 
		(cosb*y2*rTsinbMy1*DrDy3)/ (rMz3*rP2) - 
		(cosb*sinb*y2*DrDy3)/ (rMz3*r);
	 
	dg_ic[0][1][0] =
		pr3*(-((cosb*(-Dz3Dy1 + DrDy1))/ rMz3) + 
		DrDy1/rMy3) - y2P2*((cosb*(-Dz3Dy1 + DrDy1))/ (rMz3P2*r) - 
		DrDy1/ (rMy3*rP2) + (cosb*DrDy1)/ (rMz3*rP2) - 
		DrDy1/ (rMy3P2*r)); 
		
	dg_ic[0][1][1] =
		-2.*y2*(1/(rMy3*r) - cosb/(rMz3*r)) + pr3*(DrDy2/rMy3 - 
		(cosb*DrDy2)/ rMz3) - y2P2*(-(DrDy2/ (rMy3*rP2)) + 
		(cosb*DrDy2)/ (rMz3*rP2) - DrDy2/ (rMy3P2*r) + 
		(cosb*DrDy2)/ (rMz3P2*r));
		
	dg_ic[0][1][2] =
		pr3*((-1. + DrDy3)/rMy3 - (cosb*(-Dz3Dy3 + DrDy3))/ rMz3) - 
		y2P2*(-((-1. + DrDy3)/ (rMy3P2*r)) + (cosb*(-Dz3Dy3 + 
		DrDy3))/ (rMz3P2*r) - DrDy3/ (rMy3*rP2) + (cosb*DrDy3)/ (rMz3*rP2));
		 
	dg_ic[0][2][0] =
		y2*((cosb*rTcosbMy3* (-Dz3Dy1 + DrDy1))/ (rMz3P2*r) - 
		DrDy1/rP2 + (cosb*rTcosbMy3*DrDy1)/ (rMz3*rP2) - 
		(cosbP2*DrDy1)/ (rMz3*r));
		
	dg_ic[0][2][1] =
		y2*(-(DrDy2/rP2) + (cosb*rTcosbMy3*DrDy2) /(rMz3*rP2) + 
		(cosb*rTcosbMy3*DrDy2) /(rMz3P2*r) - (cosbP2*DrDy2)/ (rMz3*r)) + 
		rM1 - (cosb*rTcosbMy3)/ (rMz3*r); 
		
	dg_ic[0][2][2] =
		y2*(-((cosb*(-1. + cosb*DrDy3))/ (rMz3*r)) + 
		(cosb*rTcosbMy3* (-Dz3Dy3 + DrDy3))/ (rMz3P2*r) - 
		DrDy3/rP2 + (cosb*rTcosbMy3*DrDy3)/ (rMz3*rP2));
		
	dg_ic[1][0][0] =
		-(pr3*(-((cosb*(-Dz3Dy1 + DrDy1))/ rMz3) + DrDy1/rMy3)) + 
		(2.*y1)/(rMy3*r) + ((-1. + sinb*DrDy1)*z1)/ (rMz3*r) - 
		(rTsinbMy1* (-Dz3Dy1 + DrDy1)* z1)/(rMz3P2*r) + 
		(rTsinbMy1*Dz1Dy1)/ (rMz3*r) - (y1P2*DrDy1)/ (rMy3*rP2) - 
		(y1P2*DrDy1)/ (rMy3P2*r) - (rTsinbMy1*z1*DrDy1)/ (rMz3*rP2);
		
	dg_ic[1][0][1] =
		-(pr3*(DrDy2/rMy3 - (cosb*DrDy2)/ rMz3)) - (y1P2*DrDy2)/ (rMy3*rP2) - 
		(y1P2*DrDy2)/ (rMy3P2*r) - (rTsinbMy1*z1*DrDy2)/ (rMz3*rP2) - 
		(rTsinbMy1*z1*DrDy2)/ (rMz3P2*r) + (sinb*z1*DrDy2)/ (rMz3*r);
		
	dg_ic[1][0][2] =
		-(pr3*((-1. + DrDy3)/rMy3 - (cosb*(-Dz3Dy3 + DrDy3))/ rMz3)) - 
		(y1P2*(-1. + DrDy3))/ (rMy3P2*r) - 
		(rTsinbMy1* (-Dz3Dy3 + DrDy3)* z1)/(rMz3P2*r) + 
		(rTsinbMy1*Dz1Dy3)/ (rMz3*r) - (y1P2*DrDy3)/ (rMy3*rP2) - 
		(rTsinbMy1*z1*DrDy3)/ (rMz3*rP2) + (sinb*z1*DrDy3)/ (rMz3*r);
		 
	dg_ic[1][1][0] =
		pr2*(DrDy1* DfDr + Dz1Dy1* DfDz1 + DfDy1) + y2/(rMy3*r) + 
		(y2*(-Dz3Dy1 + DrDy1)* z1)/(rMz3P2*r) - (y2*Dz1Dy1)/ (rMz3*r) - 
		(y1*y2*DrDy1)/ (rMy3*rP2) - (y1*y2*DrDy1)/ (rMy3P2*r) + 
		(y2*z1*DrDy1)/ (rMz3*rP2); 
		
	dg_ic[1][1][1] =
		pr2*(DrDy2* DfDr + DfDy2) + y1/(rMy3*r) - z1/(rMz3*r) - 
		(y1*y2*DrDy2)/ (rMy3*rP2) - (y1*y2*DrDy2)/ (rMy3P2*r) + 
		(y2*z1*DrDy2)/ (rMz3*rP2) + (y2*z1*DrDy2)/ (rMz3P2*r);
		
	dg_ic[1][1][2] =
		pr2*(DrDy3* DfDr + Dz1Dy3* DfDz1) - 
		(y1*y2*(-1. + DrDy3))/ (rMy3P2*r) + 
		(y2*(-Dz3Dy3 + DrDy3)* z1)/(rMz3P2*r) - (y2*Dz1Dy3)/ (rMz3*r) - 
		(y1*y2*DrDy3)/ (rMy3*rP2) + (y2*z1*DrDy3)/ (rMz3*rP2);
		 
	dg_ic[1][2][0] =
		-((sinb*pr3*(-Dz3Dy1 + DrDy1))/ rMz3) - rM1 - 
		(rTcosbMy3* (-Dz3Dy1 + DrDy1)* z1)/(rMz3P2*r) + 
		(rTcosbMy3*Dz1Dy1)/ (rMz3*r) + (y1*DrDy1)/rP2 - 
		(rTcosbMy3*z1*DrDy1)/ (rMz3*rP2) + (cosb*z1*DrDy1)/ (rMz3*r); 

	dg_ic[1][2][1] =
		-((sinb*pr3*DrDy2)/ rMz3) + (y1*DrDy2)/rP2 - 
		(rTcosbMy3*z1*DrDy2)/ (rMz3*rP2) - 
		(rTcosbMy3*z1*DrDy2)/ (rMz3P2*r) + 
		(cosb*z1*DrDy2)/ (rMz3*r); 
		
	dg_ic[1][2][2] =
		-((sinb*pr3*(-Dz3Dy3 + DrDy3))/ rMz3) + 
		((-1. + cosb*DrDy3)*z1)/ (rMz3*r) - 
		(rTcosbMy3* (-Dz3Dy3 + DrDy3)* z1)/(rMz3P2*r) + 
		(rTcosbMy3*Dz1Dy3)/ (rMz3*r) + (y1*DrDy3)/rP2 - 
		(rTcosbMy3*z1*DrDy3)/ (rMz3*rP2);
		
	dg_ic[2][0][0] =
		(sinb*y2*(-1. + sinb*DrDy1))/ (rMz3*r) - 
		(sinb*y2*rTsinbMy1* (-Dz3Dy1 + DrDy1))/ (rMz3P2*r) - 
		(sinb*y2*rTsinbMy1*DrDy1)/ (rMz3*rP2);
		
	dg_ic[2][0][1] =
		(sinb*rTsinbMy1)/(rMz3*r) - (sinb*y2*rTsinbMy1*DrDy2)/
		(rMz3*rP2) - (sinb*y2*rTsinbMy1*DrDy2)/ (rMz3P2*r) + 
		(sinbP2*y2*DrDy2)/ (rMz3*r);
		
	dg_ic[2][0][2] =
		-((sinb*y2*rTsinbMy1* (-Dz3Dy3 + DrDy3))/ (rMz3P2*r)) - 
		(sinb*y2*rTsinbMy1*DrDy3)/ (rMz3*rP2) + (sinbP2*y2*DrDy3)/ (rMz3*r);
		 
	dg_ic[2][1][0] =
		(sinb*pr3*(-Dz3Dy1 + DrDy1))/ rMz3 + (sinb*y2P2*(-Dz3Dy1 + 
		DrDy1))/ (rMz3P2*r) + (sinb*y2P2*DrDy1)/ (rMz3*rP2); 
		
	dg_ic[2][1][1] =
		(-2.*sinb*y2)/(rMz3*r) + (sinb*pr3*DrDy2)/ rMz3 + 
		(sinb*y2P2*DrDy2)/ (rMz3*rP2) + (sinb*y2P2*DrDy2)/ (rMz3P2*r);
		
	dg_ic[2][1][2] =
		(sinb*pr3*(-Dz3Dy3 + DrDy3))/ rMz3 + (sinb*y2P2*(-Dz3Dy3 + 
		DrDy3))/ (rMz3P2*r) + (sinb*y2P2*DrDy3)/ (rMz3*rP2);
	 
	dg_ic[2][2][0] =
		pr2*(DrDy1* DfDr + Dz1Dy1* DfDz1 + DfDy1) - 
		(sinb*y2*rTcosbMy3* (-Dz3Dy1 + DrDy1))/ (rMz3P2*r) - 
		(sinb*y2*rTcosbMy3*DrDy1)/ (rMz3*rP2) + 
		(cosb*sinb*y2*DrDy1)/ (rMz3*r);
		
	dg_ic[2][2][1] =
		pr2*(DrDy2* DfDr + DfDy2) + 
		(sinb*rTcosbMy3)/ (rMz3*r) - (sinb*y2*rTcosbMy3*DrDy2)/ (rMz3*rP2) - 
		(sinb*y2*rTcosbMy3*DrDy2)/ (rMz3P2*r) + (cosb*sinb*y2*DrDy2)/ (rMz3*r); 
		
	dg_ic[2][2][2] =
		pr2*(DrDy3* DfDr + Dz1Dy3* DfDz1) + 
		(sinb*y2*(-1. + cosb*DrDy3))/ (rMz3*r) - 
		(sinb*y2*rTcosbMy3* (-Dz3Dy3 + DrDy3))/ (rMz3P2*r) - 
		(sinb*y2*rTcosbMy3*DrDy3)/ (rMz3*rP2); 
 

	if (half_space) {

		/* Image-Induced Displacement Gradients */
			
		dg_ic[0][0][0] +=
			pr2*(DrbDy1* DfbDrb + Dz1bDy1* DfbDz1b + DfbDy1) - 
			y2/(rbPy3b*rb) - (cosb*y2*(-1. + sinb*DrbDy1))/
			(rbPz3b*rb) + (cosb*y2*rbTsinbMy1* (Dz3bDy1 + DrbDy1))/
			(rbPz3bP2* rb) + (y1*y2*DrbDy1)/ (rbPy3b*rbP2) + 
			(cosb*y2*rbTsinbMy1* DrbDy1)/ (rbPz3b* rbP2) + 
			(y1*y2*DrbDy1)/ (rbPy3bP2*rb);
			
		dg_ic[0][0][1] +=
			pr2*(DrbDy2* DfbDrb + DfbDy2) - 
			y1/(rbPy3b*rb) - (cosb*rbTsinbMy1)/ (rbPz3b*rb) + (y1*y2*DrbDy2)/
			(rbPy3b*rbP2) + (cosb*y2*rbTsinbMy1* DrbDy2)/ (rbPz3b* rbP2) + 
			(y1*y2*DrbDy2)/ (rbPy3bP2*rb) + (cosb*y2*rbTsinbMy1* DrbDy2)/
			(rbPz3bP2* rb) - (cosb*sinb*y2*DrbDy2)/ (rbPz3b*rb);
			
		dg_ic[0][0][2] +=
			pr2*(Dy3bDy3* DrbDy3b* DfbDrb + Dz1bDy3b*Dy3bDy3* DfbDz1b) + 
			(y1*y2*(Dy3bDy3 + Dy3bDy3* DrbDy3b))/ (rbPy3bP2*rb) + 
			(cosb*y2*rbTsinbMy1* (Dz3bDy3b*Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			(rbPz3bP2* rb) + (y1*y2*Dy3bDy3* DrbDy3b)/ (rbPy3b*rbP2) + 
			(cosb*y2*rbTsinbMy1* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) - 
			(cosb*sinb*y2*Dy3bDy3* DrbDy3b)/ (rbPz3b*rb);
			 
		dg_ic[0][1][0] +=
			pr3*(-((cosb*(Dz3bDy1 + DrbDy1))/ rbPz3b) + 
			DrbDy1/ rbPy3b) - y2P2*((cosb*(Dz3bDy1 + DrbDy1))/
			(rbPz3bP2* rb) - DrbDy1/ (rbPy3b*rbP2) + (cosb*DrbDy1)/
			(rbPz3b* rbP2) - DrbDy1/ (rbPy3bP2*rb)); 
			
		dg_ic[0][1][1] +=
			-2.*y2*(1/(rbPy3b*rb) - cosb/ (rbPz3b* rb)) + 
			pr3*(DrbDy2/ rbPy3b - (cosb*DrbDy2)/ rbPz3b) - 
			y2P2*(-(DrbDy2/ (rbPy3b* rbP2)) + (cosb*DrbDy2)/ (rbPz3b* rbP2) - 
			DrbDy2/ (rbPy3bP2*rb) + (cosb*DrbDy2)/ (rbPz3bP2* rb));
			
		dg_ic[0][1][2] +=
			pr3*((Dy3bDy3 + Dy3bDy3* DrbDy3b)/ rbPy3b - 
			(cosb*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3b) - 
			y2P2*(-((Dy3bDy3 + Dy3bDy3* DrbDy3b)/ (rbPy3bP2* rb)) + 
			(cosb*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ (rbPz3bP2* rb) - 
			(Dy3bDy3* DrbDy3b)/ (rbPy3b*rbP2) + (cosb*Dy3bDy3*
			DrbDy3b)/ (rbPz3b* rbP2));
			 
		dg_ic[0][2][0] +=
			y2*(-((cosb*rbTcosbPy3b* (Dz3bDy1 + DrbDy1))/ (rbPz3bP2* rb)) + 
			DrbDy1/rbP2 - (cosb*rbTcosbPy3b* DrbDy1)/ (rbPz3b* rbP2) + 
			(cosbP2*DrbDy1)/ (rbPz3b* rb));
			
		dg_ic[0][2][1] +=
			y2*(DrbDy2/rbP2 - (cosb*rbTcosbPy3b* DrbDy2)/ (rbPz3b* rbP2) + 
			(cosbP2*DrbDy2)/ (rbPz3b* rb) - 
			(cosb*rbTcosbPy3b* DrbDy2)/ (rbPz3bP2* rb)) - rbM1 + 
			(cosb*rbTcosbPy3b)/ (rbPz3b*rb); 
			
		dg_ic[0][2][2] +=
			y2*((cosb*(Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/ (rbPz3b* rb) - 
			(cosb*rbTcosbPy3b* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			(rbPz3bP2* rb) + (Dy3bDy3* DrbDy3b)/ rbP2 - 
			(cosb*rbTcosbPy3b* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2));
			
		dg_ic[1][0][0] +=
			-(pr3*(-((cosb*(Dz3bDy1 + DrbDy1))/ rbPz3b) + DrbDy1/ rbPy3b)) + 
			(2.*y1)/(rbPy3b*rb) + ((-1. + sinb*DrbDy1)* z1b)/
			(rbPz3b*rb) - (rbTsinbMy1* (Dz3bDy1 + DrbDy1)* z1b)/
			(rbPz3bP2* rb) + (rbTsinbMy1*Dz1bDy1) /(rbPz3b* rb) - 
			(y1P2*DrbDy1)/ (rbPy3b*rbP2) - (y1P2*DrbDy1)/ (rbPy3bP2*rb) - 
			(rbTsinbMy1*z1b* DrbDy1)/ (rbPz3b* rbP2);
			
		dg_ic[1][0][1] +=
			-(pr3*(DrbDy2/ rbPy3b - (cosb*DrbDy2)/ rbPz3b)) - 
			(y1P2*DrbDy2)/ (rbPy3b*rbP2) - (y1P2*DrbDy2)/ (rbPy3bP2*rb) - 
			(rbTsinbMy1*z1b* DrbDy2)/ (rbPz3b* rbP2) - 
			(rbTsinbMy1*z1b* DrbDy2)/ (rbPz3bP2* rb) + 
			(sinb*z1b*DrbDy2)/ (rbPz3b*rb);
			
		dg_ic[1][0][2] +=
			-(pr3*((Dy3bDy3 + Dy3bDy3* DrbDy3b)/ rbPy3b - 
			(cosb*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3b)) - 
			(y1P2*(Dy3bDy3 + Dy3bDy3* DrbDy3b))/ (rbPy3bP2*rb) - 
			(rbTsinbMy1* (Dz3bDy3b*Dy3bDy3 + Dy3bDy3* DrbDy3b)*
			z1b)/ (rbPz3bP2* rb) + (rbTsinbMy1* Dz1bDy3b*Dy3bDy3)/
			(rbPz3b*rb) - (y1P2*Dy3bDy3* DrbDy3b)/ (rbPy3b*rbP2) - 
			(rbTsinbMy1*z1b* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) + 
			(sinb*z1b*Dy3bDy3* DrbDy3b)/ (rbPz3b*rb);
			 
		dg_ic[1][1][0] +=
			pr2*(DrbDy1* DfbDrb + Dz1bDy1* DfbDz1b + DfbDy1) + 
			y2/(rbPy3b*rb) + (y2*(Dz3bDy1 + DrbDy1)* z1b)/ (rbPz3bP2* rb) - 
			(y2*Dz1bDy1)/ (rbPz3b*rb) - (y1*y2*DrbDy1)/ (rbPy3b*rbP2) - 
			(y1*y2*DrbDy1)/ (rbPy3bP2*rb) + (y2*z1b*DrbDy1)/ (rbPz3b* rbP2); 
			
		dg_ic[1][1][1] +=
			pr2*(DrbDy2* DfbDrb + DfbDy2) + y1/(rbPy3b*rb) - 
			z1b/ (rbPz3b*rb) - (y1*y2*DrbDy2)/ (rbPy3b*rbP2) - 
			(y1*y2*DrbDy2)/ (rbPy3bP2*rb) + (y2*z1b*DrbDy2)/ (rbPz3b* rbP2) + 
			(y2*z1b*DrbDy2)/ (rbPz3bP2* rb);
			
		dg_ic[1][1][2] +=
			pr2*(Dy3bDy3* DrbDy3b* DfbDrb + Dz1bDy3b*Dy3bDy3* DfbDz1b) - 
			(y1*y2*(Dy3bDy3 + Dy3bDy3* DrbDy3b))/ (rbPy3bP2*rb) + 
			(y2*(Dz3bDy3b*Dy3bDy3 + Dy3bDy3* DrbDy3b)* z1b)/ (rbPz3bP2* rb) - 
			(y2*Dz1bDy3b*Dy3bDy3)/ (rbPz3b*rb) - (y1*y2*Dy3bDy3* DrbDy3b)/
			(rbPy3b*rbP2) + (y2*z1b*Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2);
			 
		dg_ic[1][2][0] +=
			(sinb*pr3*(Dz3bDy1 + DrbDy1))/ rbPz3b + rbM1 + 
			(rbTcosbPy3b* (Dz3bDy1 + DrbDy1)* z1b)/ (rbPz3bP2* rb) - 
			(rbTcosbPy3b* Dz1bDy1)/
			(rbPz3b*rb) - (y1*DrbDy1)/rbP2 + (rbTcosbPy3b*z1b* DrbDy1)/
			(rbPz3b* rbP2) - (cosb*z1b*DrbDy1)/ (rbPz3b*rb);
			
		dg_ic[1][2][1] +=
			(sinb*pr3*DrbDy2)/ rbPz3b - (y1*DrbDy2)/rbP2 + 
			(rbTcosbPy3b*z1b* DrbDy2)/ (rbPz3b* rbP2) - 
			(cosb*z1b*DrbDy2)/ (rbPz3b*rb) + (rbTcosbPy3b*z1b* DrbDy2)/
			(rbPz3bP2* rb); 
			
		dg_ic[1][2][2] +=
			(sinb*pr3*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPz3b - ((Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b)*
			z1b)/ (rbPz3b*rb) + (rbTcosbPy3b* (Dz3bDy3b*Dy3bDy3 + Dy3bDy3*
			DrbDy3b)* z1b)/ (rbPz3bP2* rb) - (rbTcosbPy3b*
			Dz1bDy3b*Dy3bDy3)/ (rbPz3b*rb) - (y1*Dy3bDy3* DrbDy3b)/ rbP2 + 
			(rbTcosbPy3b*z1b* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2);
			
		dg_ic[2][0][0] +=
			(sinb*y2*(-1. + sinb*DrbDy1))/
			(rbPz3b*rb) - (sinb*y2*rbTsinbMy1* (Dz3bDy1 + DrbDy1))/
			(rbPz3bP2* rb) - (sinb*y2*rbTsinbMy1* DrbDy1)/ (rbPz3b* rbP2);
			
		dg_ic[2][0][1] +=
			(sinb*rbTsinbMy1)/ (rbPz3b*rb) - (sinb*y2*rbTsinbMy1* DrbDy2)/
			(rbPz3b* rbP2) - (sinb*y2*rbTsinbMy1* DrbDy2)/ (rbPz3bP2* rb) + 
			(sinbP2*y2*DrbDy2)/ (rbPz3b*rb);
			
		dg_ic[2][0][2] +=
			-((sinb*y2*rbTsinbMy1* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			(rbPz3bP2* rb)) - (sinb*y2*rbTsinbMy1* Dy3bDy3* DrbDy3b)/
			(rbPz3b* rbP2) + (sinbP2*y2*Dy3bDy3* DrbDy3b)/ (rbPz3b*rb);
			 
		dg_ic[2][1][0] +=
			(sinb*pr3*(Dz3bDy1 + DrbDy1))/ rbPz3b + 
			(sinb*y2P2*(Dz3bDy1 + DrbDy1))/ (rbPz3bP2* rb) + 
			(sinb*y2P2*DrbDy1)/ (rbPz3b* rbP2); 
			
		dg_ic[2][1][1] +=
			(-2.*sinb*y2)/ (rbPz3b*rb) + (sinb*pr3*DrbDy2)/ rbPz3b + 
			(sinb*y2P2*DrbDy2)/ (rbPz3b* rbP2) + (sinb*y2P2*DrbDy2)/
			(rbPz3bP2* rb);
			
		dg_ic[2][1][2] +=
			(sinb*pr3*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPz3b + (sinb*y2P2*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			(rbPz3bP2* rb) + (sinb*y2P2*Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2);
			 
		dg_ic[2][2][0] +=
			-pr2*(DrbDy1* DfbDrb + Dz1bDy1* DfbDz1b + DfbDy1) + 
			(sinb*y2*rbTcosbPy3b* (Dz3bDy1 + DrbDy1))/ (rbPz3bP2* rb) + 
			(sinb*y2*rbTcosbPy3b* DrbDy1)/ (rbPz3b* rbP2) - 
			(cosb*sinb*y2*DrbDy1)/ (rbPz3b*rb);
			
		dg_ic[2][2][1] +=
			-pr2*(DrbDy2* DfbDrb + DfbDy2) - (sinb*rbTcosbPy3b)/
			(rbPz3b*rb) + (sinb*y2*rbTcosbPy3b* DrbDy2)/ (rbPz3b*
			rbP2) - (cosb*sinb*y2*DrbDy2)/ (rbPz3b*rb) + (sinb*y2*rbTcosbPy3b*
			DrbDy2)/ (rbPz3bP2* rb); 
			
		dg_ic[2][2][2] +=
			-pr2*(Dy3bDy3* DrbDy3b* DfbDrb + Dz1bDy3b*Dy3bDy3* DfbDz1b) - 
			(sinb*y2*(Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/
			(rbPz3b*rb) + (sinb*y2*rbTcosbPy3b* (Dz3bDy3b*Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ (rbPz3bP2* rb) + (sinb*y2*rbTcosbPy3b*
			Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2);
			
		/* Corrective Displacement Gradients */
			
		dg_ic[0][0][0] +=
			2.*(-((cosb*cotb*y2*pr3*(cosb + aQrb)* (Dz3bDy1 + DrbDy1))/
			rbPz3bP2) - 2.*cotbP2*pr3*(1. - pr)* (DrbDy1* DfbDrb + 
			Dz1bDy1* DfbDz1b + DfbDy1) + (y2*pr3*(-((pr + aQrb)/
			rbPy3b) + (y1*(pr + aQrb)* DrbDy1)/ rbPy3bP2 + (a*cotb*DrbDy1)/
			rbP2 + (a*y1*DrbDy1)/ (rbPy3b* rbP2)))/ rbPy3b - (y2*y3bMa*
			((cosb*(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b))/ rbPz3b - (a*cosb*cotb*y3b)/rbP2)* (Dz3bDy1 + 
			DrbDy1))/ (rbPz3bP2* rb) + (y2*y3bMa*
			(-((cosb*(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b)* (Dz3bDy1 + DrbDy1))/ rbPz3bP2) +
			(cosb* (2.*cosb* (-1. + sinb*DrbDy1)* (1. - pr) + cosb*cotb*
			(cosb*pr3 - aQrb)* DrbDy1 + (a*cotb*rbTcosbPy3b* DrbDy1)/
			rbP2))/ rbPz3b + (2.*a*cosb*cotb*y3b* DrbDy1)/ rbP3))/
			(rbPz3b* rb) + (y2*y3bMa* ((2.*pr + aQrb)/ rbPy3b + a/rbP2 - 
			(y1*(2.*pr + aQrb)* DrbDy1)/ rbPy3bP2 - (2.*a*y1*DrbDy1)/
			rbP3 - (a*y1*DrbDy1)/ (rbPy3b* rbP2)))/ (rbPy3b*rb) - 
			(y2*pr3*(-((y1*(pr + aQrb))/ rbPy3b) + cotb*(1. - 2.*pr - aQrb))*
			DrbDy1)/ rbPy3bP2 - (3*a*cotb*y2*y3bMa* DrbDy1)/ rbP4 - 
			(a*cosb*cotb*y2*pr3*DrbDy1) /(rbPz3b* rbP2) - (y2*y3bMa*
			((cosb*(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b))/ rbPz3b - (a*cosb*cotb*y3b)/rbP2)* DrbDy1)/
			(rbPz3b* rbP2) - (y2*y3bMa* (-(cotb*pr3) + (y1*(2.*pr + aQrb))/
			rbPy3b + (a*y1)/rbP2)* DrbDy1)/ (rbPy3b*rbP2) - (y2*y3bMa*
			(-(cotb*pr3) + (y1*(2.*pr + aQrb))/ rbPy3b + (a*y1)/rbP2)*
			DrbDy1)/ (rbPy3bP2*rb));
			
		dg_ic[0][0][1] +=
			2.*((pr3*(-((y1*(pr + aQrb))/ rbPy3b) + 
			cotb*(1. - 2.*pr - aQrb)))/ rbPy3b + 
			(cosb*cotb*pr3*(cosb + aQrb))/ rbPz3b - 
			2.*cotbP2*pr3*(1. - pr)* (DrbDy2* DfbDrb + DfbDy2) + 
			(y2*pr3*((y1*(pr + aQrb)* DrbDy2)/ rbPy3bP2 + (a*cotb*DrbDy2)/
			rbP2 + (a*y1*DrbDy2)/ (rbPy3b* rbP2)))/ rbPy3b + 
			(a*cotb*y3bMa)/rbP3 + (y3bMa*((cosb*
			(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b))/ rbPz3b - (a*cosb*cotb*y3b)/rbP2))/ (rbPz3b* rb) + 
			(y3bMa*(-(cotb*pr3) + (y1*(2.*pr + aQrb))/ rbPy3b + 
			(a*y1)/rbP2))/ (rbPy3b*rb) + (y2*y3bMa*
			((cosb*(2.*cosb*sinb*(1. - pr)* DrbDy2 + cosb*cotb*
			(cosb*pr3 - aQrb)* DrbDy2 + (a*cotb*rbTcosbPy3b* DrbDy2)/
			rbP2))/ rbPz3b - (cosb*
			(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b)* DrbDy2)/ rbPz3bP2 + (2.*a*cosb*cotb*y3b*
			DrbDy2)/ rbP3))/ (rbPz3b* rb) + (y2*y3bMa* (-((y1*(2.*pr + aQrb)*
			DrbDy2)/ rbPy3bP2) - (2.*a*y1*DrbDy2)/ rbP3 - (a*y1*DrbDy2)/
			(rbPy3b* rbP2)))/ (rbPy3b*rb) - (y2*pr3*(-((y1*(pr + aQrb))/
			rbPy3b) + cotb*(1. - 2.*pr - aQrb))* DrbDy2)/
			rbPy3bP2 - (cosb*cotb*y2*pr3*(cosb + aQrb)* DrbDy2)/
			rbPz3bP2 - (3*a*cotb*y2*y3bMa* DrbDy2)/ rbP4 - 
			(a*cosb*cotb*y2*pr3*DrbDy2) /(rbPz3b* rbP2) - (y2*y3bMa*
			((cosb*(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b))/ rbPz3b - (a*cosb*cotb*y3b)/rbP2)* DrbDy2)/
			(rbPz3b* rbP2) - (y2*y3bMa* (-(cotb*pr3) + (y1*(2.*pr + aQrb))/
			rbPy3b + (a*y1)/rbP2)* DrbDy2)/ (rbPy3b*rbP2) - (y2*y3bMa*
			((cosb*(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b))/ rbPz3b - (a*cosb*cotb*y3b)/rbP2)* DrbDy2)/
			(rbPz3bP2* rb) - (y2*y3bMa* (-(cotb*pr3) + (y1*(2.*pr + aQrb))/
			rbPy3b + (a*y1)/rbP2)* DrbDy2)/ (rbPy3bP2*rb));
			
		dg_ic[0][0][2] +=
			2.*(-((y2*pr3*(-((y1*(pr + aQrb))/ rbPy3b) + 
			cotb*(1. - 2.*pr - aQrb))* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2) - (cosb*cotb*y2*pr3*(cosb + aQrb)* (Dz3bDy3b*
			Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3bP2 - 
			2.*cotbP2*pr3*(1. - pr)* (Dy3bDy3* DrbDy3b* DfbDrb + 
			Dz1bDy3b* Dy3bDy3* DfbDz1b) + (y2*pr3*((y1*(pr + aQrb)*
			(Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPy3bP2 + (a*cotb*Dy3bDy3*
			DrbDy3b)/ rbP2 + (a*y1*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP2)))/
			rbPy3b - (y2*y3bMa* ((cosb*(pr2*cosb*rbTsinbMy1 + 
			cotb*(cosb*pr3 - aQrb)* rbTcosbPy3b))/ rbPz3b - 
			(a*cosb*cotb*y3b)/rbP2)* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ (rbPz3bP2* rb) - (y2*y3bMa* (Dy3bDy3 + Dy3bDy3*
			DrbDy3b)* (-(cotb*pr3) + (y1*(2.*pr + aQrb))/ rbPy3b + 
			(a*y1)/rbP2))/ (rbPy3bP2*rb) + (y2*y3bMa* (-((y1*(2.*pr + aQrb)*
			(Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPy3bP2) - (2.*a*y1*Dy3bDy3*
			DrbDy3b)/ rbP3 - (a*y1*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP2)))/
			(rbPy3b*rb) + (y2*y3bMa*
			(-((cosb*(pr2*cosb*rbTsinbMy1 + cotb*(cosb*pr3 - aQrb)*
			rbTcosbPy3b)* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPz3bP2) + (cosb* (cotb*(cosb*pr3 - aQrb)* (Dy3bDy3 + 
			cosb*Dy3bDy3* DrbDy3b) + 2.*cosb*sinb*(1. - pr)*Dy3bDy3*
			DrbDy3b + (a*cotb*rbTcosbPy3b* Dy3bDy3* DrbDy3b)/ rbP2))/
			rbPz3b - (a*cosb*cotb*Dy3bDy3)/ rbP2 + (2.*a*cosb*cotb*y3b*
			Dy3bDy3* DrbDy3b)/ rbP3))/ (rbPz3b* rb) + 
			(a*cotb*y2*Dy3bDy3)/rbP3 + (y2*((cosb*(pr2*cosb*rbTsinbMy1 + 
			cotb*(cosb*pr3 - aQrb)* rbTcosbPy3b))/ rbPz3b - 
			(a*cosb*cotb*y3b)/rbP2)* Dy3bDy3)/ (rbPz3b* rb) + 
			(y2*(-(cotb*pr3) + (y1*(2.*pr + aQrb))/ rbPy3b + (a*y1)/rbP2)*
			Dy3bDy3)/ (rbPy3b*rb) - (3*a*cotb*y2*y3bMa*Dy3bDy3*
			DrbDy3b)/ rbP4 - (a*cosb*cotb*y2*pr3*Dy3bDy3* DrbDy3b)/
			(rbPz3b* rbP2) - (y2*y3bMa* ((cosb*(pr2*cosb*rbTsinbMy1 + 
			cotb*(cosb*pr3 - aQrb)* rbTcosbPy3b))/ rbPz3b - 
			(a*cosb*cotb*y3b)/rbP2)* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) - 
			(y2*y3bMa* (-(cotb*pr3) + (y1*(2.*pr + aQrb))/ rbPy3b + 
			(a*y1)/rbP2)* Dy3bDy3* DrbDy3b)/ (rbPy3b*rbP2));
			 
		dg_ic[0][1][0] +=
			2.*(pr3*(-((cosb*(Dz3bDy1 + DrbDy1)*
			(1. - 2.*pr + 2.*cotbP2*(1. - pr)))/ rbPz3b) + 
			((-pr + 2.*cotbP2*(1. - pr))* DrbDy1)/ rbPy3b) - 
			(y3bMa*(Dz3bDy1 + DrbDy1)* (cosbP2 - (a*cosb + cotb*pr3*z1b)/
			rb - (cosbP2*y2P2 - (a*cotb*rbTcosbPy3b* z1b)/rb)/ (rbPz3b*
			rb) + (a*cotb*y3b*z1b)/ rbP3))/ rbPz3bP2 - 
			(pr3*(cotb*(1. - 2.*pr - aQrb) - (y2P2*(pr + aQrb)*
			DrbDy1)/ rbPy3bP2 + (a*cotb*y1*DrbDy1)/ rbP2 - 
			(a*y2P2*DrbDy1)/ (rbPy3b* rbP2)))/ rbPy3b + 
			(y3bMa*((cotb*pr3)/rb - (3*a*y2P2*DrbDy1)/
			rbP4 - (a*y2P2*DrbDy1)/ (rbPy3b* rbP3) - 
			((-a + cotb*y1*pr3)* DrbDy1)/ rbP2 - (y2P2*(2.*pr + aQrb)*
			DrbDy1)/ (rbPy3b* rbP2) - (y2P2*(2.*pr + aQrb)*
			DrbDy1)/ (rbPy3bP2* rb)))/ rbPy3b + (y3bMa*(((cosbP2*y2P2 - 
			(a*cotb*rbTcosbPy3b* z1b)/rb)* (Dz3bDy1 + DrbDy1))/
			(rbPz3bP2* rb) - (-((a*cotb* rbTcosbPy3b* Dz1bDy1)/
			rb) + (a*cotb*rbTcosbPy3b* z1b* DrbDy1)/ rbP2 - 
			(a*cosb*cotb*z1b* DrbDy1)/ rb)/ (rbPz3b* rb) - 
			(cotb*pr3*Dz1bDy1)/ rb + (a*cotb*y3b* Dz1bDy1)/ rbP3 + 
			((a*cosb + cotb*pr3*z1b)* DrbDy1)/ rbP2 + ((cosbP2*y2P2 - 
			(a*cotb*rbTcosbPy3b* z1b)/rb)* DrbDy1)/ (rbPz3b*
			rbP2) - (3*a*cotb*y3b*z1b* DrbDy1)/ rbP4))/ rbPz3b - 
			(a*cotb*y3bMa)/rbP3 + (cotb*pr3*(cosb + aQrb)* (Dz3bDy1 + 
			DrbDy1)* z1b)/ rbPz3bP2 - (cotb*pr3*(cosb + aQrb)*
			Dz1bDy1)/ rbPz3b + (pr3*(-a + (y2P2*(pr + aQrb))/ rbPy3b + 
			cotb*y1*(1. - 2.*pr - aQrb) + pr*y3b)* DrbDy1)/ rbPy3bP2 - 
			(y3bMa*(-2.*pr + (a*y2P2)/rbP3 + (-a + cotb*y1*pr3)/rb + 
			(y2P2*(2.*pr + aQrb))/ (rbPy3b* rb))* DrbDy1)/ rbPy3bP2 + 
			(3*a*cotb*y1*y3bMa* DrbDy1)/ rbP4 + (a*cotb*pr3*z1b*
			DrbDy1)/ (rbPz3b* rbP2)); 
			
		dg_ic[0][1][1] +=
			2.*(pr3*(((-pr + 2.*cotbP2*(1. - pr))* DrbDy2)/
			rbPy3b - (cosb*(1. - 2.*pr + 2.*cotbP2*(1. - pr))* DrbDy2)/
			rbPz3b) - (pr3*((2.*y2*(pr + aQrb))/ rbPy3b - 
			(y2P2*(pr + aQrb)* DrbDy2)/ rbPy3bP2 + (a*cotb*y1*DrbDy2)/
			rbP2 - (a*y2P2*DrbDy2)/ (rbPy3b* rbP2)))/ rbPy3b + 
			(y3bMa*(-((2.*cosbP2*y2 + (a*cotb*rbTcosbPy3b*
			z1b* DrbDy2)/ rbP2 - (a*cosb*cotb*z1b* DrbDy2)/ rb)/
			(rbPz3b* rb)) + ((a*cosb + cotb*pr3*z1b)* DrbDy2)/ rbP2 + 
			((cosbP2*y2P2 - (a*cotb*rbTcosbPy3b* z1b)/rb)* DrbDy2)/
			(rbPz3b* rbP2) + ((cosbP2*y2P2 - (a*cotb*rbTcosbPy3b*
			z1b)/rb)* DrbDy2)/ (rbPz3bP2* rb) - (3*a*cotb*y3b*z1b*
			DrbDy2)/ rbP4))/ rbPz3b + (y3bMa*((2.*a*y2)/rbP3 + 
			(2.*y2*(2.*pr + aQrb))/ (rbPy3b* rb) - (3*a*y2P2*DrbDy2)/
			rbP4 - (a*y2P2*DrbDy2)/ (rbPy3b* rbP3) - ((-a + cotb*y1*pr3)*
			DrbDy2)/ rbP2 - (y2P2*(2.*pr + aQrb)* DrbDy2)/ (rbPy3b*
			rbP2) - (y2P2*(2.*pr + aQrb)* DrbDy2)/ (rbPy3bP2* rb)))/
			rbPy3b + (pr3*(-a + (y2P2*(pr + aQrb))/ rbPy3b + 
			cotb*y1*(1. - 2.*pr - aQrb) + pr*y3b)* DrbDy2)/ rbPy3bP2 - 
			(y3bMa*(cosbP2 - (a*cosb + cotb*pr3*z1b)/ rb - 
			(cosbP2*y2P2 - (a*cotb*rbTcosbPy3b* z1b)/rb)/ (rbPz3b*
			rb) + (a*cotb*y3b*z1b)/ rbP3)* DrbDy2)/ rbPz3bP2 - 
			(y3bMa*(-2.*pr + (a*y2P2)/rbP3 + (-a + cotb*y1*pr3)/rb + 
			(y2P2*(2.*pr + aQrb))/ (rbPy3b* rb))* DrbDy2)/ rbPy3bP2 + 
			(3*a*cotb*y1*y3bMa* DrbDy2)/ rbP4 + (cotb*pr3*(cosb + aQrb)*
			z1b* DrbDy2)/ rbPz3bP2 + (a*cotb*pr3*z1b* DrbDy2)/
			(rbPz3b* rbP2));
			
		dg_ic[0][1][2] +=
			2.*(pr3*(((-pr + 2.*cotbP2*(1. - pr))* (Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ rbPy3b - (cosb*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b)* (1. - 2.*pr + 2.*cotbP2*(1. - pr)))/ rbPz3b) + 
			(pr3*(Dy3bDy3 + Dy3bDy3* DrbDy3b)* (-a + (y2P2*(pr + aQrb))/
			rbPy3b + cotb*y1*(1. - 2.*pr - aQrb) + 
			pr*y3b))/rbPy3bP2 - (y3bMa*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b)* (cosbP2 - (a*cosb + cotb*pr3*z1b)/ rb - 
			(cosbP2*y2P2 - (a*cotb*rbTcosbPy3b* z1b)/rb)/ (rbPz3b* rb) + 
			(a*cotb*y3b*z1b)/ rbP3))/ rbPz3bP2 - (y3bMa*(Dy3bDy3 + 
			Dy3bDy3* DrbDy3b)* (-2.*pr + (a*y2P2)/rbP3 + 
			(-a + cotb*y1*pr3)/rb + (y2P2*(2.*pr + aQrb))/ (rbPy3b*
			rb)))/ rbPy3bP2 - (pr3*(-((y2P2*(pr + aQrb)* (Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ rbPy3bP2) + pr*Dy3bDy3 + (a*cotb*y1*Dy3bDy3*
			DrbDy3b)/ rbP2 - (a*y2P2*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP2)))/
			rbPy3b + (y3bMa*(-((y2P2*(2.*pr + aQrb)* (Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ (rbPy3bP2* rb)) - (3*a*y2P2*Dy3bDy3* DrbDy3b)/
			rbP4 - (a*y2P2*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP3) - 
			((-a + cotb*y1*pr3)* Dy3bDy3* DrbDy3b)/ rbP2 - 
			(y2P2*(2.*pr + aQrb)* Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP2)))/
			rbPy3b + (y3bMa*(((cosbP2*y2P2 - (a*cotb*rbTcosbPy3b*
			z1b)/rb)* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ (rbPz3bP2*
			rb) - (-((a*cotb* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b)* z1b)/rb) - 
			(a*cotb*rbTcosbPy3b* Dz1bDy3b* Dy3bDy3)/rb + (a*cotb*
			rbTcosbPy3b* z1b*Dy3bDy3* DrbDy3b)/ rbP2)/ (rbPz3b* rb) + 
			(a*cotb*z1b* Dy3bDy3)/ rbP3 - (cotb*pr3* Dz1bDy3b*
			Dy3bDy3)/rb + (a*cotb*y3b* Dz1bDy3b* Dy3bDy3)/ rbP3 + 
			((a*cosb + cotb*pr3*z1b)* Dy3bDy3* DrbDy3b)/ rbP2 + 
			((cosbP2*y2P2 - (a*cotb*rbTcosbPy3b* z1b)/rb)* Dy3bDy3*
			DrbDy3b)/ (rbPz3b* rbP2) - (3*a*cotb*y3b*z1b* Dy3bDy3* DrbDy3b)/
			rbP4))/ rbPz3b + (cotb*pr3*(cosb + aQrb)* (Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b)* z1b)/ rbPz3bP2 + 
			((cosbP2 - (a*cosb + cotb*pr3*z1b)/ rb - (cosbP2*y2P2 - 
			(a*cotb*rbTcosbPy3b* z1b)/rb)/ (rbPz3b* rb) + 
			(a*cotb*y3b*z1b)/ rbP3)*Dy3bDy3)/rbPz3b + 
			((-2.*pr + (a*y2P2)/rbP3 + (-a + cotb*y1*pr3)/rb + 
			(y2P2*(2.*pr + aQrb))/ (rbPy3b* rb))* Dy3bDy3)/ rbPy3b - 
			(a*cotb*y1*Dy3bDy3)/rbP3 - (cotb*pr3*(cosb + aQrb)*
			Dz1bDy3b* Dy3bDy3)/ rbPz3b + (3*a*cotb*y1*y3bMa*Dy3bDy3*
			DrbDy3b)/ rbP4 + (a*cotb*pr3*z1b* Dy3bDy3* DrbDy3b)/
			(rbPz3b* rbP2));
			 
		dg_ic[0][2][0] +=
			2.*(2.*cotb*pr3*(1. - pr)*(DrbDy1* DfbDrb + Dz1bDy1*
			DfbDz1b + DfbDy1) + pr2*((cosb*y2*(cosb + aQrb)*
			(Dz3bDy1 + DrbDy1))/ rbPz3bP2 - (y2*(2.*pr + aQrb)*
			DrbDy1)/ rbPy3bP2 - (a*y2*DrbDy1)/ (rbPy3b* rbP2) + 
			(a*cosb*y2*DrbDy1)/ (rbPz3b* rbP2)) + (y2*y3bMa*
			((-2.*pr*DrbDy1)/ rbPy3bP2 - (2.*a*DrbDy1)/ rbP3))/rb - 
			(cosb*y2*y3bMa* (Dz3bDy1 + DrbDy1)*
			(1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b - 
			(a*y3b)/rbP2))/ (rbPz3bP2* rb) + (cosb*y2*y3bMa* (((cosb + aQrb)*
			rbTcosbPy3b* (Dz3bDy1 + DrbDy1))/ rbPz3bP2 - 
			(cosb*(cosb + aQrb)* DrbDy1)/ rbPz3b + (a*rbTcosbPy3b* DrbDy1)/
			(rbPz3b* rbP2) + (2.*a*y3b* DrbDy1)/ rbP3))/ (rbPz3b* rb) - 
			(y2*y3bMa* ((2.*pr)/rbPy3b + a/rbP2)* DrbDy1)/ rbP2 - 
			(cosb*y2*y3bMa* (1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b - (a*y3b)/rbP2)* DrbDy1)/ (rbPz3b* rbP2));
			
		dg_ic[0][2][1] +=
			2.*(2.*cotb*pr3*(1. - pr)*(DrbDy2* DfbDrb + DfbDy2) + 
			pr2*((2.*pr + aQrb)/ rbPy3b - (cosb*(cosb + aQrb))/ rbPz3b - 
			(y2*(2.*pr + aQrb)* DrbDy2)/ rbPy3bP2 + (cosb*y2*(cosb + aQrb)*
			DrbDy2)/ rbPz3bP2 - (a*y2*DrbDy2)/ (rbPy3b* rbP2) + 
			(a*cosb*y2*DrbDy2)/ (rbPz3b* rbP2)) + (y3bMa*((2.*pr)/rbPy3b +
			a/rbP2))/rb + (y2*y3bMa* ((-2.*pr*DrbDy2)/ rbPy3bP2 - 
			(2.*a*DrbDy2)/ rbP3))/rb + (cosb*y3bMa*
			(1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b - 
			(a*y3b)/rbP2))/ (rbPz3b* rb) + (cosb*y2*y3bMa*
			(-((cosb*(cosb + aQrb)* DrbDy2)/ rbPz3b) + ((cosb + aQrb)*
			rbTcosbPy3b* DrbDy2)/ rbPz3bP2 + (a*rbTcosbPy3b* DrbDy2)/
			(rbPz3b* rbP2) + (2.*a*y3b* DrbDy2)/ rbP3))/ (rbPz3b* rb) - 
			(y2*y3bMa* ((2.*pr)/rbPy3b + a/rbP2)* DrbDy2)/ rbP2 - 
			(cosb*y2*y3bMa* (1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b - (a*y3b)/rbP2)* DrbDy2)/ (rbPz3b* rbP2) - 
			(cosb*y2*y3bMa* (1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b - (a*y3b)/rbP2)* DrbDy2)/ (rbPz3bP2* rb)); 
			
		dg_ic[0][2][2] +=
			2.*(2.*cotb*pr3*(1. - pr)*(Dy3bDy3* DrbDy3b* DfbDrb + 
			Dz1bDy3b* Dy3bDy3* DfbDz1b) + pr2*(-((y2*(2.*pr + aQrb)*
			(Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPy3bP2) + 
			(cosb*y2*(cosb + aQrb)* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ rbPz3bP2 - (a*y2*Dy3bDy3* DrbDy3b)/ (rbPy3b*
			rbP2) + (a*cosb*y2*Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2)) + 
			(y2*y3bMa* ((-2.*pr*(Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2 - (2.*a*Dy3bDy3* DrbDy3b)/ rbP3))/rb - 
			(cosb*y2*y3bMa* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b)*
			(1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b - 
			(a*y3b)/rbP2))/ (rbPz3bP2* rb) + (cosb*y2*y3bMa*
			(-(((cosb + aQrb)* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/
			rbPz3b) + ((cosb + aQrb)* rbTcosbPy3b* (Dz3bDy3b*
			Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3bP2 - 
			(a*Dy3bDy3)/rbP2 + (a*rbTcosbPy3b* Dy3bDy3* DrbDy3b)/
			(rbPz3b* rbP2) + (2.*a*y3b*Dy3bDy3* DrbDy3b)/ rbP3))/
			(rbPz3b* rb) + (y2*((2.*pr)/rbPy3b + a/rbP2)*Dy3bDy3) /rb + 
			(cosb*y2*(1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b - 
			(a*y3b)/rbP2)* Dy3bDy3)/ (rbPz3b* rb) - (y2*y3bMa*
			((2.*pr)/rbPy3b + a/rbP2)*Dy3bDy3* DrbDy3b)/ rbP2 - 
			(cosb*y2*y3bMa* (1. - 2.*pr - ((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b - (a*y3b)/rbP2)* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2));
			
		dg_ic[1][0][0] +=
			2.*(pr3*(-((cosb*(1. + 2.*cotbP2*(1. - pr))*
			(Dz3bDy1 + DrbDy1))/ rbPz3b) + ((pr + 2.*cotbP2*(1. - pr))*
			DrbDy1)/ rbPy3b) + (cotb*pr3*(-((a*rbTsinbMy1)/ (cosb*rb)) + 
			cosb*z1b)* (Dz3bDy1 + DrbDy1))/ rbPz3bP2 - (cotb*y3bMa*
			(Dz3bDy1 + DrbDy1)* (-(cosb*sinb) + (rbTsinbMy1* (pr2*cosb - 
			((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b))/ rb + 
			(a*y1*y3b)/(cosb*rbP3)))/ rbPz3bP2 - (cotb*pr3*(-((a*
			(-1. + sinb*DrbDy1))/ (cosb*rb)) + cosb*Dz1bDy1 + (a*rbTsinbMy1*
			DrbDy1)/ (cosb*rbP2)))/ rbPz3b + (cotb*y3bMa*
			(((-1. + sinb*DrbDy1)* (pr2*cosb - ((1. + a/(cosb*rb))*
			rbTcosbPy3b)/ rbPz3b))/ rb + (rbTsinbMy1* (((1. + a/(cosb*rb))*
			rbTcosbPy3b* (Dz3bDy1 + DrbDy1))/ rbPz3bP2 - 
			(cosb*(1. + a/(cosb*rb))* DrbDy1)/ rbPz3b + (a*rbTcosbPy3b*
			DrbDy1)/ (cosb* rbPz3b* rbP2)))/rb + (a*y3b)/(cosb*rbP3) - 
			(rbTsinbMy1* (pr2*cosb - ((1. + a/(cosb*rb))* rbTcosbPy3b)/
			rbPz3b)* DrbDy1)/ rbP2 - (3*a*y1*y3b* DrbDy1)/ (cosb*rbP4)))/
			rbPz3b + (pr3*(-(cotb*pr3) + (2.*y1*(pr + aQrb))/ rbPy3b + 
			(a*cotb)/rb - (y1P2*(pr + aQrb)* DrbDy1)/ rbPy3bP2 - 
			(a*cotb*y1*DrbDy1)/ rbP2 - (a*y1P2*DrbDy1)/ (rbPy3b* rbP2)))/
			rbPy3b + (y3bMa*((-2.*a*y1)/rbP3 + (cotb*pr3)/rb - 
			(2.*y1*(2.*pr + aQrb))/ (rbPy3b* rb) + (3*a*y1P2*DrbDy1)/ rbP4 + 
			(a*y1P2*DrbDy1)/ (rbPy3b* rbP3) - ((a + cotb*y1*pr3)* DrbDy1)/
			rbP2 + (y1P2*(2.*pr + aQrb)* DrbDy1)/ (rbPy3b* rbP2) + 
			(y1P2*(2.*pr + aQrb)* DrbDy1)/ (rbPy3bP2* rb)))/ rbPy3b - 
			(a*cotb*y3bMa)/rbP3 - (y3bMa*(2.*pr - (a*y1P2)/rbP3 + 
			(a + cotb*y1*pr3)/rb - (y1P2*(2.*pr + aQrb))/ (rbPy3b* rb))*
			DrbDy1)/ rbPy3bP2 - (pr3*(-a - cotb*y1*pr3 + (y1P2*(pr + aQrb))/
			rbPy3b + (a*cotb*y1)/rb + pr*y3b)* DrbDy1)/ rbPy3bP2 + 
			(3*a*cotb*y1*y3bMa* DrbDy1)/ rbP4);
			
		dg_ic[1][0][1] +=
			2.*(pr3*(((pr + 2.*cotbP2*(1. - pr))* DrbDy2)/ rbPy3b - 
			(cosb*(1. + 2.*cotbP2*(1. - pr))* DrbDy2)/ rbPz3b) - 
			(cotb*pr3*((a*rbTsinbMy1* DrbDy2)/ (cosb*rbP2) - (a*sinb*DrbDy2)/
			(cosb*rb)))/ rbPz3b + (pr3*(-((y1P2*(pr + aQrb)* DrbDy2)/
			rbPy3bP2) - (a*cotb*y1*DrbDy2)/ rbP2 - (a*y1P2*DrbDy2)/ (rbPy3b*
			rbP2)))/ rbPy3b + (cotb*y3bMa* ((rbTsinbMy1*
			(-((cosb*(1. + a/(cosb*rb))* DrbDy2)/ rbPz3b) + 
			((1. + a/(cosb*rb))* rbTcosbPy3b* DrbDy2)/ rbPz3bP2 + 
			(a*rbTcosbPy3b* DrbDy2)/ (cosb* rbPz3b* rbP2)))/rb - (rbTsinbMy1*
			(pr2*cosb - ((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b)* DrbDy2)/
			rbP2 + (sinb* (pr2*cosb - ((1. + a/(cosb*rb))* rbTcosbPy3b)/
			rbPz3b)* DrbDy2)/ rb - (3*a*y1*y3b* DrbDy2)/ (cosb*rbP4)))/
			rbPz3b + (y3bMa*((3*a*y1P2* DrbDy2)/ rbP4 + (a*y1P2*DrbDy2)/
			(rbPy3b* rbP3) - ((a + cotb*y1*pr3)* DrbDy2)/ rbP2 + 
			(y1P2*(2.*pr + aQrb)* DrbDy2)/ (rbPy3b* rbP2) + 
			(y1P2*(2.*pr + aQrb)* DrbDy2)/ (rbPy3bP2* rb)))/ rbPy3b + 
			(cotb*pr3*(-((a*rbTsinbMy1)/ (cosb*rb)) + cosb*z1b)* DrbDy2)/
			rbPz3bP2 - (cotb*y3bMa* (-(cosb*sinb) + (rbTsinbMy1* (pr2*cosb - 
			((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b))/ rb + 
			(a*y1*y3b)/(cosb*rbP3))* DrbDy2)/ rbPz3bP2 - 
			(y3bMa*(2.*pr - (a*y1P2)/rbP3 + (a + cotb*y1*pr3)/rb - 
			(y1P2*(2.*pr + aQrb))/ (rbPy3b* rb))* DrbDy2)/ rbPy3bP2 - 
			(pr3*(-a - cotb*y1*pr3 + (y1P2*(pr + aQrb))/ rbPy3b + 
			(a*cotb*y1)/rb + pr*y3b)* DrbDy2)/ rbPy3bP2 + (3*a*cotb*y1*y3bMa*
			DrbDy2)/ rbP4);
			
		dg_ic[1][0][2] +=
			2.*(pr3*(((pr + 2.*cotbP2*(1. - pr))* (Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ rbPy3b - (cosb*(1. + 2.*cotbP2*(1. - pr))*
			(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3b) + 
			(cotb*pr3*(-((a*rbTsinbMy1)/ (cosb*rb)) + cosb*z1b)*
			(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3bP2 - 
			(cotb*y3bMa* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b)*
			(-(cosb*sinb) + (rbTsinbMy1* (pr2*cosb - ((1. + a/(cosb*rb))*
			rbTcosbPy3b)/ rbPz3b))/ rb + (a*y1*y3b)/(cosb*rbP3)))/
			rbPz3bP2 - (cotb*pr3*(cosb*Dz1bDy3b* Dy3bDy3 + (a*rbTsinbMy1*
			Dy3bDy3* DrbDy3b)/ (cosb*rbP2) - (a*sinb*Dy3bDy3*
			DrbDy3b)/ (cosb*rb)))/ rbPz3b - (y3bMa*(Dy3bDy3 + Dy3bDy3*
			DrbDy3b)* (2.*pr - (a*y1P2)/rbP3 + (a + cotb*y1*pr3)/rb - 
			(y1P2*(2.*pr + aQrb))/ (rbPy3b* rb)))/ rbPy3bP2 + 
			(pr3*(-((y1P2*(pr + aQrb)* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2) + pr*Dy3bDy3 - (a*cotb*y1*Dy3bDy3* DrbDy3b)/
			rbP2 - (a*y1P2*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP2)))/ rbPy3b - 
			(pr3*(Dy3bDy3 + Dy3bDy3* DrbDy3b)* (-a - cotb*y1*pr3 + 
			(y1P2*(pr + aQrb))/ rbPy3b + (a*cotb*y1)/rb + pr*y3b))/ rbPy3bP2 + 
			(y3bMa*((y1P2*(2.*pr + aQrb)* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			(rbPy3bP2* rb) + (3*a*y1P2*Dy3bDy3* DrbDy3b)/ rbP4 + 
			(a*y1P2*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP3) - ((a + cotb*y1*pr3)*
			Dy3bDy3* DrbDy3b)/ rbP2 + (y1P2*(2.*pr + aQrb)* Dy3bDy3*
			DrbDy3b)/ (rbPy3b* rbP2)))/ rbPy3b + (cotb*y3bMa* ((rbTsinbMy1*
			(-(((1. + a/(cosb*rb))* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/
			rbPz3b) + ((1. + a/(cosb*rb))* rbTcosbPy3b* (Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ rbPz3bP2 + (a*rbTcosbPy3b* Dy3bDy3* DrbDy3b)/
			(cosb* rbPz3b* rbP2)))/rb + (a*y1*Dy3bDy3)/ (cosb*rbP3) - 
			(rbTsinbMy1* (pr2*cosb - ((1. + a/(cosb*rb))* rbTcosbPy3b)/
			rbPz3b)* Dy3bDy3* DrbDy3b)/ rbP2 + (sinb* (pr2*cosb - 
			((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b)* Dy3bDy3* DrbDy3b)/
			rb - (3*a*y1*y3b*Dy3bDy3* DrbDy3b)/ (cosb*rbP4)))/ rbPz3b + 
			(cotb*(-(cosb*sinb) + (rbTsinbMy1* (pr2*cosb - 
			((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b))/ rb + 
			(a*y1*y3b)/(cosb*rbP3))* Dy3bDy3)/ rbPz3b + 
			((2.*pr - (a*y1P2)/rbP3 + (a + cotb*y1*pr3)/rb - 
			(y1P2*(2.*pr + aQrb))/ (rbPy3b* rb))* Dy3bDy3)/ rbPy3b - 
			(a*cotb*y1*Dy3bDy3)/rbP3 + (3*a*cotb*y1*y3bMa*Dy3bDy3*
			DrbDy3b)/ rbP4);
			 
		dg_ic[1][1][0] +=
			2.*((cotb*y2*pr3*(1. + a/(cosb*rb))* (Dz3bDy1 + DrbDy1))/
			rbPz3bP2 + 2.*cotbP2*pr3*(1. - pr)* (DrbDy1* DfbDrb + 
			Dz1bDy1* DfbDz1b + DfbDy1) + (y2*pr3*((pr + aQrb)/
			rbPy3b - (y1*(pr + aQrb)* DrbDy1)/ rbPy3bP2 - (a*cotb*DrbDy1)/
			rbP2 - (a*y1*DrbDy1)/ (rbPy3b* rbP2)))/ rbPy3b - (cotb*y2*y3bMa*
			(Dz3bDy1 + DrbDy1)* (-pr2*cosb + ((1. + a/(cosb*rb))*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/(cosb*rbP2)))/ (rbPz3bP2*
			rb) + (cotb*y2*y3bMa* (-(((1. + a/(cosb*rb))* rbTcosbPy3b*
			(Dz3bDy1 + DrbDy1))/ rbPz3bP2) + (cosb*(1. + a/(cosb*rb))*
			DrbDy1)/ rbPz3b - (a*rbTcosbPy3b* DrbDy1)/ (cosb* rbPz3b*
			rbP2) - (2.*a*y3b* DrbDy1)/ (cosb*rbP3)))/ (rbPz3b* rb) + 
			(y2*y3bMa* ((-2.*pr)/rbPy3b - (a*(rbPy3bM1 + rbM1))/ rb - 
			(a*y1* (-(DrbDy1/ rbPy3bP2) - DrbDy1/ rbP2))/rb +
			(2.*pr*y1*DrbDy1) /rbPy3bP2 + (a*y1* (rbPy3bM1 + rbM1)* DrbDy1)/
			rbP2))/ (rbPy3b*rb) - (y2*pr3*((y1*(pr + aQrb))/ rbPy3b - 
			cotb*(1. - 2.*pr - aQrb))* DrbDy1)/ rbPy3bP2 + 
			(3*a*cotb*y2*y3bMa* DrbDy1)/ rbP4 + (a*cotb*y2*pr3*DrbDy1)/
			(cosb*rbPz3b* rbP2) - (y2*y3bMa* (cotb*pr3 - (2.*pr*y1)/rbPy3b - 
			(a*y1* (rbPy3bM1 + rbM1))/ rb)* DrbDy1)/ (rbPy3b*rbP2) - 
			(cotb*y2*y3bMa* (-pr2*cosb + ((1. + a/(cosb*rb))* rbTcosbPy3b)/
			rbPz3b + (a*y3b)/(cosb*rbP2))* DrbDy1)/ (rbPz3b* rbP2) - 
			(y2*y3bMa* (cotb*pr3 - (2.*pr*y1)/rbPy3b - (a*y1* (rbPy3bM1 + 
			rbM1))/ rb)* DrbDy1)/ (rbPy3bP2*rb));
		  
		dg_ic[1][1][1] +=
			2.*((pr3*((y1*(pr + aQrb))/ rbPy3b - cotb*(1. - 2.*pr - aQrb)))/
			rbPy3b - (cotb*pr3*(1. + a/(cosb*rb)))/ rbPz3b + 
			2.*cotbP2*pr3*(1. - pr)* (DrbDy2* DfbDrb + DfbDy2) + 
			(y2*pr3*(-((y1*(pr + aQrb)* DrbDy2)/ rbPy3bP2) - (a*cotb*DrbDy2)/
			rbP2 - (a*y1*DrbDy2)/ (rbPy3b* rbP2)))/ rbPy3b - 
			(a*cotb*y3bMa)/rbP3 + (y3bMa*(cotb*pr3 - (2.*pr*y1)/rbPy3b - 
			(a*y1* (rbPy3bM1 + rbM1))/ rb))/ (rbPy3b*rb) + (cotb*y3bMa*
			(-pr2*cosb + ((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b + 
			(a*y3b)/(cosb*rbP2)))/ (rbPz3b* rb) + (y2*y3bMa*
			(-((a*y1*(-(DrbDy2/ rbPy3bP2) - DrbDy2/
			rbP2))/rb) + (2.*pr*y1*DrbDy2)/ rbPy3bP2 + (a*y1* (rbPy3bM1 + 
			rbM1)* DrbDy2)/ rbP2))/ (rbPy3b*rb) + (cotb*y2*y3bMa*
			((cosb*(1. + a/(cosb*rb))* DrbDy2)/ rbPz3b - ((1. + a/(cosb*rb))*
			rbTcosbPy3b* DrbDy2)/ rbPz3bP2 - (a*rbTcosbPy3b* DrbDy2)/ (cosb*
			rbPz3b* rbP2) - (2.*a*y3b* DrbDy2)/ (cosb*rbP3)))/ (rbPz3b*
			rb) - (y2*pr3*((y1*(pr + aQrb))/ rbPy3b - 
			cotb*(1. - 2.*pr - aQrb))* DrbDy2)/ rbPy3bP2 + 
			(cotb*y2*pr3*(1. + a/(cosb*rb))* DrbDy2)/ rbPz3bP2 + 
			(3*a*cotb*y2*y3bMa* DrbDy2)/ rbP4 + (a*cotb*y2*pr3*DrbDy2)/
			(cosb*rbPz3b* rbP2) - (y2*y3bMa* (cotb*pr3 - (2.*pr*y1)/rbPy3b - 
			(a*y1* (rbPy3bM1 + rbM1))/ rb)* DrbDy2)/ (rbPy3b*rbP2) - 
			(cotb*y2*y3bMa* (-pr2*cosb + ((1. + a/(cosb*rb))*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/(cosb*rbP2))* DrbDy2)/
			(rbPz3b* rbP2) - (y2*y3bMa* (cotb*pr3 - (2.*pr*y1)/rbPy3b - 
			(a*y1* (rbPy3bM1 + rbM1))/ rb)* DrbDy2)/ (rbPy3bP2*rb) - 
			(cotb*y2*y3bMa* (-pr2*cosb + ((1. + a/(cosb*rb))* rbTcosbPy3b)/
			rbPz3b + (a*y3b)/(cosb*rbP2))* DrbDy2)/ (rbPz3bP2*
			rb));
			
		dg_ic[1][1][2] +=
			2.*(-((y2*pr3*((y1*(pr + aQrb))/ rbPy3b - 
			cotb*(1. - 2.*pr - aQrb))* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2) + (cotb*y2*pr3*(1. + a/(cosb*rb))* (Dz3bDy3b*
			Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3bP2 + 2.*cotbP2*pr3*(1. - pr)*
			(Dy3bDy3* DrbDy3b* DfbDrb + Dz1bDy3b* Dy3bDy3* DfbDz1b) + 
			(y2*pr3*(-((y1*(pr + aQrb)* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2) - (a*cotb*Dy3bDy3* DrbDy3b)/ rbP2 - (a*y1*Dy3bDy3*
			DrbDy3b)/ (rbPy3b* rbP2)))/ rbPy3b - (y2*y3bMa* (Dy3bDy3 + 
			Dy3bDy3* DrbDy3b)* (cotb*pr3 - (2.*pr*y1)/rbPy3b - (a*y1*
			(rbPy3bM1 + rbM1))/ rb))/ (rbPy3bP2*rb) - (cotb*y2*y3bMa*
			(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b)* (-pr2*cosb + 
			((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b + (a*y3b)/(cosb*rbP2)))/
			(rbPz3bP2* rb) + (y2*y3bMa* ((2.*pr*y1*(Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ rbPy3bP2 - (a*y1* (-((Dy3bDy3 + Dy3bDy3* DrbDy3b)/
			rbPy3bP2) - (Dy3bDy3* DrbDy3b)/ rbP2))/rb + (a*y1* (rbPy3bM1 + 
			rbM1)* Dy3bDy3* DrbDy3b)/ rbP2))/ (rbPy3b*rb) + (cotb*y2*y3bMa*
			(((1. + a/(cosb*rb))* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/
			rbPz3b - ((1. + a/(cosb*rb))* rbTcosbPy3b* (Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ rbPz3bP2 + (a*Dy3bDy3)/ (cosb*rbP2) - 
			(a*rbTcosbPy3b* Dy3bDy3* DrbDy3b)/ (cosb* rbPz3b* rbP2) - 
			(2.*a*y3b*Dy3bDy3* DrbDy3b)/ (cosb*rbP3)))/ (rbPz3b* rb) - 
			(a*cotb*y2*Dy3bDy3)/rbP3 + (y2*(cotb*pr3 - (2.*pr*y1)/rbPy3b - 
			(a*y1* (rbPy3bM1 + rbM1))/ rb)*Dy3bDy3)/ (rbPy3b*rb) + 
			(cotb*y2*(-pr2*cosb + ((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b + 
			(a*y3b)/(cosb*rbP2))* Dy3bDy3)/ (rbPz3b* rb) + 
			(3*a*cotb*y2*y3bMa*Dy3bDy3* DrbDy3b)/ rbP4 + 
			(a*cotb*y2*pr3*Dy3bDy3* DrbDy3b)/ (cosb*rbPz3b* rbP2) - (y2*y3bMa*
			(cotb*pr3 - (2.*pr*y1)/rbPy3b - (a*y1* (rbPy3bM1 + rbM1))/
			rb)*Dy3bDy3* DrbDy3b)/ (rbPy3b*rbP2) - (cotb*y2*y3bMa*
			(-pr2*cosb + ((1. + a/(cosb*rb))* rbTcosbPy3b)/ rbPz3b + 
			(a*y3b)/(cosb*rbP2))* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2));
			 
		dg_ic[1][2][0] +=
			2.*(-2.*cotb*pr3*(1. - pr)*(-((cosb* (Dz3bDy1 + DrbDy1))/
			rbPz3b) + DrbDy1/ rbPy3b) - (pr2*(2.*pr + aQrb))/ rbPy3b + 
			(y3bMa*(Dz3bDy1 + DrbDy1)* (cosb*sinb + (cotb* (pr2*cosb - 
			rbTcosbPy3b/ rbPz3b)* rbTcosbPy3b)/ rb + (a*(sinb - (rbTcosbPy3b*
			z1b)/ (rbPz3b* rb) - (y3b*z1b)/ rbP2))/rb)) /rbPz3bP2 - 
			(y3bMa*((cotb* ((rbTcosbPy3b* (Dz3bDy1 + DrbDy1))/ rbPz3bP2 - 
			(cosb*DrbDy1)/ rbPz3b)* rbTcosbPy3b)/ rb + (a*((rbTcosbPy3b*
			(Dz3bDy1 + DrbDy1)* z1b)/ (rbPz3bP2* rb) - (rbTcosbPy3b*
			Dz1bDy1)/ (rbPz3b* rb) - (y3b*Dz1bDy1)/ rbP2 + (rbTcosbPy3b*
			z1b* DrbDy1)/ (rbPz3b* rbP2) - (cosb*z1b* DrbDy1)/ (rbPz3b*
			rb) + (2.*y3b*z1b* DrbDy1)/ rbP3))/rb - (cotb* (pr2*cosb - 
			rbTcosbPy3b/ rbPz3b)* rbTcosbPy3b* DrbDy1)/ rbP2 - (a*(sinb - 
			(rbTcosbPy3b* z1b)/ (rbPz3b* rb) - (y3b*z1b)/ rbP2)* DrbDy1)/
			rbP2 + (cosb*cotb* (pr2*cosb - rbTcosbPy3b/ rbPz3b)* DrbDy1)/
			rb))/ rbPz3b + (y3bMa*((-2.*pr)/rbPy3b - a/rbP2 + 
			(2.*pr*y1*DrbDy1)/ rbPy3bP2 + (2.*a*y1*DrbDy1)/ rbP3))/rb - 
			(pr2*(cosb + aQrb)* (Dz3bDy1 + DrbDy1)* z1b)/ rbPz3bP2 + 
			(pr2*(cosb + aQrb)* Dz1bDy1)/ rbPz3b + 
			(2.*y1*(1. - pr)*(2.*pr + aQrb)* DrbDy1)/ rbPy3bP2 + 
			(2.*a*y1*(1. - pr)*DrbDy1)/ (rbPy3b*rbP2) - (y3bMa*(cotb*pr3 - 
			(2.*pr*y1)/rbPy3b - (a*y1)/rbP2)* DrbDy1)/ rbP2 - 
			(2.*a*(1. - pr)*z1b* DrbDy1)/ (rbPz3b* rbP2));
			
		dg_ic[1][2][1] +=
			2.*(-2.*cotb*pr3*(1. - pr)*(DrbDy2/ rbPy3b - (cosb*DrbDy2)/
			rbPz3b) - (y3bMa*((cotb* rbTcosbPy3b* (-((cosb*DrbDy2)
			/rbPz3b) + (rbTcosbPy3b* DrbDy2)/ rbPz3bP2))/ rb + 
			(a*((rbTcosbPy3b* z1b* DrbDy2)/ (rbPz3b* rbP2) - (cosb*z1b*
			DrbDy2)/ (rbPz3b* rb) + (rbTcosbPy3b* z1b* DrbDy2)/ (rbPz3bP2*
			rb) + (2.*y3b*z1b* DrbDy2)/ rbP3))/rb - (cotb* (pr2*cosb - 
			rbTcosbPy3b/ rbPz3b)* rbTcosbPy3b* DrbDy2)/ rbP2 - (a*(sinb - 
			(rbTcosbPy3b* z1b)/ (rbPz3b* rb) - (y3b*z1b)/ rbP2)* DrbDy2)/
			rbP2 + (cosb*cotb* (pr2*cosb - rbTcosbPy3b/ rbPz3b)* DrbDy2)/
			rb))/ rbPz3b + (y3bMa*((2.*pr*y1* DrbDy2)/ rbPy3bP2 + 
			(2.*a*y1*DrbDy2)/ rbP3))/rb + (2.*y1*(1. - pr)*(2.*pr + aQrb)*
			DrbDy2)/ rbPy3bP2 + (y3bMa*(cosb*sinb + (cotb* (pr2*cosb - 
			rbTcosbPy3b/ rbPz3b)* rbTcosbPy3b)/ rb + (a*(sinb - (rbTcosbPy3b*
			z1b)/ (rbPz3b* rb) - (y3b*z1b)/ rbP2))/rb)* DrbDy2)/ rbPz3bP2 + 
			(2.*a*y1*(1. - pr)*DrbDy2)/ (rbPy3b*rbP2) - (y3bMa*(cotb*pr3 - 
			(2.*pr*y1)/rbPy3b - (a*y1)/rbP2)* DrbDy2)/ rbP2 - 
			(pr2*(cosb + aQrb)*z1b* DrbDy2)/ rbPz3bP2 - (2.*a*(1. - pr)*z1b*
			DrbDy2)/ (rbPz3b* rbP2)); 
			
		dg_ic[1][2][2] +=
			2.*(-2.*cotb*pr3*(1. - pr)*((Dy3bDy3 + Dy3bDy3* DrbDy3b)/
			rbPy3b - (cosb*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPz3b) + (2.*y1*(1. - pr)*(2.*pr + aQrb)* (Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ rbPy3bP2 + (y3bMa*(Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b)* (cosb*sinb + (cotb* (pr2*cosb - rbTcosbPy3b/
			rbPz3b)* rbTcosbPy3b)/ rb + (a*(sinb - (rbTcosbPy3b* z1b)/
			(rbPz3b* rb) - (y3b*z1b)/ rbP2))/rb)) /rbPz3bP2 - (y3bMa*((cotb*
			(-((Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b)/ rbPz3b) + (rbTcosbPy3b*
			(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPz3bP2)* rbTcosbPy3b)/
			rb + (cotb* (pr2*cosb - rbTcosbPy3b/ rbPz3b)* (Dy3bDy3 + 
			cosb*Dy3bDy3* DrbDy3b))/ rb + (a*(-(((Dy3bDy3 + cosb*Dy3bDy3*
			DrbDy3b)* z1b)/ (rbPz3b* rb)) + (rbTcosbPy3b* (Dz3bDy3b*
			Dy3bDy3 + Dy3bDy3* DrbDy3b)* z1b)/ (rbPz3bP2* rb) - (z1b*Dy3bDy3)/
			rbP2 - (rbTcosbPy3b* Dz1bDy3b* Dy3bDy3)/ (rbPz3b* rb) - 
			(y3b*Dz1bDy3b* Dy3bDy3)/ rbP2 + (rbTcosbPy3b* z1b*Dy3bDy3*
			DrbDy3b)/ (rbPz3b* rbP2) + (2.*y3b*z1b* Dy3bDy3* DrbDy3b)/
			rbP3))/rb - (cotb* (pr2*cosb - rbTcosbPy3b/ rbPz3b)* rbTcosbPy3b*
			Dy3bDy3* DrbDy3b)/ rbP2 - (a*(sinb - (rbTcosbPy3b* z1b)/ (rbPz3b*
			rb) - (y3b*z1b)/ rbP2)* Dy3bDy3* DrbDy3b)/ rbP2))/ rbPz3b + 
			(y3bMa*((2.*pr*y1* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/ rbPy3bP2 + 
			(2.*a*y1*Dy3bDy3* DrbDy3b)/ rbP3))/rb - (pr2*(cosb + aQrb)*
			(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b)* z1b)/ rbPz3bP2 - 
			((cosb*sinb + (cotb*(pr2*cosb - rbTcosbPy3b/ rbPz3b)* rbTcosbPy3b)/
			rb + (a*(sinb - (rbTcosbPy3b* z1b)/ (rbPz3b* rb) - (y3b*z1b)/
			rbP2))/rb)* Dy3bDy3)/ rbPz3b + 
			((cotb*pr3 - (2.*pr*y1)/rbPy3b - (a*y1)/rbP2)* Dy3bDy3)/rb + 
			(pr2*(cosb + aQrb)* Dz1bDy3b* Dy3bDy3)/ rbPz3b + 
			(2.*a*y1*(1. - pr)*Dy3bDy3* DrbDy3b)/ (rbPy3b*rbP2) - 
			(y3bMa*(cotb*pr3 - (2.*pr*y1)/rbPy3b - (a*y1)/rbP2)* Dy3bDy3*
			DrbDy3b)/ rbP2 - (2.*a*(1. - pr)*z1b*Dy3bDy3* DrbDy3b)/
			(rbPz3b* rbP2));
			
		dg_ic[2][0][0] +=
			2.*(pr3*((cosb*y2*(cosb + aQrb)* (Dz3bDy1 + DrbDy1))/
			rbPz3bP2 - (y2*(1. + aQrb)* DrbDy1)/ rbPy3bP2 - (a*y2*DrbDy1)/
			(rbPy3b* rbP2) + (a*cosb*y2*DrbDy1)/ (rbPz3b* rbP2)) - 
			(cosb*y2*y3bMa* (((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b + 
			(a*y3b)/rbP2)* (Dz3bDy1 + DrbDy1))/ (rbPz3bP2* rb) - (y2*y3bMa*
			(-(DrbDy1/ rbPy3bP2) - (2.*a*DrbDy1)/ rbP3))/rb + (cosb*y2*y3bMa*
			(-(((cosb + aQrb)* rbTcosbPy3b* (Dz3bDy1 + DrbDy1))/ rbPz3bP2) +
			(cosb*(cosb + aQrb)* DrbDy1)/ rbPz3b - (a*rbTcosbPy3b*
			DrbDy1)/ (rbPz3b* rbP2) - (2.*a*y3b* DrbDy1)/ rbP3))/ (rbPz3b*
			rb) + (y2*y3bMa* (rbPy3bM1 + a/rbP2)* DrbDy1)/ rbP2 - 
			(cosb*y2*y3bMa* (((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b + 
			(a*y3b)/rbP2)* DrbDy1)/ (rbPz3b* rbP2));
			
		dg_ic[2][0][1] +=
			2.*(pr3*((1. + aQrb)/ rbPy3b - (cosb*(cosb + aQrb))/ rbPz3b - 
			(y2*(1. + aQrb)* DrbDy2)/ rbPy3bP2 + (cosb*y2*(cosb + aQrb)*
			DrbDy2)/ rbPz3bP2 - (a*y2*DrbDy2)/ (rbPy3b* rbP2) + 
			(a*cosb*y2*DrbDy2)/ (rbPz3b* rbP2)) - (y3bMa*(rbPy3bM1 + 
			a/rbP2))/rb + (cosb*y3bMa* (((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b + (a*y3b)/rbP2))/ (rbPz3b* rb) - (y2*y3bMa* (-(DrbDy2/
			rbPy3bP2) - (2.*a*DrbDy2)/ rbP3))/rb + (cosb*y2*y3bMa*
			((cosb*(cosb + aQrb)* DrbDy2)/ rbPz3b - ((cosb + aQrb)*
			rbTcosbPy3b* DrbDy2)/ rbPz3bP2 - (a*rbTcosbPy3b* DrbDy2)/
			(rbPz3b* rbP2) - (2.*a*y3b* DrbDy2)/ rbP3))/ (rbPz3b* rb) + 
			(y2*y3bMa* (rbPy3bM1 + a/rbP2)* DrbDy2)/ rbP2 - (cosb*y2*y3bMa*
			(((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)* DrbDy2)/
			(rbPz3b* rbP2) - (cosb*y2*y3bMa* (((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b + (a*y3b)/rbP2)* DrbDy2)/ (rbPz3bP2* rb));
			
		dg_ic[2][0][2] +=
			2.*(pr3*(-((y2*(1. + aQrb)* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2) + (cosb*y2*(cosb + aQrb)* (Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ rbPz3bP2 - (a*y2*Dy3bDy3* DrbDy3b)/ (rbPy3b*
			rbP2) + (a*cosb*y2*Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2)) - (y2*y3bMa*
			(-((Dy3bDy3 + Dy3bDy3* DrbDy3b)/ rbPy3bP2) - (2.*a*Dy3bDy3*
			DrbDy3b)/ rbP3))/rb - (cosb*y2*y3bMa* (((cosb + aQrb)*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)* (Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ (rbPz3bP2* rb) + (cosb*y2*y3bMa*
			(((cosb + aQrb)* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/ rbPz3b - 
			((cosb + aQrb)* rbTcosbPy3b* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ rbPz3bP2 + (a*Dy3bDy3)/rbP2 - (a*rbTcosbPy3b*
			Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) - (2.*a*y3b*Dy3bDy3* DrbDy3b)/
			rbP3))/ (rbPz3b* rb) - (y2*(rbPy3bM1 + a/rbP2)*Dy3bDy3) /rb + 
			(cosb*y2*(((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)*
			Dy3bDy3)/ (rbPz3b* rb) + (y2*y3bMa* (rbPy3bM1 + a/rbP2)*Dy3bDy3*
			DrbDy3b)/ rbP2 - (cosb*y2*y3bMa* (((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b + (a*y3b)/rbP2)* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2));
			 
		dg_ic[2][1][0] +=
			2.*((y3bMa*(Dz3bDy1 + DrbDy1)* (sinb*(cosb - aQrb) - 
			(cosb*sinb*y2P2 - (a*rbTcosbPy3b* z1b)/rb)/ (rbPz3b* rb) + 
			((1. + (a*y3b)/rbP2)* z1b)/rb))/ rbPz3bP2 + pr3*(-((1. + aQrb)/
			rbPy3b) - (sinb*(Dz3bDy1 + DrbDy1))/ rbPz3b - ((cosb + aQrb)*
			(Dz3bDy1 + DrbDy1)* z1b)/ rbPz3bP2 + ((cosb + aQrb)* Dz1bDy1)/
			rbPz3b + (y1*(1. + aQrb)* DrbDy1)/ rbPy3bP2 + (a*y1*DrbDy1)/
			(rbPy3b* rbP2) - (a*z1b* DrbDy1)/ (rbPz3b* rbP2)) - 
			(y3bMa*(((cosb*sinb*y2P2 - (a*rbTcosbPy3b* z1b)/rb)* (Dz3bDy1 + 
			DrbDy1))/ (rbPz3bP2* rb) - (-((a*rbTcosbPy3b* Dz1bDy1)/ rb) + 
			(a*rbTcosbPy3b* z1b* DrbDy1)/ rbP2 - (a*cosb*z1b* DrbDy1)/ rb)/
			(rbPz3b* rb) + ((1. + (a*y3b)/rbP2)* Dz1bDy1)/ rb + 
			(a*sinb*DrbDy1)/ rbP2 + ((cosb*sinb*y2P2 - (a*rbTcosbPy3b*
			z1b)/rb)* DrbDy1)/ (rbPz3b* rbP2) - ((1. + (a*y3b)/rbP2)*
			z1b* DrbDy1)/ rbP2 - (2.*a*y3b*z1b* DrbDy1)/ rbP4))/ rbPz3b + 
			(y3bMa*(rbPy3bM1 + a/rbP2))/rb + (y1*y3bMa* (-(DrbDy1/ rbPy3bP2) - 
			(2.*a*DrbDy1)/ rbP3))/rb - (y1*y3bMa* (rbPy3bM1 + a/rbP2)*
			DrbDy1)/ rbP2); 
			
		dg_ic[2][1][1] +=
			2.*(pr3*((y1*(1. + aQrb)* DrbDy2)/ rbPy3bP2 - (sinb*DrbDy2)/
			rbPz3b + (a*y1*DrbDy2)/ (rbPy3b* rbP2) - ((cosb + aQrb)*z1b*
			DrbDy2)/ rbPz3bP2 - (a*z1b* DrbDy2)/ (rbPz3b* rbP2)) - 
			(y3bMa*(-((2.*cosb*sinb*y2 + (a*rbTcosbPy3b* z1b* DrbDy2)/ rbP2 - 
			(a*cosb*z1b* DrbDy2)/ rb)/ (rbPz3b* rb)) + (a*sinb*DrbDy2)/
			rbP2 + ((cosb*sinb*y2P2 - (a*rbTcosbPy3b* z1b)/rb)* DrbDy2)/
			(rbPz3b* rbP2) + ((cosb*sinb*y2P2 - (a*rbTcosbPy3b* z1b)/rb)*
			DrbDy2)/ (rbPz3bP2* rb) - ((1. + (a*y3b)/rbP2)* z1b* DrbDy2)/
			rbP2 - (2.*a*y3b*z1b* DrbDy2)/ rbP4))/ rbPz3b + (y1*y3bMa*
			(-(DrbDy2/ rbPy3bP2) - (2.*a*DrbDy2)/ rbP3))/rb + 
			(y3bMa*(sinb*(cosb - aQrb) - (cosb*sinb*y2P2 - (a*rbTcosbPy3b*
			z1b)/rb)/ (rbPz3b* rb) + ((1. + (a*y3b)/rbP2)* z1b)/rb)*
			DrbDy2)/ rbPz3bP2 - (y1*y3bMa* (rbPy3bM1 + a/rbP2)* DrbDy2)/ rbP2);
			
		dg_ic[2][1][2] +=
			2.*((y3bMa*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b)*
			(sinb*(cosb - aQrb) - (cosb*sinb*y2P2 - (a*rbTcosbPy3b*
			z1b)/rb)/ (rbPz3b* rb) + ((1. + (a*y3b)/rbP2)* z1b)/rb))/
			rbPz3bP2 + pr3*((y1*(1. + aQrb)* (Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPy3bP2 - (sinb*(Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			rbPz3b - ((cosb + aQrb)* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b)* z1b)/ rbPz3bP2 + ((cosb + aQrb)* Dz1bDy3b* Dy3bDy3)/
			rbPz3b + (a*y1*Dy3bDy3* DrbDy3b)/ (rbPy3b* rbP2) - (a*z1b*Dy3bDy3*
			DrbDy3b)/ (rbPz3b* rbP2)) - (y3bMa*(((cosb*sinb*y2P2 - 
			(a*rbTcosbPy3b* z1b)/rb)* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3* DrbDy3b))/
			(rbPz3bP2* rb) - (-((a* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b)*
			z1b)/rb) - (a*rbTcosbPy3b* Dz1bDy3b* Dy3bDy3)/rb + (a*rbTcosbPy3b*
			z1b*Dy3bDy3* DrbDy3b)/ rbP2)/ (rbPz3b* rb) + (((a*Dy3bDy3)/
			rbP2 - (2.*a*y3b*Dy3bDy3* DrbDy3b)/ rbP3)*z1b)/ rb + 
			((1. + (a*y3b)/rbP2)* Dz1bDy3b* Dy3bDy3)/rb + (a*sinb*Dy3bDy3*
			DrbDy3b)/ rbP2 + ((cosb*sinb*y2P2 - (a*rbTcosbPy3b* z1b)/rb)*
			Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) - ((1. + (a*y3b)/rbP2)*
			z1b*Dy3bDy3* DrbDy3b)/ rbP2))/ rbPz3b + (y1*y3bMa* (-((Dy3bDy3 + 
			Dy3bDy3* DrbDy3b)/ rbPy3bP2) - (2.*a*Dy3bDy3* DrbDy3b)/
			rbP3))/rb - ((sinb*(cosb - aQrb) - (cosb*sinb*y2P2 - 
			(a*rbTcosbPy3b* z1b)/rb)/ (rbPz3b* rb) + ((1. + (a*y3b)/rbP2)*
			z1b)/rb)* Dy3bDy3)/ rbPz3b + (y1*(rbPy3bM1 + a/rbP2)*Dy3bDy3)
			/rb - (y1*y3bMa* (rbPy3bM1 + a/rbP2)*Dy3bDy3* DrbDy3b)/ rbP2);
			 
		dg_ic[2][2][0] +=
			2.*((-2.*sinb*y2*(1. - pr)*(cosb + aQrb)* (Dz3bDy1 + DrbDy1))/
			rbPz3bP2 + pr2*(DrbDy1* DfbDrb + Dz1bDy1* DfbDz1b + DfbDy1) - 
			(sinb*y2*y3bMa* (Dz3bDy1 + DrbDy1)* (1. + ((cosb + aQrb)*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2))/ (rbPz3bP2* rb) + 
			(sinb*y2*y3bMa* (-(((cosb + aQrb)* rbTcosbPy3b* (Dz3bDy1 + 
			DrbDy1))/ rbPz3bP2) + (cosb*(cosb + aQrb)* DrbDy1)/ rbPz3b - 
			(a*rbTcosbPy3b* DrbDy1)/ (rbPz3b* rbP2) - (2.*a*y3b* DrbDy1)/
			rbP3))/ (rbPz3b* rb) - (2.*a*sinb*y2*(1. - pr)*DrbDy1)/
			(rbPz3b* rbP2) - (sinb*y2*y3bMa* (1. + ((cosb + aQrb)*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)* DrbDy1)/ (rbPz3b* rbP2));
			
		dg_ic[2][2][1] +=
			2.*((2.*sinb*(1. - pr)*(cosb + aQrb))/ rbPz3b + pr2*(DrbDy2*
			DfbDrb + DfbDy2) + (sinb*y3bMa* (1. + ((cosb + aQrb)*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2))/ (rbPz3b* rb) + 
			(sinb*y2*y3bMa* ((cosb*(cosb + aQrb)* DrbDy2)/ rbPz3b - 
			((cosb + aQrb)* rbTcosbPy3b* DrbDy2)/ rbPz3bP2 - (a*rbTcosbPy3b*
			DrbDy2)/ (rbPz3b* rbP2) - (2.*a*y3b* DrbDy2)/ rbP3))/ (rbPz3b*
			rb) - (2.*sinb*y2*(1. - pr)*(cosb + aQrb)* DrbDy2)/ rbPz3bP2 - 
			(2.*a*sinb*y2*(1. - pr)*DrbDy2)/ (rbPz3b* rbP2) - (sinb*y2*y3bMa*
			(1. + ((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)*
			DrbDy2)/ (rbPz3b* rbP2) - (sinb*y2*y3bMa* (1. + ((cosb + aQrb)*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)* DrbDy2)/ (rbPz3bP2* rb)); 
			
		dg_ic[2][2][2] +=
			2.*((-2.*sinb*y2*(1. - pr)*(cosb + aQrb)* (Dz3bDy3b* Dy3bDy3 + 
			Dy3bDy3* DrbDy3b))/ rbPz3bP2 + pr2*(Dy3bDy3* DrbDy3b* DfbDrb + 
			Dz1bDy3b* Dy3bDy3* DfbDz1b) - (sinb*y2*y3bMa* (Dz3bDy3b*
			Dy3bDy3 + Dy3bDy3* DrbDy3b)* (1. + ((cosb + aQrb)* rbTcosbPy3b)/
			rbPz3b + (a*y3b)/rbP2))/ (rbPz3bP2* rb) + (sinb*y2*y3bMa*
			(((cosb + aQrb)* (Dy3bDy3 + cosb*Dy3bDy3* DrbDy3b))/ rbPz3b - 
			((cosb + aQrb)* rbTcosbPy3b* (Dz3bDy3b* Dy3bDy3 + Dy3bDy3*
			DrbDy3b))/ rbPz3bP2 + (a*Dy3bDy3)/rbP2 - (a*rbTcosbPy3b*
			Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) - (2.*a*y3b*Dy3bDy3*
			DrbDy3b)/ rbP3))/ (rbPz3b* rb) + (sinb*y2*(1. + ((cosb + aQrb)*
			rbTcosbPy3b)/ rbPz3b + (a*y3b)/rbP2)* Dy3bDy3)/ (rbPz3b* rb) - 
			(2.*a*sinb*y2*(1. - pr)*Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2) - 
			(sinb*y2*y3bMa* (1. + ((cosb + aQrb)* rbTcosbPy3b)/ rbPz3b + 
			(a*y3b)/rbP2)* Dy3bDy3* DrbDy3b)/ (rbPz3b* rbP2)); 

	}

	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			for (k=0; k < 3; k++) {
				dg_ic[i][j][k] /= 8. * PI * (1. - pr);
			}
		}
	}
	/* Calculate the strain influence coeffs from the displacement
	 * gradient influence coeffs
	--------------------------------------------------------------*/
	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			for (k=0; k < 3; k++) {
				dummy[i][j][k] =
				strain_ic[i][j][k] = 0.5 * (dg_ic[i][j][k] + dg_ic[i][k][j]);
			}
		}
	}
}


/************************ Function: reduce_angle ***************************
* Forces angle to between +/- PI by adding or subtracting multiples of
* 2*PI as necessary.  Returns the "reduced" value.
*
* In:	angle	- the angle to be "reduced"
***************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
double   reduce_angle(double angle)
#else
double reduce_angle(angle)
double	angle;
#endif
{
	double	two_pi	= 2.0*PI;

	while (fabs(angle) > PI) {
		if (angle > 0) {
			angle -= two_pi;
		} else {
			angle += two_pi;
		}
	}

	return(angle);
}
