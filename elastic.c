/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: elastic.c
* DATE: June, 1993
* BY  : Andrew L. Thomas
*
* Functions that perform simple calculations from linear elasticity
******************************************************************************/

#include <math.h>
#include "elastic.h"
#include "matrix.h"


/************************* Function: calc_elas_consts ************************
* Calculates the values of all elastic constants given any two.  Returns
* EC_TOO_MANY or EC_TOO_FEW if the number of already defined constants is 
* not exactly two.
*
* In/Out:	sm	- Shear Modulus
*			pr	- Poisson's Ratio
*			ym	- Young's Modulus
*			bm	- Bulk Modulus
*			ll	- Lame Lambda
******************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int calc_elas_consts(elas_const_t *sm, elas_const_t *pr, elas_const_t *ym,elas_const_t *bm, elas_const_t *ll)
#else
int calc_elas_consts(sm, pr, ym, bm, ll)
elas_const_t	*sm;
elas_const_t	*pr;
elas_const_t	*ym;
elas_const_t	*bm;
elas_const_t	*ll;
#endif
{
	int num_defined;						/* # of consts already defined	*/
	elas_const_t temp1, temp2;				/* temporary variables			*/

	/* Count the number of defined constants
	----------------------------------------*/
	num_defined = 0;
	if EC_DEFINED(*sm) num_defined++;
	if EC_DEFINED(*pr) num_defined++;
	if EC_DEFINED(*ym) num_defined++;
	if EC_DEFINED(*bm) num_defined++;
	if EC_DEFINED(*ll) num_defined++;

	/* Abort if number defined is > or < 2
	--------------------------------------*/
	if (num_defined < 2) {
		return(EC_TOO_FEW);
	} else if (num_defined > 2) {
		return(EC_TOO_MANY);
	}

	/* Calculate the shear modulus (if not already defined)
	-------------------------------------------------------*/
	if (!EC_DEFINED(*sm)) {
		if EC_DEFINED(*ym) {
			if EC_DEFINED(*pr)
				(*sm) = (*ym)/(2.*(1.+(*pr)));
			else if EC_DEFINED(*bm)
				(*sm) = 3.*(*ym)*(*bm)/(9*(*bm)-(*ym));
			else if EC_DEFINED(*ll) {
				temp1 = (3.*(*ll)-(*ym));
				temp2 = sqrt(temp1*temp1 + 8.*(*ll)*(*ym));
				(*sm) = (-temp1 + temp2)/4.;
			}
		}
		else if EC_DEFINED(*pr) {
			if EC_DEFINED(*bm)
				(*sm) = 3./2. * (1.-2.*(*pr))/(1.+(*pr)) * (*bm);
			else if EC_DEFINED(*ll) 
				(*sm) = (1.-2.*(*pr))/(2.*(*pr)) * (*ll);
		}
		else {
			(*sm) = (1.-2.*(*pr))*(*ll) / (2.*(*pr));
		}
	}

	/* Calculate Lame's lambda (if not already defined)
	---------------------------------------------------*/
	if (!EC_DEFINED(*ll)) {
		if EC_DEFINED(*pr)
			(*ll) = (2.*(*pr))/(1.-2.*(*pr)) * (*sm);
		else if EC_DEFINED(*ym)
			(*ll) = (2.*(*sm)-(*ym))/((*ym)-3.*(*sm))
				* (*sm);
		else
			(*ll) = (*bm) - 2./3.*(*sm);
	}

	/* Calculate Young's modulus, Poisson's ratio, and the bulk modulus
	   (if not already defined) from the shear modulus and Lame's lambda
	--------------------------------------------------------------------*/
	if (!EC_DEFINED(*ym)) 
		(*ym) = (3.*(*ll)+2.*(*sm))/((*ll)+(*sm))
			* (*sm);
	if (!EC_DEFINED(*pr))
		(*pr) = (*ll) / (2. * ((*ll) + (*sm)));
	if (!EC_DEFINED(*bm))
		(*bm) = (*ll) + 2./3.*(*sm);

	return(0);
}

/************************ Function: strain_to_stress ************************
* Uses Hooke's law to convert the given strain tensor into a stress tensor
*
* In:	strain		- the strain tensor to convert
*		shear_mod	- shear modulus
*		lame_lambda - Lame lambda
* Out:	stress		- the calculated stress tensor
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void strain_to_stress(double strain[3][3], double shear_mod,double lame_lambda, double stress[3][3])
#else
void strain_to_stress(strain, shear_mod, lame_lambda, stress)
double	strain[3][3];
double	shear_mod;
double	lame_lambda;
double	stress[3][3];
#endif

{
	int		i,j;
	double	strain_copy[3][3];
	double	vol_strain = strain[0][0] + strain[1][1] + strain[2][2];

	copy_matrix(strain,strain_copy);
	for (i=0; i < 3; i++) {
		stress[i][i] = lame_lambda * vol_strain + 
			2.* shear_mod * strain_copy[i][i];
		for (j=0; j < 3; j++) {
			if (i != j) {
				stress[i][j] = 2. * shear_mod * strain_copy[i][j];
			}
		}
	}
}

/************************ Function: stress_to_strain ************************
* Uses Hooke's law to convert the given strain tensor into a stress tensor
*
* In:	stress		- the stress tensor to convert
*		youngs_mod	- Young's modulus
*		psn_ratio	- Poisson's ratio
* Out:	strain		- the calculated strain tensor
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void stress_to_strain(double stress[3][3], double youngs_mod,double psn_ratio, double strain[3][3])
#else
void stress_to_strain(stress, youngs_mod, psn_ratio, strain)
double	stress[3][3];
double	youngs_mod;
double	psn_ratio;
double	strain[3][3];
#endif

{
	int		i,j,k;
	double	stress_copy[3][3];
	double	shear_mod = youngs_mod/(2.*(1.+psn_ratio));

	copy_matrix(stress,stress_copy);
	for (i=0; i < 3; i++) {
		j = i+1;
		if (j > 2) j -= 3;
		k = i+2;
		if (k > 2) k -= 3;
		strain[i][i] = (stress_copy[i][i] - psn_ratio * (stress_copy[j][j]
			+ stress_copy[k][k]))/youngs_mod;
		for (j=0; j < 3; j++) {
			if (i != j) {
				strain[i][j] = stress_copy[i][j] / (2. * shear_mod);
			}
		}
	}
}
