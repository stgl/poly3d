/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/****************************************************************************
* FILE: elastic.h
* DATE: June, 1993
* BY  : Andrew L. Thomas
*
* Header file for elastic.c
*****************************************************************************/

/******************************* Defines/Macros *****************************/
#define		EC_DEFINED(A)	((A) >= 0)	/* Macro to test whether elastic	*/
										/*  constant is defined				*/
#define		EC_UNDEF_VAL	-1.0		/* Value for undefined elas consts	*/
#define		EC_TOO_FEW		-1			/* Return if < 2 consts predefined	*/
#define		EC_TOO_MANY		-2			/* Return if > 2 consts predefined	*/

typedef	double	elas_const_t;			/* Type definition for elas consts	*/

/************************* Function Declarations ****************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int		calc_elas_consts(elas_const_t *sm, elas_const_t *pr, elas_const_t *ym,
		elas_const_t *bm, elas_const_t *ll);
void	strain_to_stress(double strain[3][3], double shear_mod,
		double lame_lambda, double stress[3][3]);
void	stress_to_strain(double stress[3][3], double youngs_mod,
		double psn_ratio, double strain[3][3]);
#else
int		calc_elas_consts();
void	strain_to_stress();
void	stress_to_strain();
#endif
