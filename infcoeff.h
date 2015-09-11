/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: infcoeff.h
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* Header file for infcoeff.c
*****************************************************************************/

/********************************* Constants ********************************/
#define COMNINOU_TINY 1.0e-14

/******************************* Function Declarations **********************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void	comninou_displ_ics(double y[3], double a, double beta,
		double pr, int half_space, double displ_ic[3][3]);
void	comninou_strain_ics(double y[3], double a, double beta,
		double pr, int half_space, double strain_ic[3][3][3]);
#else
void	comninou_displ_ics();
void	comninou_strain_ics();
#endif
