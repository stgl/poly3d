/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: safetan.c
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* The standard tan(alpha), atan(y/x), and atan2(y,x) functions will return 
* NaN if...
*			(1) alpha is very close to N*PI/2, N = +/- 1,2,3...
*			(2) x = 0
*			(3) x << y (so that y/x -> NaN)
*
* The following "safe" tangent and arctangent functions check for and
* correct these problems
*****************************************************************************/

#include <math.h>
#include "safetan.h"
#include "pi.h"

/*************************** function: safe_tan ***************************
* Computes tan(alpha) avoiding the singularity at alpha = N*PI/2,
* N = +/- 1,2,3....
**************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
safetan_t safe_tan(safetan_t alpha)
#else
safetan_t safe_tan(alpha)
safetan_t alpha;
#endif

{

	safetan_t	half_pi = PI/2.0;

	/* Force alpha between +/- pi/2 (Note: If alpha is exactly
	   +/- pi, this algorithm only changes it's sign, but this is
	   inconsequential for the purposes of this function).
	---------------------------------------------------------------*/
	if (fabs(alpha) > half_pi) {
		alpha += (alpha > 0) ?	  PI * (((int) alpha/PI) + 1) :
								-(PI * (((int) alpha/PI) + 1));
	}

	/* Don't let alpha be exactly +- pi/2 or tan() will blow up
	-------------------------------------------------------------*/
	if ((half_pi - fabs(alpha)) < SAFE_TAN_LIMIT) {
        alpha = (alpha > 0) ?	  half_pi - SAFE_TAN_LIMIT :
								-(half_pi - SAFE_TAN_LIMIT);
	}

	return(tan(alpha));
}

/*************************** function: safe_atan **************************
* Computes atan(y/x) avoiding the singularities at x = 0 and x << y.
**************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
safetan_t safe_atan(safetan_t	y,
					safetan_t	x)
#else
safetan_t safe_atan(y, x)
safetan_t	y;
safetan_t	x;
#endif

{
	/* Avoid huge y/x ratios
	------------------------*/
	if (fabs(x*SAFE_ATAN_LIMIT) < fabs(y))
		return((y >= 0.0) ? PI/2.0 : -(PI/2.0)); 

	/* Avoid atan2 domain error for x=y=0
	-------------------------------------*/
	if (y == 0.0 && x == 0.0)
		return(0.0);

	return(atan(y/x));
}

/************************** function: safe_atan2 ****************************
* Computes atan2(y,x) avoiding the singularities at x = 0 and x << y.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
safetan_t safe_atan2(	safetan_t	y,
						safetan_t	x)
#else
safetan_t safe_atan2(y, x)
safetan_t y;
safetan_t x;
#endif

{

	/* Avoid huge y/x ratios
	------------------------*/
	if (fabs(x*SAFE_ATAN_LIMIT) < fabs(y))
		return((y > 0) ? PI/2.0 : -(PI/2.0)); 

	/* Avoid atan2 domain error for x=y=0
	-------------------------------------*/
	if (y == 0.0 && x == 0.0)
		return((x >= 0.0) ? 0.0 : PI);

	return(atan2(y,x));

}
