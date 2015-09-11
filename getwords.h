/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: getwords.h
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* Header file for getwords.c
*****************************************************************************/


/****************************** Includes/Defines ****************************/
#define GW_EOF_ERR			-1
#define GW_MALLOC_ERR		-2
#define GW_MAXWORDS_ERR		-3


/**************************** Function Declarations *************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int getwords(FILE *ifp, char *line, int maxline, char *word[], int maxwords,
	char continue_char);
#else
int getwords();
#endif
