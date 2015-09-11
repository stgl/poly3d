/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: getwords.c
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* A function for reading input lines and setting an array of pointers to
* each word on the line.
*****************************************************************************/


/***************************** Includes/Defines *****************************/
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <ctype.h>
#include <string.h>
#include "getwords.h"

#define TRUE				1
#define FALSE				0


/**************************** function: getwords ***********************
*
* Reads a line of input into the character array line from the file ifp
* using fgets.  Copies line to the static array wline.  Then sets the 
* character pointers word[0]..word[n-1] to point to the first letters
* of the 1st..nth words in line.  Changes the space following each word 
* in line to '\0', so that each word[i] will be a null terminated, 
* one-word string.
*
* Returns the number of words on the line, or one of the following error
* codes:
*				GW_EOF_ERR		-	End of File
*				GW_MALLOC_ERR	-	Malloc() error
*				GW_MAXWORDS_ERR	-	Max number of words exceeded error
*
* In:	ifp				- input file pointer
*		maxline			- maximum number of characters on a line
*		maxword			- maximum number of words on a line
*		continue_char	- line continuation character
* Out:	line			- input line
*		word			- array of words on input line
************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int getwords(FILE *ifp, char *line, int maxline, char *word[],int maxwords, char continue_char)
#else
int getwords(ifp, line, maxline, word, maxwords, continue_char)
FILE	*ifp;
char	*line;
int 	maxline;
char	*word[];
int		maxwords;
char	continue_char;
#endif

{
	int		i;
	int		j			= 0;					/* index for word[]			*/
	int		inword		= FALSE;
	int		inquotes	= FALSE;
	char	*start_at;

	static char	*wline = NULL;
	static int	oldmaxline;

	/* Free/allocate memory for wline if necessary
	----------------------------------------------*/
	if ((wline == NULL) || (maxline > oldmaxline)) {
		if (wline != NULL)
			free(wline);
		if ((wline = (char *) malloc((size_t) maxline)) == NULL)
			return(GW_MALLOC_ERR);
		oldmaxline = maxline;
	}

	/* Read lines from ifp until a line not ending in continue_char is
	   found.  Concatonate these lines together into a single getwords() line
	-------------------------------------------------------------------------*/
	start_at = line;
	while (start_at != NULL) {
		if (start_at != line) {
			*start_at = ' ';
			start_at++;
		}
		if (fgets(start_at,maxline-(start_at-line),ifp) == NULL)
			return(GW_EOF_ERR);
		start_at = strchr(line,continue_char);
	}

	strcpy(wline,line);

	/* Step through line and set word pointers
	------------------------------------------*/
	for (i = 0; wline[i] != '\0'; i++) {
		if (wline[i]=='"') {
			if (!inquotes) {
				if (inword) {
					if (j >= maxwords-1) return (GW_MAXWORDS_ERR);
					word[j++] = &wline[i+1];
				} /*if*/
				else
					if (j >= maxwords-1) return (GW_MAXWORDS_ERR);
					word[j++] = &wline[i+1];
				inquotes = TRUE;
			} /*if*/
			else {
				wline[i] = '\0';
				inquotes = FALSE;
			} /*else*/
		} /*if*/
		else if (!inquotes && !inword && !isspace(wline[i])) {
			if (j >= maxwords-1) return (GW_MAXWORDS_ERR);
			word[j++] = &wline[i];
			inword = TRUE;
		} /*if*/
		else if (inword && isspace(wline[i])) {
			wline[i] = '\0';
			inword = FALSE;	
		} /*else if*/
	} /*for*/

	return j;
}
