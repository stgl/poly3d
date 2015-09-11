/*==================================================
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/*****************************************************************************
* FILE: getoptPoly3D.c
* DATE: June, 1993
* BY:   Andrew L. Thomas
*
* Functions for parsing command line arguments, checking them for
* compatibility with allowed program options.
*****************************************************************************/


/***************************** Includes/Defines *****************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "getoptPoly3D.h"


/**************************** External Variables ***************************/
char    *getopt_arg_E;


/*************************** Function Declarations *************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
struct arg_rec *read_args(int argc, char *argv[], char *arg_string);
void add_arg(struct arg_rec **top, int i, char arg_type,
                                char *argv[]);
#else
struct arg_rec *read_args();
void add_arg();
#endif


/**************************** Function: getoptPoly3D ******************************
*
* Creates a linked list of arguments to the program based on the arg_string
* passed in.  If an argument has a ':' after it the next argv argument belongs
* with that flag.  It tags the structures as 4 types: FILE_ARG, meaning that
* there is no flag that goes with it; FLAG, a -x argument, with no file that
* goes with it; a FLAG_AND_FILE, which has a -x argument followed by a file
* argument; and finally an error argument NO_SUCH_ARG.  The structure of an
* argument list is defined in getoptPoly3D.h as follows:
* 
*		struct arg_rec{
*			char flag_name;
*			char arg_type;
*			struct arg_rec *next_rec;
*			char *file_arg;
*		};
*
* Adapted from a similar function by Dan Kolkowitz (Nov, 1990) written as
* part of a solution set for CS040, Stanford University.
*
* In:	arg_string	- the argument string (see above)
*		argc		- number of command line arguments
*		argv		- array of commane line arguments
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int               getoptPoly3D(char *arg_string, int argc, char *argv[])
#else
int getoptPoly3D(arg_string, argc, argv)
char *arg_string;
int argc;
char *argv[];
#endif
              
{
	int i=1;
	static int args_parsed=0;
	static struct arg_rec *next_arg, *current_arg;

	/* The first time getoptPoly3D is called, call read_args() to set up
	   the linked list of arguments. Also set the next_arg pointer.
	---------------------------------------------------------------*/
	if (!args_parsed)
		if ((next_arg = read_args(argc,argv,arg_string)) 
			== (struct arg_rec *) NO_SUCH_ARG)
			return (NO_SUCH_ARG);
	args_parsed++;

	if (!next_arg)
		return (NO_MORE_ARGS); 

	current_arg = next_arg;
	next_arg = next_arg->next_rec;

	switch (current_arg->arg_type) {
		case FILE_ARG:
			getopt_arg_E = current_arg->file_arg;
			return (FILE_ARG);
		case FLAG:
			getopt_arg_E = NULL;
			return (current_arg->flag_name);
		case FLAG_AND_FILE:
			getopt_arg_E = current_arg->file_arg;
			return (current_arg->flag_name);
		case NO_SUCH_ARG:
			return (NO_SUCH_ARG);
	} /*switch*/
}             


/************************** Function: read_args ************************
* Parses each argument and determines what type to add.
* Calls add_arg() to allocate the actual structure and set the
* actual field types.
* 
* In:	arg_string	- the argument string (see above)
*		argc		- number of command line arguments
*		argv		- array of commane line arguments
************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
struct arg_rec*   read_args(int argc, char *argv[], char *arg_string)
#else
struct arg_rec *read_args(argc, argv, arg_string)
int argc;
char *argv[];
char *arg_string;
#endif

{
	int i=1;
	char *c_ptr, newarg;
	static struct arg_rec *top;
	
	while (i < argc) {
		if (argv[i][0] == '-') {
			newarg = argv[i][1];
			c_ptr = (char *) strchr(arg_string, newarg);
			/* c_ptr = (char *) index(arg_string, newarg);  */
			if (!c_ptr)
				return ((struct arg_rec *) NO_SUCH_ARG);
			if (*(c_ptr+1) == ':') {
				add_arg(&top, i, FLAG_AND_FILE, argv);
				i+=2;
			} else add_arg(&top, i++, FLAG, argv);
		} else add_arg(&top, i++, FILE_ARG, argv);
	}
	return (top);
}

/**************************** Function: add_arg **************************
* Allocates an actual node.
* Mallocs the space for the structure.
* Sets a pointer to the argument if it is required.
*
* In:		i			- command line argument number
*			arg_type	- argument type
*			argv		- array of command line arguments
* In/Out:	top 		- first member of argument list
**************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void              add_arg(struct arg_rec **top, int i, char arg_type, char *argv[])
#else
void add_arg(top, i, arg_type, argv)
struct arg_rec **top;
int i;
char arg_type;
char *argv[];
#endif

{
	static struct arg_rec *last_arg;
	struct arg_rec *test_var;
	int *test_var2;

	/* The first time add_arg is called, top points nowhere, and the
	   the first node in the linked list of arguments must be created.
	   Allocate the memory needed for the new argument record and
	   reset the last_arg pointer.
	-----------------------------------------------------------------*/
	if (!*top) {
		*top = (struct arg_rec *) malloc(sizeof(struct arg_rec));
		test_var = (struct arg_rec *) malloc(sizeof(struct arg_rec));
		test_var2 = (int *) malloc(sizeof(int));
		last_arg = *top;
	} /*if*/
	else {
		last_arg->next_rec = (struct arg_rec *) malloc(sizeof(struct arg_rec));
		last_arg = last_arg->next_rec;
	} /*else*/

	/* Put the appropriate info into the new (last) argument record
	---------------------------------------------------------------*/
	last_arg->next_rec = NULL;
	last_arg->arg_type = arg_type;

	switch (arg_type) {
		case FILE_ARG:
			last_arg->file_arg = argv[i];
			break;
		case FLAG_AND_FILE:
			last_arg->file_arg = argv[i+1];
		case FLAG:
			last_arg->flag_name = argv[i][1];
			break;
	} /*switch*/
}

