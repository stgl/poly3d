/*================================================== 
SET TABSTOPS AT EVERY FOUR SPACES FOR PROPER DISPLAY
====================================================*/

/***************************************************************************** 
* FILE: getoptPoly3D.c
* DATE: June, 1993 
* BY:   Andrew L. Thomas 
*
* Header file for getopt.c
*****************************************************************************/

#ifndef getoptPoly3D_h_


/***************************** Includes/Defines ******************************/
#define getoptPoly3D_h_

#define FILE_ARG 0
#define FLAG 1
#define FLAG_AND_FILE 2
#define NO_MORE_ARGS -1
#define NO_SUCH_ARG -2


/********************************* Structures ********************************/
struct arg_rec {
	char flag_name;
	char arg_type;
	struct arg_rec *next_rec;
	char *file_arg;
};


/******************************* External Variables **************************/
extern char *getopt_arg_E;


/****************************** Function Declarations ************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
int getoptPoly3D(char *arg_string, int argc, char *argv[]);
#else
int getoptPoly3D();
#endif

#endif
