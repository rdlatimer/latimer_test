#ifdef matrixreading/*  FUNCTIONS TO READ MATRICES & ARRAYS AND RETURN INFORMATION ABOUT THOSE MATRICES */	#include <stdlib.h>	#include <stdio.h>	#include <time.h>	#include <math.h>	#include <string.h>	int *numberstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	int *mincharstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	int *maxcharstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	int *numb_states(long **omat, int *nstates, int notu, int nchars, int UNKNOWN, int INAP);	int no_states(long **omat, int ch, int notu, int UNKNOWN, int INAP);	int *autapomorphies(long **omat, int *nstates, int notu, int nchars, int UNKNOWN, int INAP);	int *tallyautaps(long **omat, int *nstates, int notu, int nchars);	int autapo_char(long **omat, int notu, int ch, int unknown, int inap);	int *maxchanges(long **matrix, int nchars, int notu, int *nstates, int *ctype, int INAP, int UNKNOWN);	int *unknownstates(long **omat, int notu, int nchars, int unknown, int inap);	long **scalematrix(long **matrix, int notu, int nchars, int unknown, int inap, int scale);	int *applicabletaxa(long **mat, int nchars, int notu, int inap);	int inapplicables(long **omat, int notu, int nchars, int inap);	int *inaplist(long **omat, int notu, int nchars, int inap, int depend);	int *findindependents(long **cmat, int notu, int inap, int unknown, int depend, int *dependents);	int findinlmatrix(long **omat, int rows, int col, int find, int start, int end);	int maxentryivector(int *v, int n);	int countentryivector(int *v, int n);	int maxinlmatrix(long **mat, int r, int c);	int maxinclmatrix(long **cmat, int notu, int nchar, int unknown, int inap);	int compareivector(int *v1, int *v2, int n);	int maxstringlength(char **strings, int n, char eol);	int *stringlengths(char **strings, int n, char eol);	int	colminclmatrix(long **mat, int notu, int c, int unknown, int inap);	int	colmaxclmatrix(long **mat, int notu, int c, int unknown, int inap);	int minclmatrix(long **mat, int notu, int c, int unknown, int inap);	int maxclmatrix(long **mat, int notu, int c, int unknown, int inap);	int countuniqivector(int *v, int n);	int countuniqlvector(long *v, int n);	int countstatescharvector(long *chvector, int notu, int UNKNOWN, int INAP);	long **staterichness(long **mat, int notu, int nchar, int inap, int unknown);	int readicube(int **mat, int i, int j, int k, int mxk);	long readlcube(long **mat, int i, int j, int k, int mxk); 	unsigned long readucube(unsigned long **mat, int i, int j, int k, int mxk); 	float readfcube(float **mat, int i, int j, int k, int mxk); 	double readdcube(double **mat, int i, int j, int k, int mxk); 	void **assignicube(int **mat, int i, int j, int k, int mxk, int x); 	void **assignlcube(long **mat, int i, int j, int k, int mxk, long x); 	void **assignucube(unsigned long **mat, int i, int j, int k, int mxk, unsigned long x); 	void **assignfcube(float **mat, int i, int j, int k, int mxk, float x); 	void **assigndcube(double **mat, int i, int j, int k, int mxk, double x);	long **charstateunkncombos (int *nstates, int *chuns, int nchars);	void *extractlongcolumn(long **mat, long *column, int i, int rows);	long **charactersbystates(int *states, int nchars, int maxst);	unsigned long *categorizeulmatrix(unsigned long **Data, int Cat, int Mem, int N);	int *countunscoredtaxa(long **chmat, int notu, int nchars, int UNKNOWN, int INAP);	int tally_unscored_for_taxon(long **chmat, int sp, int nchars, int UNKNOWN, int INAP);	#else	extern int *numberstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	extern int *mincharstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	extern int *maxcharstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	extern int *numb_states(long **omat, int *nstates, int notu, int nchars, int UNKNOWN, int INAP);	extern int no_states(long **omat, int ch, int notu, int UNKNOWN, int INAP);	extern int *autapomorphies(long **omat, int *nstates, int notu, int nchars, int UNKNOWN, int INAP);	extern int *tallyautaps(long **omat, int *nstates, int notu, int nchars);	extern int autapo_char(long **omat, int notu, int ch, int UNKNOWN, int INAP);	extern int *maxchanges(long **matrix, int nchars, int notu, int *nstates, int *ctype, int INAP, int UNKNOWN);	extern int *unknownstates(long **omat, int notu, int nchars, int UNKNOWN, int INAP);	extern long **scalematrix(long **matrix, int notu, int nchars, int UNKNOWN, int INAP, int scale);	extern int *applicabletaxa(long **mat, int nchars, int notu, int inap);	extern int inapplicables(long **omat, int notu, int nchars, int inap);	extern int *inaplist(long **omat, int notu, int nchars, int inap, int depend);	extern int *findindependents(long **cmat, int notu, int inap, int unknown, int depend, int *dependents);	extern int findinlmatrix(long **omat, int rows, int col, int find, int start, int end);	extern int maxentryivector(int *v, int n);	extern int countentryivector(int *v, int n);	extern int maxinclmatrix(long **cmat, int notu, int nchar, int unknown, int inap);	extern int compareivector(int *v1, int *v2, int n);	extern int maxstringlength(char **strings, int n, char eol);	extern int *stringlengths(char **strings, int n, char eol);	extern int colminclmatrix(long **mat, int notu, int c, int unknown, int inap);	extern int colmaxclmatrix(long **mat, int notu, int c, int unknown, int inap);	extern int minclmatrix(long **mat, int notu, int c, int unknown, int inap);	extern int maxclmatrix(long **mat, int notu, int c, int unknown, int inap);	extern int countuniqivector(int *v, int n);	extern int countuniqlvector(long *v, int n);	extern int countstatescharvector(long *chvector, int notu, int UNKNOWN, int INAP);	extern long **staterichness(long **mat, int notu, int nchar, int inap, int unknown);	extern int readicube(int **mat, int i, int j, int k, int mxk);	extern long readlcube(long **mat, int i, int j, int k, int mxk); 	extern unsigned long readucube(unsigned long **mat, int i, int j, int k, int mxk); 	extern float readfcube(float **mat, int i, int j, int k, int mxk); 	extern double readdcube(double **mat, int i, int j, int k, int mxk); 	extern void **assignicube(int **mat, int i, int j, int k, int mxk, int x); 	extern void **assignlcube(long **mat, int i, int j, int k, int mxk, long x); 	extern void **assignucube(unsigned long **mat, int i, int j, int k, int mxk, unsigned long x); 	extern void **assignfcube(float **mat, int i, int j, int k, int mxk, float x); 	extern void **assigndcube(double **mat, int i, int j, int k, int mxk, double x); 	extern long **charstateunkncombos (int *nstates, int *chuns, int nchars);	extern void *extractlongcolumn(long **mat, long *column, int i, int rows);	extern long **charactersbystates(int *states, int nchars, int maxst);	extern unsigned long *categorizeulmatrix(unsigned long **Data, int Cat, int Mem, int N);	extern int *countunscoredtaxa(long **chmat, int notu, int nchars, int UNKNOWN, int INAP);	extern int tally_unscored_for_taxon(long **chmat, int sp, int nchars, int UNKNOWN, int INAP);#endif