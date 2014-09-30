//#define minmax #include "minmax.h"#include <stdlib.h>#include <stdio.h>#include <time.h>#include <math.h>#include <string.h>/* find the maximum of two integers */int imax(int i, int j){int	k;if (i>j)	k=i;else		k=j;return k;}/* find the minimum of two integers */int imin(int i, int j){int	k;if (i<j)	k=i;else		k=j;return k;}/* find the minimum of two long integers */long lmin(long i, long j){int	k;if (i<j)	k=i;else		k=j;return k;}/* find the maximum of two real numbers */double dmax(double i, double j){double	k;if (i>j)	k=i;else		k=j;return k;}/* find the minimum of two integers */double dmin(double i, double j){double	k;if (i<j)	k=i;else		k=j;return k;}/* find the maximum of two unsigned long */int ulmax(unsigned long i, unsigned long j){int	k;if (i>j)	k=i;else		k=j;return k;}/* find the minimum of two unsigned long */int ulmin(unsigned long i, unsigned long j){int	k;if (i<j)	k=i;else		k=j;return k;}/* maxiarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/int maxiarray(int *i, int N){int	a,b;b=-1000000;for (a=0; a<N; ++a)	if (i[a]>b)	b=i[a];return b;}/* miniarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/int miniarray(int *i, int N){int j, min;min=RAND_MAX;for (j=0; j<N; ++j)	if (i[j]<min)	min=i[j];return min;}/* maxiarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/int maxlarray(long *i, int N){int	a,b;b=-1000000;for (a=0; a<N; ++a)	if (i[a]>b)	b=i[a];return b;}/* miniarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/int minlarray(long *i, int N){int j, min;min=RAND_MAX;for (j=0; j<N; ++j)	if (i[j]<min)	min=i[j];return min;}/* maxiarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/int maxularray(unsigned long *i, int N){int	a,b;b=0;for (a=0; a<N; ++a)	if (i[a]>b)	b=i[a];return b;}/* minularray - finds the maximum number in an unsigned long array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/int minularray(unsigned long *i, int N){int j, min;min=RAND_MAX;for (j=0; j<N; ++j)	if (i[j]<min)	min=i[j];return min;}/* maxdarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/double maxdarray(double *i, int N){int		a;double	max;max=-1000000;for (a=0; a<N; ++a)	if (i[a]>max)	max=i[a];return max;}/* mindarray - finds the maximum number in an array/*/* Requires: i - an array to be examined;/*			 N - the length of the array/*/* Returns: b (the maximum number)/**************************************************************************/double mindarray(double *i, int N){int j;double	min;min=RAND_MAX;for (j=0; j<N; ++j)	if (i[j]<min)	min=i[j];return min;}/* maxulmatrixcol - finds the maximum number in a particular column of an unsigned long matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/int maxulmatrixcol(unsigned long **i, int r, int c){int	j, max=0;for (j=0; j<r; ++j)	{	if (i[j][c]>max)	{		max=i[j][c];		}	}return max; }/* minulmatrixcol - finds the minimum number in a particular column of an unsigned long matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/int minulmatrixcol(unsigned long **i, int r, int c){int	j, min=RAND_MAX;for (j=0; j<r; ++j)	if (i[j][c]<min)	min=i[j][c];return min;}/* maxlmatrixcol - finds the maximum number in a particular column of a long matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/int maxlmatrixcol(long **i, int r, int c){int	j, max=0;for (j=0; j<r; ++j)	{	if (i[j][c]>max)	{		max=i[j][c];		}	}return max; }/* minlmatrixcol - finds the minimum number in a particular column of a long matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/int minlmatrixcol(long **i, int r, int c){int	j, min=RAND_MAX;for (j=0; j<r; ++j)	if (i[j][c]<min)	min=i[j][c];return min;}/* maximatrixcol - finds the maximum number in a particular column of an integer matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/int maximatrixcol(int **i, int r, int c){int	j, max=0;for (j=0; j<r; ++j)	{	if (i[j][c]>max)	{		max=i[j][c];		}	}return max; }/* minimatrixcol - finds the minimum number in a particular column of an integer matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/int minimatrixcol(int **i, int r, int c){int	j, min=RAND_MAX;for (j=0; j<r; ++j)	if (i[j][c]<min)	min=i[j][c];return min;}/* maxdmatrixcol - finds the maximum number in a particular column of a long matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/double maxdmatrixcol(double **i, int r, int c){int	j;double max=0.0f;for (j=0; j<r; ++j)	{	if (i[j][c]>max)	{		max=i[j][c];		}	}return max; }/* mindmatrixcol - finds the minimum number in a particular column of a long matrix/*/* Requires: i - an matrix to be examined;/*			 r - the number of rows;/*			 c - the particular column being examined;/* Returns: max (the maximum number)/**************************************************************************/double mindmatrixcol(double **i, int r, int c){int	j;double min=RAND_MAX;for (j=0; j<r; ++j)	if (i[j][c]<min)	min=i[j][c];return min;}int maxulmatrix(unsigned long **i, int r, int c){int	a, j, k;int max=(-1*RAND_MAX);for (j=0; j<r; ++j)	{	for (k=0; k<c; ++k)	{		a=((int) i[j][k]);		if (a>max)	max=a;		}	}return max;}