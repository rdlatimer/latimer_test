#define mak_sort_functions #include "sort.h"#define memory#include "memory.h"#include <limits.h>				/* needed for RAND_MAX calls	*//*	copied almost directly from K&R p62/*	special thanks to D. L. Shell and K&R for reproducing it./*	vector v is sorted from v[0] ... v[n-1] into increasing order/*	Needs to be sent array (same will be returned sorted) and length n/*  Modified by P. Wagner to include sorting in decreasing order (_sort_dec) 09/01/2001/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//* integer version/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//* sort from low to high */int *ishellsort_inc (int *v, int n){int gap, i, j, temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}return v;}/*  sort from high to low */int *ishellsort_dec (int *v, int n){int gap, i, j, temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}return v;}void *ishellsort_inc_command (int *v, int n){int gap, i, j,  temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}}/*  sort from high to low */void *ishellsort_dec_command (int *v, int n){int gap, i, j, temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}}/* long version/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */long *lshellsort_inc (long *v, int n){int gap, i, j;long temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}return v;}void *lshellsort_inc_command (long *v, int n){int gap, i, j;long temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}}long *lshellsort_dec (long *v, int n){int gap, i, j;long temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}return v;}/* unsigned long version/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */unsigned long *ulshellsort_inc (unsigned long *v, int n){int gap, i, j;unsigned long temp;for (gap=n/2 ; gap>0 ; gap/=2)	for (i=gap ; i<n ; i++)		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;		}return v;}void *ulshellsort_inc_command (unsigned long *v, int n){int gap, i, j;unsigned long temp;for (gap=n/2 ; gap>0 ; gap/=2)	for (i=gap ; i<n ; i++)		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;		}}unsigned long *ulshellsort_dec (unsigned long *v, int n){int gap, i, j;unsigned long temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}return v;}void *ulshellsort_dec_command (unsigned long *v, int n){int gap, i, j;unsigned long temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}}/* float version/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */float *fshellsort_inc (float *v, int n){	int gap, i, j;	float temp;		for (gap=n/2 ; gap>0 ; gap/=2)		for (i=gap ; i<n ; i++)			for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {				temp = v[j];				v[j] = v[j+gap];				v[j+gap] = temp;			}		return v;}float *fshellsort_dec (float *v, int n){int gap, i, j;float temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}return v;}/* double version/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */double *dshellsort_inc (double *v, int n){	int gap, i, j;	double temp;		for (gap=n/2 ; gap>0 ; gap/=2)		for (i=gap ; i<n ; i++)			for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {				temp = v[j];				v[j] = v[j+gap];				v[j+gap] = temp;			}		return v;}double *dshellsort_dec (double *v, int n){	int gap, i, j;	double temp;		for (gap=n/2 ; gap>0 ; gap/=2)		for (i=gap ; i<n ; i++)			for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {				temp = v[j];				v[j] = v[j+gap];				v[j+gap] = temp;			}				for (i=0; i<n/2; ++i)	{		temp=v[(n-1)-i];		v[(n-1)-i]=v[i];		v[i]=temp;		}		return v;}void *dshellsort_inc_command (double *v, int n){int gap, i, j;double temp;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;			}		}	}}void *dshellsort_dec_command (double *v, int n){int gap, i, j;double temp;for (gap=n/2 ; gap>0 ; gap/=2)	for (i=gap ; i<n ; i++)		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			temp = v[j];			v[j] = v[j+gap];			v[j+gap] = temp;		}		for (i=0; i<n/2; ++i)	{	temp=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=temp;	}return v;}/* sort_decintbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/int	*sort_decintbydouble (int *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_decintbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/void	*sort_decintbydoublecommand (int *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}}/* sort_incintbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/void	*sort_incintbydoublecommand (int *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}}/* sort_deculongbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of unsigned long, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/unsigned long *sort_inculongbydouble (unsigned long *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_deculongbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of unsigned long, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/unsigned long *sort_deculongbydouble (unsigned long *i, double *k, int N){int	a, b, c, d;unsigned long	*j, *p;j=ulvector(N);j[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[j[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		j[b+1]=j[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	j[d+1]=c;	}p=ulvector(N);for (a=0; a<N; ++a)	p[a]=i[j[a]];for (a=0; a<N; ++a)	i[a]=p[a];free_ulvector(j);free_ulvector(p);return i;}/* sort_deculongbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of unsigned long, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/void *sort_inculongbydoublecommand (unsigned long *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}//return i;}/* sort_deculongbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of unsigned long, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/void *sort_deculongbydoublecommand (unsigned long *i, double *k, int N){int	a, b, c, d;unsigned long	*j, *p;j=ulvector(N);j[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[j[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		j[b+1]=j[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	j[d+1]=c;	}p=ulvector(N);for (a=0; a<N; ++a)	p[a]=i[j[a]];for (a=0; a<N; ++a)	i[a]=p[a];free_ulvector(j);free_ulvector(p);//return i;}/* sort_deculongbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of unsigned long, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/void *sort_inculongbyulongcommand (unsigned long *i, unsigned long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}//return i;}/* sort_deculongbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of unsigned long, the order of which gives the rank order/*			of elements on array k./*		k: an array of unsigned long numbers that will be used to rank i./*		N: the length of arrays i and k/***************************************************************************/void *sort_deculongbyulongcommand (unsigned long *i, unsigned long *k, int N){int	a, b, c, d;unsigned long	*j, *p;j=ulvector(N);j[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[j[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		j[b+1]=j[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	j[d+1]=c;	}p=ulvector(N);for (a=0; a<N; ++a)	p[a]=i[j[a]];for (a=0; a<N; ++a)	i[a]=p[a];free_ulvector(j);free_ulvector(p);//return i;}/* sort_decintbyint - produces an array sorted by an array of integers./* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/int	*sort_decintbyint (int *i, int *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_decintbylong - produces an array sorted by an array of long./* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/int	*sort_decintbylong (int *i, long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_declongbylong - produces an long array sorted by an array of long in descending order./* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/long *sort_declongbylong (long *i, long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_declongbydec - produces an long array sorted by an array of long in descending order./* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/long *sort_declongbydouble (long *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}long **sort_declongrowbylongrow(long **i, long **k, int r1, int r2, int n){int	a, b, c, d;i[r1][0]=0;for (a=1; a<n; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[r2][c]>k[r2][i[r1][b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[r1][b+1]=i[r1][b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[r1][d+1]=c;	}return i;}/* sort_declongbydec - produces an long array sorted by an array of long in descending order./* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/void *sort_declongbydoublecommand (long *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}void **sort_decullongrowbylongrow(unsigned long **i, long **k, int r1, int r2, int n){int	a, b, c, d;i[r1][0]=0;for (a=1; a<n; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[r2][c]>k[r2][i[r1][b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[r1][b+1]=i[r1][b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[r1][d+1]=c;	}//return i;}/* sort_decintbyint - produces an array sorted by an array of integers./* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/int	*sort_incintbyint (int *i, int *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_incintbydouble - produces an array sorted by an array of real numbers/* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of real numbers./*		N: the length of arrays i and k/***************************************************************************/int	*sort_incintbydouble (int *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_incintbylong - produces an array sorted by an array of long with smallest numbers first/* Requires:/*		i: an array of integers, the order of which gives the rank order/*			of elements on array k./*		k: an array of long numbers./*		N: the length of arrays i and k/***************************************************************************/int	*sort_incintbylong (int *i, long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_inclongbylong - sorts array k of unsigned long by array i of unsigned long./* Requires:/*		i: an array of long, the order of which gives the rank order/*			of elements on array k./*		k: an array of long, on which i is sorted./*		N: the length of arrays i and k/***************************************************************************/long	*sort_inclongbylong (long *i, long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_inclongbydouble - sorts array k of unsigned long by array i of unsigned long./* Requires:/*		i: an array of long, the order of which gives the rank order/*			of elements on array k./*		k: an array of double, on which i is sorted./*		N: the length of arrays i and k/***************************************************************************/long	*sort_inclongbydouble (long *i, double *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_inclongbyluong - sorts unsigned long array i by unsigned long array k./* Requires:/*		i: an array of long./*		k: an array of long, the order of which gives the rank order/*			of elements on array k./*		N: the length of arrays i and k/***************************************************************************/unsigned long *sort_inculongbyulong (unsigned long *i, unsigned long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort_deculongbyulong - sorts unsigned long array k by unsigned long array i./* Requires:/*		i: an array of long./*		k: an array of long, the order of which gives the rank order/*			of elements on array k./*		N: the length of arrays i and k/***************************************************************************/unsigned long *sort_deculongbyulong (unsigned long *i, unsigned long *k, int N){int	a, b, c, d;i[0]=0;for (a=1; a<N; ++a)	{	c=a;		d=a-1;	for (b=a-1; (b>=0 && k[c]>k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}return i;}/* sort columns within rows such that:	m[0] = 9 6 4 2	m[1] = 6 5 5 4	m[2] = 7 7 5 5	m[3] = 5 4 4 3goes to:	m[0] = 5 4 4 2	m[1] = 6 5 4 3	m[2] = 7 6 5 4	m[3] = 9 7 5 5/* Requires:/*		m: an matrix of long./*		row: the number of rows./*		col: the number of columns./***************************************************************************/long **sortcolinlmatrix_inc(long **m, int row, int col){int		c, r;long	*dummy;dummy=lvector(row);for (c=0; c<col; ++c)	{	for (r=0; r<row; ++r)		dummy[r]=m[r][c];		lshellsort_inc(dummy,row);		for (r=0; r<row; ++r)		m[r][c]=dummy[r];	}free_lvector(dummy);return m;}/* sortlmatrixoncol_dec - sorts long matrix m[r][c] by the elements in row m[srow]sortlmatrix(m,1,4,5) with m initially:	m[0] = 6 5 5 4 3	m[1] = 7 7 5 5 4	m[2] = 9 6 4 2 2	m[3] = 5 4 4 3 2yields:	m[0] = 5 4 4 3 2	m[1] = 6 5 5 4 3	m[2] = 7 7 5 5 4	m[3] = 9 6 4 2 2	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */void **sortlmatrixoncol_inc (long **m, int scol, int row, int col){int a, c, r;long *temp;temp=lvector(col); 	/* this stores the values presently being sorted, with temp[srow] 	*/					/* 	the value used to sort temp relative to other rows in m			*/for (r=1; r<row; ++r)	{	for (c=0; c<col; ++c)	temp[c]=m[r][c];	for (a=r-1; temp[scol]<m[a][scol] && a>=0; --a)	{		for (c=0; c<col; ++c)	m[a+1][c]=m[a][c];		}	if (a<(r-1))	for (c=0; c<col; ++c)		m[a+1][c]=temp[c];	}free_lvector(temp);//return m;}/* sortlmatrixoncol_dec - sorts long matrix m[r][c] by the elements in m[�][col]sortlmatrix(m,1,4,5) with m initially:	m[0] = 6 5 5 4 3	m[1] = 7 7 5 5 4	m[2] = 9 6 4 2 2	m[3] = 5 4 4 3 2yields:	m[0] = 9 6 4 2 2	m[1] = 7 7 5 5 4	m[2] = 6 5 5 4 3	m[3] = 5 4 4 3 2	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */void **sortlmatrixoncol_dec (long **m, int scol, int row, int col){int a, c, r;long *temp;temp=lvector(col); 	/* this stores the values presently being sorted, with temp[scol] 	*/					/* 	the value used to sort temp relative to other rows in m			*/for (r=1; r<row; ++r)	{	for (c=0; c<col; ++c)	temp[c]=m[r][c];	for (a=r-1; temp[scol]>m[a][scol] && a>=0; --a)	{		for (c=0; c<col; ++c)	m[a+1][c]=m[a][c];		}	if (a<(r-1))	for (c=0; c<col; ++c)		m[a+1][c]=temp[c];	}free_lvector(temp);//return m;}/* sortulmatrixoncol_inc - sorts long matrix m[r][c] by the elements in row m[srow]sortulmatrix(m,1,4,5) with m initially:	m[0] = 6 5 5 4 3	m[1] = 7 7 5 5 4	m[2] = 9 6 4 2 2	m[3] = 5 4 4 3 2yields:	m[0] = 5 4 4 3 2	m[1] = 6 5 5 4 3	m[2] = 7 7 5 5 4	m[3] = 9 6 4 2 2input:	m: unsigned long matrix to be sorted	scol: column on which to be sorted (0�[col-1])	row: number of rows	col: number of columns	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */void **sortulmatrixoncol_inc (unsigned long **m, int scol, int row, int col){int a, c, r;unsigned long *temp;temp=ulvector(col); 	/* this stores the values presently being sorted, with temp[srow] 	*/					/* 	the value used to sort temp relative to other rows in m			*/for (r=1; r<row; ++r)	{	for (c=0; c<col; ++c)	temp[c]=m[r][c];	for (a=r-1; temp[scol]<m[a][scol] && a>=0; --a)	{		for (c=0; c<col; ++c)	m[a+1][c]=m[a][c];		}	if (a<(r-1))	for (c=0; c<col; ++c)		m[a+1][c]=temp[c];	}free_ulvector(temp);}/* sortulmatrixoncol_dec - sorts long matrix m[r][c] by the elements in m[�][col]sortlmatrix(m,1,4,5) with m initially:	m[0] = 6 5 5 4 3	m[1] = 7 7 5 5 4	m[2] = 9 6 4 2 2	m[3] = 5 4 4 3 2yields:	m[0] = 9 6 4 2 2	m[1] = 7 7 5 5 4	m[2] = 6 5 5 4 3	m[3] = 5 4 4 3 2	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */void **sortulmatrixoncol_dec (unsigned long **m, int scol, int row, int col){int a, c, r;unsigned long *temp;temp=ulvector(col); 	/* this stores the values presently being sorted, with temp[scol] 	*/					/* 	the value used to sort temp relative to other rows in m			*/for (r=1; r<row; ++r)	{	for (c=0; c<col; ++c)	temp[c]=m[r][c];	for (a=r-1; temp[scol]>m[a][scol] && a>=0; --a)	{		for (c=0; c<col; ++c)	m[a+1][c]=m[a][c];		}	if (a<(r-1))	for (c=0; c<col; ++c)		m[a+1][c]=temp[c];	}free_ulvector(temp);}/* sortdmatrixoncol_dec - sorts double matrix m[r][c] by the elements in row m[srow]sortdmatrix(m,1,4,5) with m initially:	m[0] = 6 5 5 4 3	m[1] = 7 7 5 5 4	m[2] = 9 6 4 2 2	m[3] = 5 4 4 3 2yields:	m[0] = 5 4 4 3 2	m[1] = 6 5 5 4 3	m[2] = 7 7 5 5 4	m[3] = 9 6 4 2 2	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */double **sortdmatrixoncol_inc (double **m, int scol, int row, int col){int a, c, r;double *temp;temp=dvector(col); 	/* this stores the values presently being sorted, with temp[srow] 	*/					/* 	the value used to sort temp relative to other rows in m			*/for (r=1; r<row; ++r)	{	for (c=0; c<col; ++c)	temp[c]=m[r][c];	for (a=r-1; temp[scol]<m[a][scol] && a>=0; --a)	{		for (c=0; c<col; ++c)	m[a+1][c]=m[a][c];		}	if (a<(r-1))	for (c=0; c<col; ++c)		m[a+1][c]=temp[c];	}free_dvector(temp);return m;}/* sortdmatrixoncol_dec - sorts double matrix m[r][c] by the elements in m[�][col]sortdmatrix(m,1,4,5) with m initially:	m[0] = 6 5 5 4 3	m[1] = 7 7 5 5 4	m[2] = 9 6 4 2 2	m[3] = 5 4 4 3 2yields:	m[0] = 9 6 4 2 2	m[1] = 7 7 5 5 4	m[2] = 6 5 5 4 3	m[3] = 5 4 4 3 2	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */double **sortdmatrixoncol_dec (double **m, int scol, int row, int col){int a, c, r;double *temp;temp=dvector(col); 	/* this stores the values presently being sorted, with temp[scol] 	*/					/* 	the value used to sort temp relative to other rows in m			*/for (r=1; r<row; ++r)	{	for (c=0; c<col; ++c)	temp[c]=m[r][c];	for (a=r-1; temp[scol]>m[a][scol] && a>=0; --a)	{		for (c=0; c<col; ++c)	m[a+1][c]=m[a][c];		}	if (a<(r-1))	for (c=0; c<col; ++c)		m[a+1][c]=temp[c];	}free_dvector(temp);		return m;}void sortdec_ulvectorsVWbyV (unsigned long *v, unsigned long *w, int n){int gap, i, j;unsigned long tempv, tempw;for (gap=n/2 ; gap>0 ; gap/=2)	{	for (i=gap ; i<n ; i++)	{		for (j=i-gap ; j>=0 && v[j]>v[j+gap] ; j-=gap) {			tempv = v[j];			v[j] = v[j+gap];			v[j+gap] = tempv;			tempw = w[j];			w[j] = w[j+gap];			w[j+gap] = tempw;			}		}	}		for (i=0; i<n/2; ++i)	{	tempv=v[(n-1)-i];	v[(n-1)-i]=v[i];	v[i]=tempv;	tempw=w[(n-1)-i];	w[(n-1)-i]=w[i];	w[i]=tempw;	}//return ();}/* sort_inclongbydouble - sorts array k of unsigned long by array i of unsigned long./* Requires:/*		i: an array of long, the order of which gives the rank order/*			of elements on array k./*		k: an array of double, on which i is sorted./*		N: the length of arrays i and k/***************************************************************************/void *sort_inclvectorbydvector (long *i, double *k, int N){int	a, b, c, d;long	*j;j=lvector(N);for (a=0; a<N; ++a)	j[a]=i[a];i[0]=0;for (a=1; a<N; ++a)	{	c=a;	d=a-1;	for (b=a-1; (b>=0 && k[c]<k[i[b]]); --b)	{		d=b-1;		if (b==-1)	break;	/* this should not be necessary but the computer is stupid */		i[b+1]=i[b];		if (b==0)	{			break;	/* this should not be necessary but the computer is stupid */			}		}	i[d+1]=c;	}for (a=0; a<N; ++a)	i[a]=j[i[a]];free_lvector(j);}