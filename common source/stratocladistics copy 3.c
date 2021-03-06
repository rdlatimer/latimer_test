/* Routines to date nodes and calculate the gaps implied by trees     */
/*		Written by Peter Wagner 01/95
/*		Updated:	03/96
/*					04/96
/*					02/97
/*					09/97
/*					08/98
/*					09/02
/*					01/03
/*					06/03 - interger arrays converted to long for dates
/*					03/05 - stratocompatibility added
/*					09/05 - more stratocompatibility added
/*					02/09 - final assault on stratocompatibility launched!
/**********************************************************************/

#define compatibility_functions
#include "compatibility_functions.h"
#define historicaldiversity
#include "historicaldiversity.h"
#define matrixanalysis
#include "matrixanalysis.h"
#define matrixchange
#include "matrixchange.h"
#define matrixreading
#include "matrixreading.h"
#define memory
#include "memory.h"
#define minmax
#include "minmax.h"
#define MonteCarloPhylogenyFunctions
#include "MonteCarloPhylogenyFunctions.h"
#define Optimization
#include "optimization.h"
#define sort
#include "sort.h"
#define stratocladistics
#include "stratocladistics.h"
#define TreeRead
#include "tree_read.h"

/** dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
long *dateclade(long *fa, long **tree, int clades, int notu)
{
int	node, b, sp, htu;
long *F;

F=lvector(notu+clades);

for (sp=0; sp<notu; ++sp)				F[sp]=fa[sp];

for (node=clades-1; node>=0; --node)	{
	F[htu=node+notu]=F[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (F[htu]>F[sp])	F[htu]=F[sp];
		}
	}
	
return F;
}

/** dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
long **datecladefull(long **range, long **tree, int clades, int notu)
{
int	b;
long *fa, *dates, **clrange;

fa=lvector(notu);
for (b=0; b<notu; ++b)	fa[b]=range[b][0];
dates=dateclade(fa,tree,clades,notu);
free_lvector(fa);

clrange=lmatrix(notu+clades,2);
for (b=0; b<notu+clades; ++b)	{
	if (b<notu)	{
		clrange[b][0]=range[b][0];
		clrange[b][1]=range[b][1];
		}
	else
		clrange[b][1]=dates[b];
	}

/* next step: drage clrange[b][0] to clrange[anc][1] for clades, where anc is the ancestral node of b	*/
return clrange;
}

/** datecladesim - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	tree - a matrix containing the tree structure plus first-last appearance data;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/

/** dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
double *datecladereal(double *fa, long **tree, int clades, int notu)
{
int	node, b, sp, htu;
double *F;

F=dvector(notu+clades);
b=maxdarray(fa,notu);

for (sp=0; sp<notu; ++sp)				F[sp]=fa[sp];

for (node=clades-1; node>=0; --node)	{
	F[htu=node+notu]=F[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (F[htu]>F[sp])	F[htu]=F[sp];
		}
	}
	
return F;
}


/** datecladeFile - Finds age of clades based on the oldest member in the clade and an fa file
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	fa - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
long* datecladefile(char *Origins, long **tree, int clades, int notu)
{
int	 node, b, sp, htu;
long *fa;
FILE *fopen();
FILE *FKA;

fa=lvector(notu+clades);
FKA = fopen(Origins,"r");
for (b=0; b<notu; ++b)
	fscanf(FKA,"%i",&fa[b]);
fclose(FKA);

b=maxlarray(fa,notu);

for (node=clades-1; node>=0; --node)	{
	fa[htu=node+notu]=fa[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (fa[htu]>fa[sp])	fa[htu]=fa[sp];
		}
	}
	
return fa;
}


/** dateclade - Finds age of clades based on the oldest member in the clade
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	F - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
double *datecladerealfile(char *Origins, long **tree, int clades, int notu)
{
int	 node, b, sp, htu;
double *fa;
FILE *fopen();
FILE *FKA;

fa=dvector(notu+clades);
FKA = fopen(Origins,"r");
for (b=0; b<notu; ++b)
	fscanf(FKA,"%lf",&fa[b]);
fclose(FKA);

for (node=clades-1; node>=0; --node)	{
	fa[htu=node+notu]=fa[sp=tree[node][1]];
	for (b=2; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (fa[htu]>fa[sp])	fa[htu]=fa[sp];
		}
	}
	
return fa;
}

/** datecladeAdd - Finds minimum divergence dates assuming monophyletic groups
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	fa - First Appearances (these are modified along the way)
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
****************************************************************************/
void datecladeadd(long *fa, long **tree, int clades, int notu)
{
int	node, b, sp;

for (node=clades-1; node>=0; --node)	{
	fa[node+notu]=100000;
	for (b=1; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (fa[node+notu]>fa[sp])	fa[node+notu]=fa[sp];
		}
	}
}

/** datecladeExtra - Finds minimum divergence dates allowing for ancestors
/* Requires:
/*	Diverg - A matrix of divergence times between taxa
/*	tree - a matrix containing the tree structure;
/*	div2 - diversity of each clade;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	fa - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
long *datecladeextra(int **Diverg, long *fa, long **tree, int clades, int notu)
{
int		a, node, sp1, sp2, both;
long	**VennTree;	/* tree listing all taxa in the node */

VennTree=lmatrix(notu, notu);
for (node=0; node<clades; ++node)
	for (a=0; a<=tree[node][0]; ++a)	VennTree[node][a]=tree[node][a];
VennTree=clademember(VennTree, notu, clades);

for (sp1=0; sp1<notu-1; ++sp1)	{
	for (sp2=sp1+1; sp2<notu; ++sp2)	{
		if (Diverg[sp1][sp2]!=0)	{
			/* work down the tree */
			/* the highest node with both sp1 and sp2 now gets fa[node+notu]=Diverg[sp1][sp2]	*/
			both=0;
			for (node=clades-1; both<2; --node)	{
				both=0;
				/* increment both whenever there is a match 	*/
				/* if it equals two, then there are two matches */
				for (a=1; a<=VennTree[node][0]; ++a)	{
					if (VennTree[node][a]==sp1 || VennTree[node][a]==sp2)	++both;
					if (both==2)	{
						a=VennTree[node][0];
						if (fa[node+notu]>Diverg[sp1][sp2])
							fa[node+notu]=Diverg[sp1][sp2];
						}	/* register fa for node, and end search */
					}	/* end search of node */
				}
			}	/* only bother if we have a divergence date of note */
		}
	}
	
free_lmatrix(VennTree, notu, notu);
return fa;
}

long *datecladediverge(long **DT, long *fa, long **tree, int clades, int notu)
{
int sp1, sp2, max, anc;
int **LCA;

max=DT[0][0];
LCA=lastcommonancestormatrix(tree,clades,notu);

for (sp1=0; sp1<notu+clades; ++sp1)	fa[sp1]=max;
for (sp1=0; sp1<notu-1; ++sp1)	{
	for (sp2=sp1+1; sp2<notu; ++sp2)	{
		while (LCA[sp1][sp2]==-1 && sp2<notu)	++sp2;
		if (sp2>=notu)	break;
		anc=LCA[sp1][sp2]+notu;
		if (fa[anc]>DT[sp1][sp2])	fa[anc]=DT[sp1][sp2];
		}
	}
free_imatrix(LCA,notu,notu);
return fa;
}


/** BranchAgeBin - Finds age of branches using a discrete scale assuming no ancestors
/*		NOTE: if you have continuous time, use "BranchAgeCont"
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - First Appearances
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
****************************************************************************/
long *blbinnaive(long **tree, long *fa, int clades, int notu)
{
int	 a, b, sp, anc;
long *bt;

bt=lvector(notu+clades);

for (a=0; a<clades; ++a)	{
	anc=a+notu;
	for (b=1; b<=tree[a][0]; ++b)	{
		sp=tree[a][b];
		bt[sp]=fa[sp]-fa[anc];
		}
	}
	
return bt;
}

/** BranchAgeBin - Finds age of branches using a continuous scale assuming no ancestors
/*		NOTE: if you have time bins (e.g., stages), use "BranchAgeBin"
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - First Appearances
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
****************************************************************************/
double *blmanaive(long **tree, double *fa, int clades, int notu)
{
int	 a, b, sp, anc;
double	 *bt;

bt=dvector(notu+clades);

for (a=0; a<clades; ++a)	{
	anc=a+notu;
	for (b=1; b<=tree[a][0]; ++b)	{
		sp=tree[a][b];
		bt[sp]=fa[sp]-fa[anc];
		}
	}
	
return bt;
}

/** blma - Tallies branch lengths in real time allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	bt - range extension on each branch.
****************************************************************************/
double *branchlngma(long **tree, double **ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
double *bt;
double *F, *fa;

fa=dvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
bt=dvector(notu+clades);

F=datecladereal(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		bt[sp]=F[anc]-F[sp];
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	bt[sp]=(F[sp]-ranges[anc][1])-1;
			}
		else
			bt[sp]=F[sp]-F[notu+nd];
		}
	}
free_dvector(F);
free_dvector(fa);
return bt;
}


/** CladeSurvive - Determines how long clades last
/* Requires:
/*	la - Last Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	L - modified array la, with elements notu�notu+clades now filled.
****************************************************************************/
long *CladeSurvive(long *la, long **tree, int clades, int notu)
{
int node, b, sp;
long *L;

L=lvector(notu+clades);

for (sp=0; sp<notu; ++sp)				L[sp]=la[sp];
for (sp=notu; sp<(notu+clades); ++sp)	L[sp]=0;

for (node=clades-1; node>=0; --node)	{
	for (b=1; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (L[node+notu]<L[sp])	L[node+notu]=L[sp];
		}
	}

return L;
}
/** CladeSurviveAdd - Finds minimum divergence dates assuming monophyletic groups
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	la - Last Appearances
/*	tree - a matrix containing the tree structure;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	L - modified array la, with elements notu�notu+clades now filled.
****************************************************************************/
void CladeSurviveAdd(long *la, long **tree, int clades, int notu)
{
int node, b, sp;

for (node=clades-1; node>=0; --node)	{
	for (b=1; b<=tree[node][0]; ++b)	{
		sp=tree[node][b];
		if (la[node+notu]<la[sp])	la[node+notu]=la[sp];
		}
	}
}

/** temporalbranchlength - Routine to get basic branch lengths in time units.  
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	bt - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
double *temporalbranchlength(long **tree, double *fa, int clades, int notu)
{
int	cl, sp, sc;
double	*bt;

bt=dvector(notu+clades);

for (cl=clades-1; cl>=0; --cl)	{
	for (sc=1; sc<=tree[cl][0]; ++sc)	{
		sp=tree[cl][sc];

		bt[sp]=fa[sp]-fa[cl+notu];
		if (bt[sp]<0)	bt[sp]*=-1;

		}
	}

return bt;
}

/** temporalbranchlength - Routine to get basic branch lengths in time units.  
/*		NOTE: if you have branch lengths, use "datecladesExtra - this allows for ancestors
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	clades - number of nodes (clades)
/*	notu - number of taxa;
/* Notes:
/*	clades are treated as taxa in array fa and given an htu number of node # + #notu;
/* Returns:
/*	bt - modified array fa, with elements notu�notu+clades now filled.
****************************************************************************/
double **temporalbranchlengthequal(long **tree, long *fba, double *byr, int clades, int notu)
{
int		a,b,c,d,g,s;
int		anc, cl, sp, ties;
int		*lumped;
long	*ancestors, *height, *taxa, *cladebins;
long	**vtree;
double	en, f, l;
double	*bt;
double	**phy;

/* get list of ancestors */
ancestors=((long*) listancestor(tree,clades,notu));
height=((long*) patristicheight(tree,clades,notu));
bt=dvector(notu+clades);
vtree=lmatrix(clades,notu);
equallmatrix(vtree,tree,clades,notu);
vtree=clademember(vtree,notu,clades);

/* initially set birth to middle of bin	*/
/*for (sp=0; sp<notu; ++sp)		bt[sp]=(byr[fba[sp]]+byr[1+fba[sp]])/2;	*/

cladebins=dateclade(fba, tree, clades, notu);
/*cladeages=datecladereal(fa,tree,clades,notu);	*/
for (sp=0; sp<clades+notu; ++sp)		bt[sp]=(byr[cladebins[sp]]+byr[1+cladebins[sp]])/2;

lumped=ivector(clades);

taxa=lvector(notu);
for (sp=0; sp<notu; ++sp)	taxa[sp]=sp;
taxa=sort_declongbylong(taxa,height,notu);
taxa=sort_declongbylong(taxa,fba,notu);

for (s=0; s<notu; ++s)	{
	sp=taxa[s];
	cl=(anc=ancestors[sp])-notu;
	while (cl<0 && s<notu)	{
		sp=taxa[++s];
		cl=(anc=ancestors[sp])-notu;
		}
	if (s>=notu)	break;
	f=bt[sp];
	l=bt[anc];
	
	/* If a taxon's clade already has been adjusted, then the ancestor might be younger		*/
	/* make sure that all descendants of that ancestor are at least as old as that ancestor	*/
	if (f<l)	{
		bt[sp]=f=l;
/*		cl=anc-notu;
/*		for (d=1; d<=vtree[cl][0]; ++d)	{
/*			g=vtree[cl][d];
/*			if (bt[g]<bt[anc]
/*			}	*/
		}
	ties=0;
	lumped[0]=sp;
	while (f==l && cl>0)	{
		lumped[++ties]=anc;
		anc=ancestors[cl+notu];
		cl=anc-notu;
		f=l;

		if (f<l)	{
			f=l;
			for (a=0; a<=ties; ++a)	{
				b=lumped[a];
				bt[b]=l;
				}
			}

		if (anc>=0)	l=bt[anc];
		else		l=0;
		}
	lumped[ties+1]=anc;

	if (ties>0)	{
		/* reset birth to the end of the bin	*/
		bt[sp]=byr[1+fba[sp]];
		
		a=fba[sp]+1;
		en=(byr[a]-bt[anc])/((double) (ties+1));
		
		/* if the first clade already has been adjusted, then start from that; otherwise, start from its ancestor	*/
		a=lumped[ties];
		if (bt[a]!=(byr[fba[a]]+byr[fba[a]+1])/2)	--ties;
		
		for (a=ties; a>=0; --a)	{
			b=lumped[a];
			c=lumped[a+1];
			/* if b is a clade,then make sure all of its descendants have their ages modified accordingly	*/
			if (b>notu && a>1)	{
				cl=b-notu;
				for (d=1; d<=vtree[cl][0]; ++d)	{
					g=vtree[cl][d];
					if (bt[g]==bt[b])	bt[g]=bt[c]+en;
					}
				}
			bt[b]=bt[c]+en;
			}
		}
	}
	
/*for (sp=0; sp<notu+clades; ++sp)	{
	if (ancestors[sp]>=notu)	{
		anc=ancestors[sp];
		bt[1][sp]=bt[0][sp]-bt[0][anc];
		}
	}	*/

phy=dmatrix(notu+clades,2);
for (a=0; a<notu+clades; ++a)	{
	phy[a][0]=bt[a];
	anc=ancestors[a];
	phy[a][1]=bt[a]-bt[anc];
	}

free_dvector(bt);
free_lvector(height);
free_lvector(taxa);
free_lmatrix(vtree,clades,notu);
free_ivector(lumped);
free_lvector(ancestors);
free_lvector(cladebins);
return phy;
}


double *temporalbranchlengthApo(long **tree, long *fa, int *bl, int clades, int notu)
{
int	cl, sp, sc;
int	*anc;
double	*bt;

anc=ivector(clades);
bt=dvector(notu+clades);

for (cl=clades-1; cl>=0; --cl)	{
	for (sc=1; sc<=tree[cl][0]; ++sc)	{
		sp=tree[cl][sc];
		if (bl[sp]==0 && sp<notu)	anc[cl]=1;
		}
	}

for (cl=clades-1; cl>=0; --cl)	{
/* if none of the sister taxa are ancestral, then it is the same as temporalbranchlength */
	if (anc[cl]==0)	{
		for (sc=1; sc<=tree[cl][0]; ++sc)	{
			sp=tree[cl][sc];
			if (sp<notu)
				bt[sp]=0-fa[cl+notu];
			else
				bt[sp]=fa[sp]-fa[cl+notu];
			}
		}
	/* if species is ancestral, then it can acquire apomorphies from the clade's
	origin until the derived taxa first appear */
	else	{
		for (sc=1; sc<=tree[cl][0]; ++sc)	{
			sp=tree[cl][sc];
			if (sp<notu && bl[sp]==0)
				bt[sp]=0;
			else if (sp<notu && bl[sp]>0)
				bt[sp]=0-fa[sp];
			}
		}
	}

return bt;
}
/*
double *AdjustTemporalBLAp(long **tree, long *fa, int *bl, double *bt,double rate, int clades, int notu)
{
int	sp;
double	t;
double *Clock;

Clock=dvector(clades+notu);

for (sp=0; sp<clades+notu; ++sp)
	Clock[sp]=(t=bl[sp])/rate;


free_dvector(Clock);
}	*/

/** naivestratdebt - Tallies stratigraphic gaps assuming no ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use naiveMIG for real numbers).
****************************************************************************/
int naivestratdebt(long **tree, long *fa, int notu)
{
int d, nd, sp, clades;
int	sd=0;
long *F;

clades=cladecountbytaxa(tree,notu);
F=dateclade(fa,tree,clades,notu);

for (nd=0; nd<clades; ++nd)	{
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		sd=sd+(F[sp]-F[nd+notu]);
		}
	}
return sd;
}

/** stratdebt - Tallies stratigraphic gaps in bins (e.g., stages) allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use MIGcalc for real numbers).
****************************************************************************/
int stratdebt(long **tree, long**ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
int	sd=0;
long *F, *fa;

fa=lvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
F=dateclade(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		sd=sd+(F[anc]-F[sp]);
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	sd=sd+((F[sp]-ranges[anc][1])-1);
			}
		else
			sd=sd+(F[sp]-F[notu+nd]);
		}
	}
free_lvector(F);
free_lvector(fa);
return sd;
}

/** stratdebtnodes - Tallies stratigraphic gaps in bins (e.g., stages) allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use MIGcalc for real numbers).
****************************************************************************/
long *stratdebtnodes(long **tree, long**ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
long	sd=0;
long *F, *fa;

fa=lvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
F=dateclade(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		sd=sd+(F[anc]-F[sp]);
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	sd=sd+((F[sp]-ranges[anc][1])-1);
			}
		else
			sd=sd+(F[sp]-F[notu+nd]);
		}
	}
free_lvector(F);
free_lvector(fa);
return fa;		/* CHANGE THIS!!!	*/
}

/** stratdebtreal - Tallies stratigraphic gaps in real time allowing for ancestors.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	la - an array of last appearance dates;
/*	bl - an array of branch lengths;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	sd - Number of gaps (in integers - use MIGcalc for real numbers).
****************************************************************************/
double stratdebtreal(long **tree, double **ranges, int *bl, int notu)
{
int d, nd, sp, clades, anc;
double	sd=0.0;
double *F, *fa;

fa=dvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];

clades=cladecountbytaxa(tree,notu);
F=datecladereal(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	anc=-1;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp<notu && bl[sp]==0)	{
			anc=sp;
			d=tree[nd][0];
			}
		}
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		/* if ancestor present and ranges overlap or abut, then no gap */
		if (anc>-1)	{
			/* if descendent precedes ancestor, then there is a gap    */
			if (F[sp]<F[anc])		sd=sd+(F[anc]-F[sp]);
			/* if ancestor disappears before desc appears, then a gap  */
			else if (F[sp]>ranges[anc][1])	sd=sd+((F[sp]-ranges[anc][1])-1);
			}
		else
			sd=sd+(F[sp]-F[notu+nd]);
		}
	}
free_dvector(F);
free_dvector(fa);
return sd;
}

/** SCI - Calculates Stratigraphic Consistency Index after Huelsenbeck 1994.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	SCI - consistent nodes / total nodes.
****************************************************************************/
double calcSCI(long **tree, long *fa, int notu)
{
int	d, dcl, sp, nd, d1, d2;
int	clades;
int	clfirst, spfirst, sim;
long *F;
double	con=0.0f, incon=0.0f, SCI=0.0f;

clades=cladecountbytaxa(tree,notu);

F=dateclade(fa,tree,clades,notu);

for (nd=clades-1; nd>=0; --nd)	{
	dcl=0;
	for (d=1; d<=tree[nd][0]; ++d)	{
		sp=tree[nd][d];
		if (sp>notu)	++dcl;
		}
	if (dcl>0)	{
		/* if it is a bifurcation, then its easy */
		if (tree[nd][0]==2)	{
			if (dcl==2)	{
				d1=tree[nd][1];
				d2=tree[nd][2];
				if (F[d1]==F[d2])	con=con+2.0;
				else	{
					con=con+1.0;
					incon=incon+1.0;
					}
				}
			else if (dcl==1)	{
				if (tree[nd][1]<notu)	{
					d1=tree[nd][1];
					d2=tree[nd][2];
					}
				else	{
					d1=tree[nd][2];
					d2=tree[nd][1];
					}
				if (F[d1]<=F[d2])	con=con+1;
				else				incon=incon+1;
				}
			}	/* end routine for bifurcation */

		/* polytomies are trickier */
		/* if a clade is part of a polytomy and there is any taxon as old	*/
		/*		as or older than that clade, then it is consistent			*/
		/*		NOTE: only one clade can be inconsistent; 			 		*/
		else	{
			clfirst=spfirst=0;
			for (d=1; d<=tree[nd][0]; ++d)
				spfirst=clfirst=clfirst+F[tree[nd][0]];
			
			for (d=1; d<=tree[nd][0]; ++d)	{
				sp=tree[nd][d];
				/* keep track of which TU's appear first */
				if (sp<notu && F[sp]<spfirst)	{
					spfirst=F[sp];
					}
				else if (sp>=notu)	{
					if (F[sp]<clfirst)	{
						clfirst=F[sp];
						sim=1;
						}
					else if (F[sp]==clfirst)	++sim;
					}
				}
			if (clfirst<spfirst && sim==1)	{
				con=con+(dcl-1);			/* at most, one node attached to a polytomy can be inconsistent */
				++incon;				
				}
			else	con=con+dcl;
			}
		}
	}

SCI=con/(con+incon);

free_lvector(F);
return SCI;
}

/** SCIm - Calculates modified Stratigraphic Consistency Index after Habib & Gittelman.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/* Returns:
/*	SCI - consistent nodes / total nodes.
****************************************************************************/
double calcSCIm(long **tree, long *fa, int notu)
{
int	b, c, cl, d, dcl, nd, scl, mxdc;
int	clades, *f1, *dc;
long *F, *fn, **treecl;
double	denom=0.0f, SCIm=0.0f;

clades=cladecountbytaxa(tree,notu);
treecl=cladesinclades(tree,clades,notu);
dc=ivector(clades);
mxdc=0;
for (c=0; c<clades; ++c)	{
	dc[c]=0;
	for (d=1; d<=tree[c][0]; ++d)	if (tree[c][d]>=notu)	++dc[c];
	if (dc[c]>mxdc)	mxdc=dc[c];
	}
f1=ivector(mxdc);

F=dateclade(fa,tree,clades,notu);
fn=lvector(clades);
for (c=0; c<clades; ++c)	fn[c]=F[c+notu];
free_lvector(F);

for (nd=0; nd<clades; ++nd)	{
	/* only examine nodes with sister nodes */
	while (dc[nd]<2 && nd<clades)	++nd;
	if (nd>=clades)	break;
	/* go through node and find daughter clades (converting to clade # as you go) */
	dcl=-1;
	for (d=1; d<=tree[nd][0]; ++d)	if (tree[nd][d]>=notu)	f1[++dcl]=tree[nd][d]-notu;

	/* Now, go through nodes derived from sister node and see how many appear later	*/
	for (d=0; d<dc[nd]; ++d)	{
		cl=f1[d];
		for (c=0; c<dc[nd]; ++c)	{
			if (c!=d)	{
				scl=f1[c];						/* sister clade */
				denom+=treecl[scl][0];
				if (fn[cl]<=fn[scl])	SCIm+=treecl[scl][0];
				else	{
					/* tally all of sister clades descendant nodes that are as old as or younger */
					for (b=1; b<=treecl[scl][0]; ++b)	{
						dcl=treecl[scl][b];		/* daughter clade of sister clade */
						if (fn[dcl]>=fn[cl])	++SCIm;
						}
					}
				/* tally all of sister clades descendant nodes */
				}	/* finish test on condition that c and d are not the same */
			}	/* finish search of sister clades */
		}	/* finish comparisons of clades within node nd */
	}	/* finish search of tree */

SCIm/=denom;

free_ivector(f1);
free_ivector(dc);
free_lvector(fn);
free_lmatrix(treecl,clades,clades);
return SCIm;
}

/* GER - Calculates Gap Excess Ratio after Wills 1998, 1999.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use GERcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	GER: 1 - (implied - minimum gaps)/(max gaps - min gaps).
****************************************************************************/
double calcGER(long **tree, long *fa, int notu)
{
int		d;
int		clades;
double	mngap=0.0f, mxgap=0.0f, gaps=0.0f;
double	GER=0.0f;

clades=cladecountbytaxa(tree,notu);

mngap=maxlarray(fa,notu)-minlarray(fa,notu);
mxgap=larraytotal(fa,notu)-(minlarray(fa,notu)*notu);

d=naivestratdebt(tree,fa,notu);
gaps=d;

GER=(mxgap-gaps)/(mxgap-mngap);
return GER;
}

/** RCI - Calculates Relative Completeness Index after Benton 1993.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	RCI - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcRCI(long **tree, long**ranges, int notu)
{
int		d, sp;
int		clades;
long	*fa;
double	gaps=0.0f, SRL=0.0f;
double	RCI=0.0f;

fa=lvector(notu);
for (d=0; d<notu; ++d)	fa[d]=ranges[d][0];
clades=cladecountbytaxa(tree,notu);
for (sp=0; sp<notu; ++sp)	SRL=SRL+(1+(ranges[sp][1]-ranges[sp][0]));

d=naivestratdebt(tree,fa,notu);
gaps=d;

RCI=1-(gaps/SRL);
free_lvector(fa);
return RCI;
}

/** MSM - Manhattan Stratigraphic Metric after Siddall ????.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	MSM - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcMSM(long **tree, long *fa, int notu)
{
int		d;
int		clades;
double	gaps=0.0f, mngap=0.0f;
double	MSM=0.0f;

clades=cladecountbytaxa(tree,notu);

mngap=maxlarray(fa,notu)-minlarray(fa,notu);
	
d=naivestratdebt(tree,fa,notu);
gaps=d;
MSM=(mngap/gaps);

return MSM;
}

/* GER - Calculates Gap Excess Ratio after Wills 1998, 1999.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a PHYLOGENY;
/* Returns:
/*	GER: 1 - (implied - minimum gaps)/(max gaps - min gaps).
****************************************************************************/
double calcGERcor(long **tree, long **ranges, int *bl, int notu)
{
int		d;
int		clades;
int		onset, end;
long	*richness;
double	mngap=0.0f, mxgap=0.0f, gaps=0.0f;
double	GER=0.0f;

clades=cladecountbytaxa(tree,notu);

mxgap=sumlmatrixcol(ranges,notu,0)-(minlmatrixcol(ranges,notu,0)*notu);

onset=minlmatrixcol(ranges,notu,0);
end=maxlmatrixcol(ranges,notu,1);

/* necessary gap exists only if no taxa are sampled over an interval */
richness=richnesstally(ranges, notu);
for (d=onset; d<=end; ++d)
	if (richness[d]==0)	++mngap;

d=stratdebt(tree,ranges,bl,notu);
gaps=d;

free_lvector(richness);

GER=(mxgap-gaps)/(mxgap-mngap);
return GER;
}

/** RCI - Calculates Relative Completeness Index after Benton 1993.  
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	RCI - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcRCIcor(long **tree, long**ranges, int *bl, int notu)
{
int		d, sp;
int		clades;
double	gaps=0.0f, SRL=0.0f;
double	RCI=0.0f;

clades=cladecountbytaxa(tree,notu);
for (sp=0; sp<notu; ++sp)	SRL=SRL+(1+(ranges[sp][1]-ranges[sp][0]));

d=stratdebt(tree,ranges,bl,notu);
gaps=d;

RCI=1-(gaps/SRL);
return RCI;
}

/** MSM - Manhattan Stratigrahic Metric after Siddall 1997.  
/* NOTE: this metric is completely nonsensical and should not be used.
/* Requires:
/*	tree - a matrix containing the tree structure;
/*	fa - an array of first appearance dates;
/*	notu - number of taxa;
/* Note: pay attention to how this handles polytomies
/*		Also, a "consistent" node is a clade with an older sister taxon
/*		NOTE: this calculates for a CLADOGRAM, not a phylogeny;
/*			Use RCIcor to calculate GER for a phylogeny (with branch lengths);
/* Returns:
/*	MSM - difference bn. observed and minimum gaps over difference bn. max and min.
****************************************************************************/
double calcMSMcor(long **tree, long**ranges, int *bl, int notu)
{
int		d;
int		clades;
int		onset, end;
long	*richness;
double	gaps=0.0f, mngap=0.0f;
double	MSM=0.0f;

clades=cladecountbytaxa(tree,notu);

onset=minlmatrixcol(ranges,notu,0);
end=maxlmatrixcol(ranges,notu,1);

richness=richnesstally(ranges, notu);

for (d=onset; d<=end; ++d)
	if (richness[d]==0)	++mngap;
	
d=stratdebt(tree,ranges,bl,notu);
gaps=d;
if (d==0)	MSM=0.0;
else		MSM=(mngap/gaps);

free_lvector(richness);
return MSM;
}


/** stratcompatibilty - determines whether each character in a matrix is compatible with stratigraphy or not.  
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	types - type for each character (0: ordered; 1: unordered);
/*	nchars - number of characters;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stratcomp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
****************************************************************************/
unsigned long *stratcompatibility (long **ranges, long **matrix, int *types, int nchars, int notu, int UNKNOWN, int INAP)
{
int				a, c, q, s, t;
int				fs, ls;
int				mnst, mxst, sts;
int				onset, end;
int				*diverse, **stord;
long			*rich;
unsigned long	*stratcomp;

onset=minlmatrixcol(ranges,notu,0);
end=maxlmatrixcol(ranges,notu,1);

rich=richnesstally(ranges,notu);

diverse=ivector(end+1);
stratcomp=ulvector(nchars);

sts=maxinclmatrix(matrix,notu,nchars,UNKNOWN,INAP);
stord=imatrix(sts+1,3);

for (c=0; c<nchars; ++c)	{
	/* find minimum and maximum states	*/
	mnst=colminclmatrix(matrix,notu,c,UNKNOWN,INAP);
	mxst=colmaxclmatrix(matrix,notu,c,UNKNOWN,INAP);
	
	stratcomp[c]=1;
	for (s=mnst; s<=mxst && stratcomp[c]==1; ++s)	{
		clearivector(diverse,end+1,0);
		for (t=0; t<notu; ++t)	{
			if (matrix[t][c]==s)	for (a=ranges[t][0]; a<=ranges[t][1]; ++a)	++diverse[a];
			}
		q=0;
		for (a=onset; a<=end; ++a)	{
			if (diverse[a]>0 && q==0)						q=1;		/* first find of character 			*/
			else if ((diverse[a]==0 && rich[a]>0) && q==1)	q=2;		/* first absence after being found	*/
			else if (diverse[a]>0 && q==2)	{							/* a gap is demonstrated			*/
				a=end;
				stratcomp[c]=0;
				}
			}	/* search diverse array for gaps between first and last sightings of state s */
		}	/* search each state s in character c for gaps	*/ 
	
	/* one more test for ordered multistates	*/
	if (stratcomp[c]==1 && (types[c]==0 && (mxst-mnst)>1))	{
		/* find the oldest state(s)	*/
		clearimatrix(stord,sts,3,RAND_MAX);
		for (s=0; s<=sts; ++s)	{
			for (t=0; t<notu; ++t)	{
				if (matrix[t][c]==s)	{
					if (stord[s][0]==RAND_MAX)		{
						stord[s][0]=s;
						stord[s][1]=ranges[t][0];
						stord[s][2]=ranges[t][1];
						}
					else {
						if (ranges[t][0]<stord[s][1])	stord[s][1]=ranges[t][0];
						if (ranges[t][1]>stord[s][2])	stord[s][2]=ranges[t][1];
						}
					}	/* end case where taxon has the right state */
				}	/* end search of taxa	*/
			
			if (fs==RAND_MAX && stord[s][1]==onset)	fs=s;
			}	/* end search through states	*/
		
		/* now make sure that everything is adjacent	*/
		/* fs gives the first appearing state; all states above it and below it must appear after that BUT without gap between it and fs	*/
		/*    Also, each state moving away from the starting point has to appear later BUT without gap between it and �prior� state			*/
		ls=fs;	/* because we might skip states, look at the last coded state; this ignores �gaps� in continuous character sequence	*/
		for (s=fs+1; s<=mxst && stratcomp[c]==1; ++s)	{
			while (stord[s][0]==RAND_MAX && s<=mxst)	++s;
			if (s>mxst)	break;
			if (stord[s][1]<stord[ls][1])				stratcomp[c]=0;
			else if (stord[s][1]>(stord[ls][2]+1))		stratcomp[c]=0;
			ls=s;
			}

		ls=fs;	/* because we might skip states, look at the last coded state; this ignores �gaps� in continuous character sequence	*/
		for (s=fs-1; s>=mnst && stratcomp[c]==1; --s)	{
			while (stord[s][0]==RAND_MAX && s>=mnst)	--s;
			if (s<mnst)	break;
			if (stord[s][1]<stord[ls][1])				stratcomp[c]=0;
			else if (stord[s][1]>(stord[ls][2]+1))		stratcomp[c]=0;
			ls=s;
			}
		}	/* end case of ordered multistate	*/
	}	/* test each character for stratigraphic compatibility	*/

free_imatrix(stord,sts,3);
free_ivector(diverse);
free_lvector(rich);
return stratcomp;
}

/** stratcompatfull - determines whether a compatible character pair are compatible with stratigraphy or not.  
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	types - type for each character (0: ordered; 1: unordered);
/*	states - number of states for ach character;
/*	ch1 - 1st character being compared;
/*	ch2 - 2nd character being compared;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	strcmp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		strcmp[0]: character 1
		strcmp[1]: character 2
		strcmp[2]: character 1 states
		strcmp[3]: character 2 states
		strcmp[4]: state comparisons (2 for binary; possibly more for multistate)
		strcmp[5]: first appearances compatible with stratigraphy
		strcmp[6]: pair synoptic ranges completely compatible with stratigraphy
		strcmp[7]: taxon synoptic ranges completely compatible with stratigraphy
		strcmp[8]: divergent (0) vs. hierarchical (1) stratigraphic compatibility
****************************************************************************/
unsigned long *stratcompatfull(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, f;
int	cc, sp, max, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs;
int cl[20], cr[20];
int	keyst1, keyst2, sw1, sw2;				/* these are the "swing" states between characters + important info	*/
unsigned long *strcmp;
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **pairrng, **prfl, **pairfnd, *gaps;

max=maxlmatrixcol(ranges, notu, 1);

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

strcmp=ulvector(9);
clearulvector(strcmp,9,0);

strcmp[0]=ch1;
strcmp[1]=ch2;
strcmp[2]=states[ch1];
strcmp[3]=states[ch2];

combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		if (pairrng[a][0]<max && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			}
		}
	}

sc1=sc2=0;
for (st1=0; st1<states[ch1]; ++st1)	if (combos[0][st1]>1)	++sc1;
for (st2=0; st2<states[ch2]; ++st2)	if (combos[1][st2]>1)	++sc2;

pairs=lmatrix(prs=(states[ch1]+states[ch2]),2);	/* state pairs between two characters		*/
prfl=lmatrix(prs,2);								/* first/last appearance of character pairs	*/
gaps=lvector(prs);

/*	Find each state paired with multiple states: these are the crux	*/
keyst1=-1;

/* play off of the character with the most swing states	*/
if (sc2>sc1)	{
	free_lmatrix(pairrng,states[ch1]*states[ch2],2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
	free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
	c=ch1;
	ch1=ch2;
	ch2=c;
	c=sc1;
	sc1=sc2;
	sc2=c;
	pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
	combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			if (pairrng[a][0]<max && pairrng[a][0]>=0)	{
				++combos[0][st1];
				++combos[1][st2];
				}
			}
		}
	}	

/*if (mxcm>1)	{	*/
st1=0;
for (c=0; c<sc1; ++c)	{
	clearlmatrix(pairs,states[ch1]+states[ch2],2,-1);
	while (combos[0][st1]<2)	++st1;
	if (st1>=states[ch1])	break;
	keyst1=st1;			/* the swing state is found!	*/
	/* find the other matches for keyst1	*/
	prs=0;		/* p is the number of pairs	*/
	clearlvector(gaps,(states[ch1]+states[ch2]),0);
	for (sp=0; sp<notu; ++sp)	{
		if (matrix[sp][ch1]==keyst1 && (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))	{
			b=0;
			/* make sure that this pair has not been found yet	*/
			for (d=0; d<=prs && b==0; ++d)	if (matrix[sp][ch2]==pairs[d][1])	b=1;	
			if (b==0)	{
				pairs[prs][0]=keyst1;					/* pair #p st. 1 found	*/
				pairs[prs][1]=matrix[sp][ch2];		/* pair #p st. 2 found	*/
				prfl[prs][0]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][0];	/* FA of pair p	*/
				prfl[prs][1]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][1];	/* LA of pair p	*/
				/* look for gaps in the range of pair prs	*/
				for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
					if (pairfnd[(keyst1*states[ch2])+matrix[sp][ch2]][f]==0)	++gaps[prs];
				++prs;	/* increment pairs	*/
				}
			}	/* end search for additional 2nd char state pairs	*/
		}

	keyst2=-1;
	for (st2=0; st2<states[ch2] && keyst2==-1; ++st2)	{
		if (combos[1][st2]>1)	{
			for (sp=0; sp<notu; ++sp)	{
				/* look at this only if this swing is paired with the swing from character 1: 
					this will not always be true for multistates; 2009.02.10						*/
				if (matrix[sp][ch2]==st2 && matrix[sp][ch1]==keyst1)	{
					keyst2=st2;
					sp=notu;		/* keystate 2 found, so kill the loop	*/
					}
				}

			if (keyst2>-1)	{
				rc=1;	/* remaining combos: when this hits combo[1][keyst2], then end	*/
				for (sp=0; sp<notu & rc<combos[1][keyst2]; ++sp)	{
					if (matrix[sp][ch2]==keyst2 && (matrix[sp][ch1]!=UNKNOWN && matrix[sp][ch1]!=INAP))	{
						b=0;
						/* make sure that we don't add pairs twice	*/ 
						for (d=0; d<prs && b==0; ++d)
							if (matrix[sp][ch1]==pairs[d][0])	b=1;

						if (b==0)	{
							pairs[prs][0]=matrix[sp][ch1];	/* pair #p st. 1 found	*/
							pairs[prs][1]=keyst2;				/* pair #p st. 2 found	*/
							prfl[prs][0]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][0];	/* FA of pair p	*/
							prfl[prs][1]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][1];	/* LA of pair p	*/
							/* look for gaps in the range of pair p	*/
							for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
								if (pairfnd[(pairs[prs][0]*states[ch2])+keyst2][f]==0)	++gaps[prs];
							++prs;	/* increment pairs	*/
							++rc;	/* one more of the "linked" pairs found	*/
							}	/* end pair	*/
						}	/* end case where we've found part of the pair	*/
					}	/*	end search for pairs	*/
				/* find cc, the "swing" pair; for 00, 10, 11, this would be 10: 							*/
				/*		it must be either intermediate (00->10->11 or 11->01->00) or ancestral (00<-10->11)	*/
				/*		thus, only one state pair can be older												*/
				/* ADD cl & cr: the "left" and "right" combos, to which this one will be compared; damn, I'm good..... */
				cc=-1;
				scf1=scf2=sw1=sw2=0;
				for (d=0; d<prs && cc==-1; ++d)	{
					if (pairs[d][0]==keyst1 && pairs[d][1]==keyst2)	cc=d;	/* "swing" pair	*/
					}
				/* now examine other state pairs with keystate 1	*/
				scu1=1;		/* this notes a gap between swing pairs and end pairs: it takes only one to set this to 0	*/
				for (d=0; d<prs; ++d)	{
					if (pairs[d][0]==keyst1 && d!=cc)	{
						cl[sw1]=d;											/* other pairs with key state 1	*/
						++sw1;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf1=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu1=0;
							/* d=p;	*/
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu1=0;
						}	/* end examination of "swing" pair involving keystate 1	*/
					}	/* end search for "swing" pairs involving keystate 1	*/
				/* now examine other state pairs with keystate 2	*/
				scu2=1;
				for (d=0; d<prs; ++d)	{
					if (pairs[d][1]==keyst2 && d!=cc)	{
						cr[sw2]=d;											/* other pairs with key state 2	*/
						++sw2;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf2=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu2=0;
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu2=0;
						}	/* end examination of "swing" pair involving keystate 2	*/
					}	/* end search for "swing" pairs involving key state 2	*/					
				}	/* end case where 2nd keystate is found	*/
				
			/* now, tally up stratigraphic compatibility	*/
			++strcmp[4];						/* another comparison made				*/
			if ((a=scf1+scf2)<2)	{
				++strcmp[5];									/* general stratigraphic compatibility			*/
				if (scu1==0 && scu2==0)			++strcmp[6];	/* range stratigraphic compatibility			*/
				if ((b=maxlarray(gaps,prs))==0)	++strcmp[7];	/* no gaps in the relevant character pairs		*/
				if (a==1)						++strcmp[8];	/* 10->00->01 instead of 10<-00->01				*/
				}
			}	/* end test of whether "swing state" is key state 2 	*/
		}	/* end search for 2nd key state	*/
	++st1;	/* if there are multiple swing states, then move on to the next possible one!	*/
	}	/* end survey of "swing" states matched	*/

free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
free_lmatrix(prfl,states[ch1]+states[ch2],2);
free_lmatrix(pairs,(states[ch1]+states[ch2]),2);
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
free_lvector(gaps);
return strcmp;
}

/** stratcompatfullplus - determines whether a compatible character pair are compatible with stratigraphy or not.  
/* Requires:
/*	ranges - a matrix of first and last appearances where:
/*		ranges[x][0]=first appearance;
/*		ranges[x][1]=last appearanc;
/*	matrix - cladistic character matrix;
/*	states - number of states for ach character;
/*	charcomps - compatibility of each character;
/*	ch1 - 1st character being compared;
/*	ch2 - 2nd character being compared;
/*	notu - number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	strcmp - an array telling you whether the character shows gaps in its  stratigraphic range or not.
		strcmp[0]: character 1
		strcmp[1]: character 2
		strcmp[2]: character 1 states
		strcmp[3]: character 2 states
		strcmp[4]: state comparisons (2 for binary; possibly more for multistate)
		strcmp[5]: first appearances compatible with stratigraphy
		strcmp[6]: pair synoptic ranges completely compatible with stratigraphy
		strcmp[7]: taxon synoptic ranges completely compatible with stratigraphy
		strcmp[8]: divergent (0) vs. hierarchical (1) stratigraphic compatibility
		strcmp[9]: hierarchical stratigraphic compatibility consistent with anagenesis
		strcmp[10]: hierarchical stratigraphic compatibility consistent with budding
		strcmp[11]: compatibility & stratigraphy suggests 00->11->10 transition rather than 10<-00->11
		strcmp[12]: compatibility & stratigraphy indecisive about 00->11->10 vs 10<-00->11
***************************************************************************************************************************/
unsigned long *stratcompatfullplus(long **ranges, long **matrix, int *states, unsigned long *charcomps, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int	a, b, c, d, f, g, h;
int	cc, sp, max, st1, st2, sc1, sc2, rc, scf1, scf2, scu1, scu2, prs, m1, m2;
int cl[20], cr[20];
int	keyst1, keyst2, sw1, sw2, swing;				/* these are the "swing" states between characters + important info, with "swing" giving pair #	*/
unsigned long *strcmp;
unsigned long **combos;
/*long **pairfnd;*/
long **pairs, **pairrng, **prfl, **pairfnd, *gaps;

max=maxlmatrixcol(ranges, notu, 1);

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/

strcmp=ulvector(13);
clearulvector(strcmp,13,0);

strcmp[0]=ch1;
strcmp[1]=ch2;
strcmp[2]=states[ch1];
strcmp[3]=states[ch2];

combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
for (st1=0; st1<states[ch1]; ++st1)	{
	for (st2=0; st2<states[ch2]; ++st2)	{
		a=st1*states[ch2]+st2;
		if (pairrng[a][0]<max && pairrng[a][0]>=0)	{
			++combos[0][st1];
			++combos[1][st2];
			}
		}
	}

sc1=sc2=0;
for (st1=0; st1<states[ch1]; ++st1)	if (combos[0][st1]>1)	++sc1;
for (st2=0; st2<states[ch2]; ++st2)	if (combos[1][st2]>1)	++sc2;

pairs=lmatrix(prs=(states[ch1]+states[ch2]),2);	/* state pairs between two characters		*/
prfl=lmatrix(prs,2);							/* first/last appearance of character pairs	*/
gaps=lvector(prs);

/*	Find each state paired with multiple states: these are the crux	*/
keyst1=-1;

/* play off of the character with the most swing states	*/
if (sc2>sc1)	{
	free_lmatrix(pairrng,states[ch1]*states[ch2],2);
	free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
	free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
	c=ch1;
	ch1=ch2;
	ch2=c;
	c=sc1;
	sc1=sc2;
	sc2=c;
	pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
	pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
	combos=ulmatrix(2,1+(a=imax(states[ch1],states[ch2])));
	for (st1=0; st1<states[ch1]; ++st1)	{
		for (st2=0; st2<states[ch2]; ++st2)	{
			a=st1*states[ch2]+st2;
			if (pairrng[a][0]<max && pairrng[a][0]>=0)	{
				++combos[0][st1];
				++combos[1][st2];
				}
			}
		}
	}	

/*if (mxcm>1)	{	*/
st1=0;
for (c=0; c<sc1; ++c)	{
	clearlmatrix(pairs,states[ch1]+states[ch2],2,-1);
	while (combos[0][st1]<2)	++st1;
	if (st1>=states[ch1])	break;
	keyst1=st1;			/* the swing state is found!	*/
	/* find the other matches for keyst1	*/
	prs=0;		/* p is the number of pairs	*/
	clearlvector(gaps,(states[ch1]+states[ch2]),0);
	for (sp=0; sp<notu; ++sp)	{
		if (matrix[sp][ch1]==keyst1 && (matrix[sp][ch2]!=UNKNOWN && matrix[sp][ch2]!=INAP))	{
			b=0;
			/* make sure that this pair has not been found yet	*/
			for (d=0; d<=prs && b==0; ++d)	if (matrix[sp][ch2]==pairs[d][1])	b=1;	
			if (b==0)	{
				pairs[prs][0]=keyst1;					/* pair #p st. 1 found	*/
				pairs[prs][1]=matrix[sp][ch2];		/* pair #p st. 2 found	*/
				prfl[prs][0]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][0];	/* FA of pair p	*/
				prfl[prs][1]=pairrng[(keyst1*states[ch2])+pairs[prs][1]][1];	/* LA of pair p	*/
				/* look for gaps in the range of pair prs	*/
				for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
					if (pairfnd[(keyst1*states[ch2])+matrix[sp][ch2]][f]==0)	++gaps[prs];
				++prs;	/* increment pairs	*/
				}
			}	/* end search for additional 2nd char state pairs	*/
		}

	keyst2=-1;
	for (st2=0; st2<states[ch2] && keyst2==-1; ++st2)	{
		if (combos[1][st2]>1)	{
			for (sp=0; sp<notu; ++sp)	{
				/* look at this only if this swing is paired with the swing from character 1: 
					this will not always be true for multistates; 2009.02.10						*/
				if (matrix[sp][ch2]==st2 && matrix[sp][ch1]==keyst1)	{
					keyst2=st2;
					sp=notu;		/* keystate 2 found, so kill the loop	*/
					}
				}

			if (keyst2>-1)	{
				rc=1;	/* remaining combos: when this hits combo[1][keyst2], then end	*/
				for (sp=0; sp<notu & rc<combos[1][keyst2]; ++sp)	{
					if (matrix[sp][ch2]==keyst2 && (matrix[sp][ch1]!=UNKNOWN && matrix[sp][ch1]!=INAP))	{
						b=0;
						/* make sure that we don't add pairs twice	*/ 
						for (d=0; d<prs && b==0; ++d)
							if (matrix[sp][ch1]==pairs[d][0])	b=1;

						if (b==0)	{
							pairs[prs][0]=matrix[sp][ch1];	/* pair #p st. 1 found	*/
							pairs[prs][1]=keyst2;				/* pair #p st. 2 found	*/
							prfl[prs][0]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][0];	/* FA of pair p	*/
							prfl[prs][1]=pairrng[(pairs[prs][0]*states[ch2])+keyst2][1];	/* LA of pair p	*/
							/* look for gaps in the range of pair p	*/
							for (f=prfl[prs][0]+1; f<prfl[prs][1]; ++f)
								if (pairfnd[(pairs[prs][0]*states[ch2])+keyst2][f]==0)	++gaps[prs];
							++prs;	/* increment pairs	*/
							++rc;	/* one more of the "linked" pairs found	*/
							}	/* end pair	*/
						}	/* end case where we've found part of the pair	*/
					}	/*	end search for pairs	*/
				/* find cc, the "swing" pair; for 00, 10, 11, this would be 10: 							*/
				/*		it must be either intermediate (00->10->11 or 11->01->00) or ancestral (00<-10->11)	*/
				/*		thus, only one state pair can be older												*/
				/* ADD cl & cr: the "left" and "right" combos, to which this one will be compared; damn, I'm good..... */
				cc=-1;
				scf1=scf2=sw1=sw2=0;
				for (d=0; d<prs && cc==-1; ++d)	{
					if (pairs[d][0]==keyst1 && pairs[d][1]==keyst2)	cc=d;	/* "swing" pair	*/
					}
				/* now examine other state pairs with keystate 1	*/
				scu1=1;		/* this notes a gap between swing pairs and end pairs: it takes only one to set this to 0	*/
				for (d=0; d<prs; ++d)	{
					if (pairs[d][0]==keyst1 && d!=cc)	{
						cl[sw1]=d;											/* other pairs with key state 1	*/
						++sw1;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf1=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu1=0;
							/* d=p;	*/
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu1=0;
						}	/* end examination of "swing" pair involving keystate 1	*/
					}	/* end search for "swing" pairs involving keystate 1	*/
				/* now examine other state pairs with keystate 2	*/
				scu2=1;
				for (d=0; d<prs; ++d)	{
					if (pairs[d][1]==keyst2 && d!=cc)	{
						cr[sw2]=d;											/* other pairs with key state 2	*/
						++sw2;
						/* if "end" pair appears before "swing" pair	*/
						if (prfl[d][0]<prfl[cc][0])	{
							scf2=1;						/* flag that end pair precedes swing pair	*/
							/* if there is no gap between the "end" pair's LA and the "swing" pair's FA	*/
							/* NOTE: in the case of multistates, only ONE of the "end" pairs need overlap with swing pair	*/
							if (prfl[d][1]>=(prfl[cc][0]-1))	scu2=0;
							}	/* end case where "end" pair appears before "swing" pair	*/
						/* if there is not gap between the "swing" pair's LA and the "end" pair's FA, then erase gap	*/
						else if (prfl[cc][1]>=(prfl[d][0]-1))	scu2=0;
						}	/* end examination of "swing" pair involving keystate 2	*/
					}	/* end search for "swing" pairs involving key state 2	*/					
				}	/* end case where 2nd keystate is found	*/
				
			/* now, tally up stratigraphic compatibility	*/
			++strcmp[4];						/* another comparison made				*/
			if ((a=scf1+scf2)<2)	{
				++strcmp[5];									/* general stratigraphic compatibility			*/
				if (scu1==0 && scu2==0)			++strcmp[6];	/* range stratigraphic compatibility			*/
				if ((b=maxlarray(gaps,prs))==0)	++strcmp[7];	/* no gaps in the relevant character pairs		*/
				if (a==1 && prs>=3)	{
					++strcmp[8];	/* 10->00->01 instead of 10<-00->01				*/
					/* 2011-05-16 Determine whether ranges overlap	*/
					/* go through relevant state pairs and determine whether LA[older]<=FA[younger]	*/
					/* this uses prfl for the 3 pairs in question: just need to figure out which 3 they are all of these years later....	*/
					/* OK, try this: 
							1) find swing pair based on keyst1 & keyst2.  IF we are here, then it is intermediate in age.
							2) Find the older pair with one of the swing states: did it disappear at same time or earlier?
								a) repeat if there are multiple pairs with the swing state
							3) Find the younger pair based on keyst1 & keyst2: did it appear at same time or later?	*/
					swing=-1;
					for (d=0; d<prs; ++d)	{
						if (pairs[d][0]==keyst1 && pairs[d][1]==keyst2)	{
							swing=d;
							d=prs;
							}
						}
					if (swing!=-1)	{
						g=h=0;
						for (d=0; d<prs; ++d)	{
							if (d!=swing && (pairs[d][0]==keyst1 || pairs[d][1]==keyst2))	{
								/* if state pair appears after or at same time swing pair disappears, then consistent with anagenesis	*/
								if (prfl[d][0]>=prfl[swing][1])	{
									if (g==0)	++strcmp[9];	/* consistent with anagenesis			*/
									else		++strcmp[10];	/* only one case of anagenesis possible	*/
									g=1;
									}
								/* if swing pair appears after or at same time state pair disappears, then consistent with anagenesis	*/
								else if (prfl[swing][0]>=prfl[d][1])	{
									if (h==0)	++strcmp[9];	/* consistent with anagenesis	*/	
									else		++strcmp[10];	/* only one case of anagenesis possible	*/
									h=1;
									}
								/* otherwise, we have overlapping ranges: consistent with budding	*/
								else			++strcmp[10];	/* consistent with budding		*/
								}	/* end comparisons with other state pairs	*/
							}	/* end search of other state pairs	*/
						}	/* double check: but there really has to be a swing pair in here!	*/
					}	/* end case of hierarchically compatible pair	*/
				}	/* end case of stratigraphically compatible characters	*/
			
			else	{
				/* 2011-05-16: Find swing pair and determine which character "reverse" to make this a literal interpretation of the fossil record	*/
				/* Try this: 
					1) find swing pair;
					2) find which of the two states is older	*/
				m1=m2=RAND_MAX;
				for (d=0; d<prs; ++d)	{
					if (pairs[d][0]==keyst1 && prfl[d][0]<m1)	m1=prfl[d][0];
					if (pairs[d][1]==keyst2 && prfl[d][0]<m2)	m2=prfl[d][0];
					}
				/* if the older state is from the less compatible character, then increment	*/
				if ((m1<m2 && charcomps[ch2]>charcomps[ch1]) || (m2<m1 && charcomps[ch1]>charcomps[ch2]))	++strcmp[11];
				else if (m2==m1 || charcomps[ch2]==charcomps[ch1])											++strcmp[12];
				}
			}	/* end test of whether "swing state" is key state 2 	*/
		}	/* end search for 2nd key state	*/
	++st1;	/* if there are multiple swing states, then move on to the next possible one!	*/
	}	/* end survey of "swing" states matched	*/

free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
free_lmatrix(prfl,states[ch1]+states[ch2],2);
free_lmatrix(pairs,(states[ch1]+states[ch2]),2);
free_ulmatrix(combos,2,1+(a=imax(states[ch1],states[ch2])));
free_lvector(gaps);
return strcmp;
}

/** statepairranges - determines the first and last appearance of each pair.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	pairrng - first and last appearance of each pair
****************************************************************************/
long **statepairranges(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int  a, st1, st2, sp, pairs;
long **pairrng;

pairs=1+(states[ch1]*states[ch2])+states[ch2];
pairrng=lmatrix(pairs,2);

for (a=0; a<pairs; ++a)	{
	pairrng[a][0]=RAND_MAX;
	pairrng[a][1]=-1*RAND_MAX;
	}

for (sp=0; sp<notu; ++sp)	{
	if (((st1=matrix[sp][ch1])!=UNKNOWN && (st2=matrix[sp][ch2])!=UNKNOWN) && (matrix[sp][ch1]!=INAP && matrix[sp][ch2]!=INAP))	{
		a=st1*states[ch2]+st2;
		if (ranges[sp][0]<pairrng[a][0])	pairrng[a][0]=ranges[sp][0];
		if (ranges[sp][1]>pairrng[a][1])	pairrng[a][1]=ranges[sp][1];
		}
	}
return	pairrng;
}

/** statepairfinds - determines whether each pair is found in each bin.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	pairfnd - whether each pair is found in each bin
****************************************************************************/
long **statepairfinds (long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int  a, st1, st2, sp, st, pairs;
long **pairfnd;

pairs=1+(states[ch1]*states[ch2])+states[ch2];
a=maxlmatrixcol(ranges,notu,1);
pairfnd=lmatrix(pairs,1+a);

for (sp=0; sp<notu; ++sp)	{
	if (((st1=matrix[sp][ch1])!=UNKNOWN && (st2=matrix[sp][ch2])!=UNKNOWN) && (matrix[sp][ch1]!=INAP && matrix[sp][ch2]!=INAP))	{
		a=(st1*states[ch2])+st2;			/* for 2 binary: 0<-00; 1<-01; 2<-10; 3<-11	*/
		for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
/*			if (pairfnd[a][st]==0)			pairfnd[a][st]=1;	*/
			pairfnd[a][st]=1;
		}
	}

return pairfnd;
}


/** stateranges - determines the first and last bin for each state.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	staterng - range of each state
****************************************************************************/
long **stateranges (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int  a, s, sp;
long **staterng;

staterng=lmatrix(states,2);

for (s=0; s<states; ++s)	{
	staterng[s][0]=RAND_MAX;
	staterng[s][1]=-1*RAND_MAX;
	}

for (sp=0; sp<notu; ++sp)	{
	if ((s=matrix[sp][ch])!=UNKNOWN && matrix[sp][ch]!=INAP)	{
		if (s<0)	{
			s/=-1;
			while (s>0)	{
				a=s%10;
				if (ranges[sp][0]<staterng[a][0])	staterng[a][0]=ranges[sp][0];
				if (ranges[sp][1]>staterng[a][1])	staterng[a][1]=ranges[sp][1];
				s/=10;
				}
			}
		else	{
			if (ranges[sp][0]<staterng[s][0])	staterng[s][0]=ranges[sp][0];
			if (ranges[sp][1]>staterng[s][1])	staterng[s][1]=ranges[sp][1];
			}
		}
	}
return	staterng;
}

/** statefinds - determines whether each state occurs in each bin.  
/* Requires:
/*	ranges - ranges of each taxon:
/*		range[x][y] = 1 if taxon x is found in bin y, 0 if absent;
/*	matrix - matrix of characters for notu taxa
/*	states - the number of states for character ch;
/*	ch - the character examined;
/*	notu - the number of taxa;
/*	UNKNOWN - coding for an unknown character;
/*	INAP - coding for an inapplicable character;
/*
/* Returns:
/*	stfnd - whether each state is found in each bin
****************************************************************************/
long **statefinds (long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int a, s, sp, st;
long **stfnd;

stfnd=lmatrix(states,1+(a=maxlmatrixcol(ranges,notu,1)));

for (sp=0; sp<notu; ++sp)	{
	if ((s=matrix[sp][ch])!=UNKNOWN && matrix[sp][ch]!=INAP)	{
		/* if polymorphic	*/
		if (s<0)	{
			s/=-1;
			while (s>0)	{
				a=s%10;
				for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
					if (stfnd[a][st]==0)			stfnd[a][st]=1;
				s/=10;
				}
			}
		/* if not polymorphic	*/
		else	{
			for (st=ranges[sp][0]; st<=ranges[sp][1]; ++st)
				if (stfnd[s][st]==0)			stfnd[s][st]=1;
			}
		}
	}

return stfnd;
}


/** stratconsiststate - determines whether all state in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	stfnd - whether a state is found in each bin
/*		stfnd[x][y] = 1 if state x is found in bin y, 0 if absent;
/*	s - the state being examined;
/*
/* Returns:
/*	scs - 0 if not continuous, 1 if continuous
****************************************************************************/
unsigned long stratconsiststate(long **staterng, long **stfnd, int s)
{
int a, st;
unsigned long scs=1;

for (st=staterng[s][0]+1; st<staterng[s][1] && scs==1; ++st)	{
	a=stfnd[s][st];
	if (stfnd[s][st]==0)	scs=0;
	}
return scs;
}


/** stratconsistchar - determines whether all state in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	states - the number fo states for the character;
/*
/* Returns:
/*	scc - 0 if not all are continuous, 1 if all are continuous
****************************************************************************/
unsigned long stratconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int	a, s;
unsigned long scc=1;
long **staterng, **stfnd;

stfnd=statefinds(ranges, matrix, states, ch, notu, UNKNOWN, INAP);		/* pres/abs of states within strat bins	*/
staterng=stateranges(ranges, matrix, states, ch, notu, UNKNOWN, INAP);	/* find ranges of individual states	*/

for (s=0; s<states && scc==1; ++s)	{
	a=stratconsiststate(staterng, stfnd, s);
	if (a==0)	scc=0;
	}
free_lmatrix(staterng,states,2);
free_lmatrix(stfnd,states,2);
return scc;
}


/** stratconsiststpair - determines whether a state pair for two characters has a continuous ranges.  
/* Requires:
/*	pairrng - ranges of each pair:
/*		pairrng[x][0]=first appearance of pair x;
/*		pairrng[x][1]=last appearance of pair x;
/*	pr - the pair number, given by (states_char_1 x states_char_2) + state_of_char_2;
/*
/* Returns:
/*	scsp - 0 if incompatible, 1 if compatible
****************************************************************************/
unsigned long stratconsiststpair(long **pairrng, long **pairfnd, int pr)
{
int st;
unsigned long scsp=1;

/* determine whether state pairs have consistent ranges	*/
for (st=pairrng[pr][0]+1; st<pairrng[pr][1]; ++st)	{
	if (pairfnd[pr][st]==0)	scsp=0;
	} 
return scsp;
}

/** stratconsistchpair - determines whether all state pairs for two characters have continuous ranges.  
/* Requires:
/*	pairrng - ranges of each pair:
/*		pairrng[x][0]=first appearance of pair x;
/*		pairrng[x][1]=last appearance of pair x;
/*	st1 - states for 1st character;
/*	st2 - states for 2nd character;
/*	max - oldest stratigraphic bin;
/*
/* Returns:
/*	scp - 0 if incompatible, 1 if compatible
****************************************************************************/
unsigned long stratconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int pr, cs1, cs2;
unsigned long scp=1;
long max;
long **pairrng, **pairfnd;

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
max=maxlmatrixcol(ranges, notu, 1);

for (cs1=0; cs1<states[ch1] && scp==1; ++cs1)	{
	for (cs2=0; cs2<states[ch2] && scp==1; ++cs2)	{
		pr=(cs1*states[ch2])+cs2;
		if (pairrng[pr][0]>=0 && pairrng[pr][0]<=max)	
			scp=stratconsiststpair(pairrng, pairfnd,pr);
		}
	}
free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
return scp;
}


/** stratfrconsistchpair - determines the proportion of state pairs for two characters have continuous ranges.  
/* Requires:
/*	pairrng - ranges of each pair:
/*		pairrng[x][0]=first appearance of pair x;
/*		pairrng[x][1]=last appearance of pair x;
/*	st1 - states for 1st character;
/*	st2 - states for 2nd character;
/*	max - oldest stratigraphic bin;
/*
/* Returns:
/*	scp - 0 if incompatible, 1 if compatible
****************************************************************************/
double stratfrconsistchpair(long **ranges, long **matrix, int *states, int ch1, int ch2, int notu, int UNKNOWN, int INAP)
{
int a=0, pr, cs1, cs2;
long max;
long **pairrng, **pairfnd;
double y=0.0f, scp=0.0f;

pairrng=statepairranges(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get ranges for state pairs	*/
pairfnd=statepairfinds(ranges, matrix, states, ch1, ch2, notu, UNKNOWN, INAP);	/* get pre/abs for state pairs	*/
max=maxlmatrixcol(ranges, notu, 1);

for (cs1=0; cs1<states[ch1]; ++cs1)	{
	for (cs2=0; cs2<states[ch2]; ++cs2)	{
		pr=(cs1*states[ch2])+cs2;
		if (pairrng[pr][0]>=0 && pairrng[pr][0]<=max)	{
			a=a+stratconsiststpair(pairrng, pairfnd,pr);
			++y;
			}
		}
	}
scp=((double) a)/y;

free_lmatrix(pairrng,states[ch1]*states[ch2],2);
free_lmatrix(pairfnd,states[ch1]*states[ch2],1+max);
return scp;
}

/** stratconsistchar - determines the proportion of  all states in a character have continuous ranges.  
/* Requires:
/*	staterng - ranges of each state:
/*		staterng[x][0]=first appearance of state x;
/*		staterng[x][1]=last appearance of state x;
/*	states - the number fo states for the character;
/*
/* Returns:
/*	scc - 0 if not all are continuous, 1 if all are continuous
****************************************************************************/
double stratfrconsistchar(long **ranges, long **matrix, int states, int ch, int notu, int UNKNOWN, int INAP)
{
int	a=0, s;
double scc=0.0f;
long **staterng, **stfnd;

stfnd=statefinds(ranges, matrix, states, ch, notu, UNKNOWN, INAP);		/* pres/abs of states within strat bins	*/
staterng=stateranges(ranges, matrix, states, ch, notu, UNKNOWN, INAP);	/* find ranges of individual states	*/

for (s=0; s<states; ++s)	{
	a=a+stratconsiststate(staterng, stfnd, s);
	scc+=a;
	}

scc=((double) a)/((double) states);
free_lmatrix(staterng,states,2);
free_lmatrix(stfnd,states,2);
return scc;
}

/* cleancladerangedata: remove gaps from a clade's range
****************************************************************************/
void cleancladerangedata(long **ranges, int notu)
{
int	b, c, sp;
int	lfa, lla;
long	*rich;

lfa=maxlmatrixcol(ranges,notu,0);
lla=maxlmatrixcol(ranges,notu,1);
rich=richnesstally(ranges,notu);		/* Sepkoski style synoptic richness per bin	*/

for (c=0; c<=lla; ++c)	{
	if (rich[c]==0)	{
		for (sp=0; sp<notu; ++sp)	{
			if (ranges[sp][0]>c)	{
				--ranges[sp][0];
				--ranges[sp][1];
				}
			else if (ranges[sp][1]>c)
				--ranges[sp][1];
			}
		for (b=c; b<lla; ++b)	rich[b]=rich[b+1];
		--c;
		--lla;						/* decrement the final last appearance								*/
		if (c<lfa)	--lfa;			/* decrement the final first appearance only if it is after bin c	*/
		}	/* end case of gap in clade range */
	}	/* end search for gaps in clade range	*/
free_lvector(rich);
/*return (0);*/
}

/* nu_stratcompat: number of stratigraphic compatibilities for each character
Requires:
	charmat: character matrix;
	compmat: compatibilty matrix;
	ranges: stratigraphic ranges of taxa;
	nstats: number of states for each character
	notu: number of taxa
	strong: whether stratigraphic compatibilty is "weak" (0=based only on first appearances) or "strong" (1=using whole ranges)
	hier: whether we count all (0=divergent+hierarchical) or only hierarchical (1=00->01->11) stratigraphic compatibilities
	UNKNOWNN: code for unknown states
	INAP: code for inapplicable state
****************************************************************************
unsigned long *nu_stratcompat(unsigned long **charmat, unsigned long **compmat, long **ranges, int nchars, int notu, int strong, int hier, int UNKNOWN, int INAP)
{
int c1, c2;


for (c1=0; c1<(nchars-1); ++c1)	{
	for (c2=(c1+1); c2<nchars; ++c2)	{
		}
	}
}	*/


/* stratcompatparambstest: simulations for stratigraphic compatibility with parameters determined by best fit for stratigraphic data
Requires:
	taxonname: group analyzed
	citation: source of data
	summary: summary of observed data
	mbl: simulation parameters
	omatrix: observed matrix
	ctype: character types
	nstates: states per character
	bias: character evolution bias
	maxch: maximum change per character
	depend: character dependencies
	notu: number of taxa
	nchars: number of characters
	compat: compatibility type (general or hierarchical)
	RUNS: simulations
	debug: 1 to seed based on replications, 0 for time
	UNKNOWN: code for unknown states
	INAP: code for inapplicable states
Returns:
	simmary[0]: p observed general stratigraphic compatibility
	simmary[1]: p observed strict stratigraphic compatibility
	simmary[2]: p observed hierarchical stratigraphic compatibility
	simmary[3]: expected general stratigraphic compatibility
	simmary[4]: expected strict stratigraphic compatibility
	simmary[5]: expected hierarchical stratigraphic compatibility
	simmary[6]: p proportional hierarchical stratigraphic compatibility
	simmary[7]: expected proportional hierarchical stratigraphic compatibility
	simmary[8]: p hierarchical compatibility consistent with anagenesis
	simmary[9]:	expected hierarchical compatibility consistent with anagenesis
	simmary[10]: p observed �reversal� character being less compatible
	simmary[11]: expected number of stratigraphically incompatible pairs in which �reversed� character is less compatible
	simmary[12]: p observed average compatibility of stratigraphically compatible pairs
	simmary[13]: expected compatibility of stratigraphically compatible pairs
	simmary[14]: p observed average compatibility of stratigraphically incompatible pairs
	simmary[15]: expected compatibility of stratigraphically incompatible pairs
NOTE: Should I add expected compatibility of SIN & GSC characters?  Ummm.... sure, why not?
****************************************************************************/
double *stratcompatparambstest(char taxonname[90], char citation[90], double *summary, double *mbl/*, long **omatrix*/, int *ctype, int *nstates, int *bias, int *maxch, int *depend, int notu, int nchars, int compat, int RUNS, int debug, int UNKNOWN, int INAP)
{
/* part one - matrix properties */
int		r, s;
int		c1, c2;
int		keep;
int 	clades;
//int	uninf, out=0, uninfcompat, posscompat, auts, ttlcompat, comptype=0;
long  	**chmatrix, **ranges;
long	**tree, **vtree;
//long	**dbtree;
unsigned long	*sstrcmp;
unsigned long **compmat;
unsigned long	*charcomps;
double	*genstratcomp, *strstratcomp, *hierstratcmp, *hierstratprp,*scha,*sih;
double	*avecompscp, *avecompsin;
double	*simmary;
double	cp;
double	tiegs=0.0f, tiess=0.0f, tiehs=0.0f, tiehp=0.0f, tiescha=0.0f, tiesih=0.0f, tieavecscp=0.0f, tieavecsin=0.0f;
double	sci, scgw;
//double	**genrstcmpcor, **strcstcmpcor, **hierstcmpcor;
//long	**histories;
long	secs;
char	outfile[120];
FILE	*output;

if (RUNS%2==0)	++RUNS;			/* it's easier to use an odd number....	*/

time(&secs);
srand((unsigned int) secs);

genstratcomp=dvector(RUNS);
strstratcomp=dvector(RUNS);
hierstratcmp=dvector(RUNS);
hierstratprp=dvector(RUNS);
scha=dvector(RUNS);
sih=dvector(RUNS);
avecompscp=dvector(RUNS);
avecompsin=dvector(RUNS);
ranges=lmatrix(notu,2);
chmatrix=lmatrix(notu,nchars);
simmary=dvector(16);

for (r=0; r<RUNS; ++r)	{
	sci=scgw=0.0f;
	/* 2011-03-28: BLOWOUT at r=5	*/
	genstratcomp[r]=strstratcomp[r]=hierstratcmp[r]=hierstratprp[r]=0.0f;
	if (debug==1)	{
		srand((unsigned int) (r+1)*notu*nchars);
		}
	
	/* 	trees[notu-1]: branch length of otus	
		trees[notu]:   branch length of clades
		trees[notu+1]: first appearances of otus
		trees[notu+2]: last appearances of otus
	**********************************************/
	tree=evolvetree(notu,mbl,1);
	clades=cladecountbytaxa(tree,notu);
	/* pull range data out of back of tree matrix	*/
	for (s=0; s<notu; ++s)	{
		ranges[s][0]=tree[notu+1][s];
		ranges[s][1]=tree[notu+2][s];
		}
	cleancladerangedata(ranges,notu);

	vtree=VennTreePlus(tree,clades,notu,notu);
	/* evolve character matrix, chmatrix	*/
//	long **evolvetocompat(long **tree, int tcomp, int notu, long **matrix, int nchars, int *nstates, int *ctype, int *bias, int *maxch, int *depend, int comptype, int UNKNOWN, int INAP)
//	chmatrix=evolvetocompat(tree,empcompat,notu,matrix,nchars,nstates,ctype,bias,chmax,depend,comptype,UNKNOWN,INAP);
	chmatrix=evolvetocompatibility(vtree, compat, notu, clades, chmatrix, nchars, nstates, ctype, bias, maxch, depend, 0, UNKNOWN, INAP);
//	chmatrix=evolvematrix(vtree, notu, clades, chmatrix, nchars, nstates, ctype, bias, maxch, 2*nchars, UNKNOWN, INAP);
/*		sstates=numberstates(chmatrix,notu,nchars,UNKNOWN,INAP);	*/
//	compat=nu_comp(nstates, notu, chmatrix, ctype, nchars, 0, 0, UNKNOWN, INAP);
	compmat=compatible(nstates,notu,chmatrix,ctype,nchars,0,0,UNKNOWN,INAP);
	charcomps=char_comp(nstates,notu,chmatrix,ctype,nchars,0,0,UNKNOWN,INAP);

/*	for (c1=0; c1<nchars; ++c1)	++examples[compmat[c1][c1]];	*/
	cp=0.0f;
	for (c1=0; c1<(nchars-1); ++c1)	{
		for (c2=(c1+1); c2<nchars; ++c2)	{
			if (compmat[c1][c2]==1)	{
				/*sstrcmp=stratcompatfull(ranges, chmatrix, nstates, c1, c2, notu, UNKNOWN, INAP);	*/
				/* 2011-10-14: changed to "plus"													*/
				sstrcmp=stratcompatfullplus(ranges, chmatrix, nstates, charcomps, c1, c2, notu, UNKNOWN, INAP);
				cp+=sstrcmp[4];								/* compatible comparisons					*/
				genstratcomp[r]+=sstrcmp[5];				/* general stratigraphic consistency		*/
				strstratcomp[r]+=sstrcmp[6];				/* strict stratigraphic consistency			*/
				hierstratcmp[r]+=sstrcmp[8];				/* hierarchical stratigraphic consistency	*/
				scha[r]+=sstrcmp[9];						/* hierarchical compatible pairs consistent with anagenesis			2011-10-14	*/
				sih[r]+=sstrcmp[11];						/* stratigraphically incompatible pairs consistent with hierarchy	2011-10-14	*/
				/* 2011-10-14: tally compatibilities of stratigraphically compatible & incompatible character pairs	*/
				if (sstrcmp[5]==1)	{
					avecompscp[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					scgw+=1.0;													/* this will tally characters, not state pairs						*/
					}
				else	{
					avecompsin[r]+=((double) (charcomps[c1]+charcomps[c2]))/2;
					sci+=1.0;
					}
	
				/* this will give the proportion of compatibilities that are stratigraphically compatible by compatibility	*/
/*				genrstcmpcor[compmat[c1][c1]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c1][c1]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c1][c1]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
/*				genrstcmpcor[compmat[c2][c2]]+=sstrcmp[5];			/* general stratigraphic consistency		*/
/*				strcstcmpcor[compmat[c2][c2]]+=sstrcmp[6];			/* strict stratigraphic consistency			*/
/*				hierstcmpcor[compmat[c2][c2]]+=sstrcmp[8];			/* hierarchical stratigraphic consistency	*/
				free_ulvector(sstrcmp);
				}
			}
		}
	/* cp will vary a little from simulation to simulation because it is based on state numbers, so report things in proportions	*/
	genstratcomp[r]/=((double) cp);
	strstratcomp[r]/=((double) cp);
	hierstratcmp[r]/=((double) cp);
	hierstratprp[r]=hierstratcmp[r]/genstratcomp[r];
	scha[r]/=((double) cp);																	/* added 2011-10-14	*/
	sih[r]/=((double) cp);																	/* added 2011-10-14	*/
	avecompscp[r]/=scgw;																	/* added 2011-10-14	*/
	avecompsin[r]/=sci;																		/* added 2011-10-14	*/

	/* tally cases where simulated values are lower than observed values		*/
	if (genstratcomp[r]<summary[0])		simmary[0]+=(1/((double) RUNS));
	if (strstratcomp[r]<summary[1])		simmary[1]+=(1/((double) RUNS));
	if (hierstratcmp[r]<summary[2])		simmary[2]+=(1/((double) RUNS));
	if (hierstratprp[r]<(summary[2]/summary[0]))
										simmary[6]+=(1/((double) RUNS));
	if (scha[r]<summary[3])				simmary[8]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (sih[r]<summary[4])				simmary[10]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompscp[r]<summary[5])		simmary[12]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]<summary[6])		simmary[14]+=(1/((double) RUNS));					/* added 2011-10-14	*/
	/* tally cases where simulated values are exactly equal to observed values	*/
	if (genstratcomp[r]==summary[0])	tiegs+=(1/((double) RUNS));
	if (strstratcomp[r]==summary[1])	tiess+=(1/((double) RUNS));
	if (hierstratcmp[r]==summary[2])	tiehs+=(1/((double) RUNS));
	if (hierstratprp[r]==(summary[2]/summary[0]))
										tiehp+=(1/((double) RUNS));
	if (scha[r]==summary[3])			tiescha+=(1/((double) RUNS));						/* added 2011-10-14	*/
	if (sih[r]==summary[4])				tiesih+=(1/((double) RUNS));						/* added 2011-10-14	*/
	if (avecompscp[r]==summary[5])		tieavecscp+=(1/((double) RUNS));					/* added 2011-10-14	*/
	if (avecompsin[r]==summary[6])		tieavecsin+=(1/((double) RUNS));					/* added 2011-10-14	*/

	if ((r%10)==9)	{
//			printf("Doing Rate = %5.4f",mbl[2]);
		if (r>10)	{
			printf("\b\b\b\b\b");					/* clear ", Tree %d"		*/
			}
		if (r>18 && r<100) 			printf("\b\b\b");
		else if (r>99 && r<1000)	printf("\b\b\b\b");
		else if (r>999)				printf("\b\b\b\b\b");
		printf("Tree %d\n", r+1);
		}
	
	free_ulmatrix(compmat,nchars,nchars);
	free_lmatrix(vtree,clades+1,notu);
	free_lmatrix(tree,clades+3,notu);
	}

/* sort for medians	*/
genstratcomp=dshellsort_inc(genstratcomp,RUNS);
strstratcomp=dshellsort_inc(strstratcomp,RUNS);
hierstratcmp=dshellsort_inc(hierstratcmp,RUNS);
hierstratprp=dshellsort_inc(hierstratprp,RUNS);
scha=dshellsort_inc(scha,RUNS);						/* added 2011-10-14	*/
sih=dshellsort_inc(sih,RUNS);						/* added 2011-10-14	*/
avecompscp=dshellsort_inc(avecompscp,RUNS);			/* added 2011-10-14	*/
avecompsin=dshellsort_inc(avecompsin,RUNS);			/* added 2011-10-14	*/

/* get alpha values: if observed values are "low" then we want p[observed or lower]: so, add ties to total	*/
/* for "high" values, we do not need to do this: 1-simmary gives p[observed or more extreme] as we tallied only sims less than observed	*/
if (simmary[0]<0.5)		simmary[0]+=tiegs;
if (simmary[1]<0.5)		simmary[1]+=tiess;
if (simmary[2]<0.5)		simmary[2]+=tiehs;
if (simmary[6]<0.5)		simmary[6]+=tiehp;
if (simmary[8]<0.5)		simmary[8]+=tiescha;
if (simmary[10]<0.5)	simmary[10]+=tiesih;
if (simmary[12]<0.5)	simmary[12]+=tieavecscp;
if (simmary[14]<0.5)	simmary[14]+=tieavecscp;
/* tally median expected stratigraphic compatibilities	*/
simmary[3]=genstratcomp[RUNS/2];
simmary[4]=strstratcomp[RUNS/2];
simmary[5]=hierstratcmp[RUNS/2];
simmary[7]=hierstratprp[RUNS/2];
simmary[9]=scha[RUNS/2];
simmary[11]=sih[RUNS/2];
simmary[13]=avecompscp[RUNS/2];
simmary[15]=avecompsin[RUNS/2];

keep=1;
//keep=0;
//printf("Enter '1' if you want output for this taxon: ");
//scanf("%i",&keep);	

if (keep==1)	{
	strcpy(outfile,taxonname);
	strcat(outfile,"_");
	strcat(outfile,citation);
	strcat(outfile,"_expectations.xls");
	
	output=fopen(outfile,"w");
	fprintf(output,"Run\tGeneral Stratigraphic Compatibility\tStrict Stratigraphic Compatibility\tHierarchical Stratigraphic Compatibility\tprop. Hierarchical Stratigraphic Compatibility\tPoss Anagenetic Hierarchical Stratigraphic Compatibility\n");
	for (r=0; r<RUNS; ++r)
		fprintf(output,"%d\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n",r+1,genstratcomp[r],strstratcomp[r],hierstratcmp[r],hierstratprp[r],scha[r],sih[r],avecompscp[r],avecompsin[r]);
	fclose(output);
	}

free_dvector(genstratcomp);
free_dvector(strstratcomp);
free_dvector(hierstratcmp);
free_dvector(hierstratprp);
free_dvector(scha);
free_dvector(sih);
free_dvector(avecompscp);
free_dvector(avecompsin);
free_lmatrix(ranges,notu,2);

return simmary;
}