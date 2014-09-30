/* 
 * File:   main.cpp
 * Author: Randy
 *
 * Created on September 30, 2014, 10:22 AM
 */

#include "PaleoMCMC.h"
#include <iostream>


#define UNKNOWN			-88		/* Code for unknown characters	*/
#define	INAP			-99		/* Code for inapplicables		*/
#define TREE_END		-1
#define	OUTGROUP		0		/* Number of outgroup taxon 	*/
#define SW_TREE			0
#define SW_BRANCH		1
#define SW_ANCESTOR		2
#define SW_ALPHA		3
#define SW_BETA			4
#define SW_POLYT		5

/***************************************************************************************************
chmatrix[i][j]		: states for taxon i character j.
nstates[i]   		: real number of states for character i
ctype[i]      		: ordered [0], unordered [1];
comatrix[i][j]		: compatibility of character i and character j.
*****************************************************************************************************/

/* Priorities:
	1. Establish loop of branch swapping (done 2014-09-12)
	2. Put branch length modification in (assuming exponential priors)  (done 2014-09-15)
	2. Put in initial starting estimate (done 2014-09-15)
	3. Put Mk into branch swapping to execute pully principle (mostly done)
	4. Put anagenetic ancestry into continuous change model & estimate likelihoods given stasis + tree
	5. Put stratigraphic likelihood into swapping
	6. Put geographic stratigraphic likelihood into swappng
	7. Add punctuated change option for swapping	(partly done: need to add polytomies)
	8. Add pulsed speciation change option for swapping
	9. Get starter tree assuming one of the three models
*****************************************************************************************************/

/* notes:
	A. Use oldest divergence for both pulsed and continuous change: that makes P[X species | time, lambda] always included	
	B. set up continuous (alpha) and pulsed (beta) rates all the time; set beta=0 for purely continuous change	
	C. set up lambda[0] for continuous speciation & lambda[1] for pulsed speciation rates
		then set up array giving pulses
	D. To do ancestor-descendant hypotheses, make ancestral species the ancestor in the ape vector;
		However, assign likelihoods to node including ancestral species.  (That is P[anc->desc states])
		To do this:
			just reset marginal likelihoods of node to those for ancestral species.
			then modify the likelihood of the observed states.  (If Unknown, then it's like a normal node).
	E. On anagenetic descendants with contemporaneous ancestors (implying cryptic species), do not use entire divergence time for
		stratigraphic likelihood; instead, base it on the minimum number of changes so that if there is:
			1 change: half the branch;
			2 changes: 2/3rds the branch;
			3 chanegs: 3/4ths the branch; etc.
		NOTE: if anagenetic descendant comes AFTER ancestor disappears, then use full divergence time
*****************************************************************************************************/

/* MCMC options
	  I. swap branches (use pulley principle)
	 II. alter branch lengths (use pulley principle; need ranges & divergences to make sure that we don't make divergences too late!
	III. make ancestors (use pulley principle; need ranges; need anagenetic descendant vector
	 IV. alter alpha (recalculate likelihood of whole tree)
	  V. alter beta (recalculate likelihood of whole tree)
*****************************************************************************************************/



#include <cstdlib>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
        int		a, c, i, r, s, t;
int		ch;
//double	*margin, *marginch, **margins;
int		debug=1;
//long	secs;
int *info;								/* dummy for getting matrix dimensions	*/
int	SWITCH;
double	x, y;
char	input1[200], input2[200];
char	inapcode;
char	polytomies='y';
char	output1[90];
FILE	*outfile1;
int		runs=100, rates;
long	**tree;
long	**omatrix, **chmatrix;
int		*ctype, *nstates, *novels;
long	*charswstates, **statechars;
int	notu, ochars, nchars;
//char	***instructions, commands[120], taxonname[120];
long	*ape, *last_ape;
long	apeT[100];
unsigned long	*minbr, *ancestral, *last_minbr, *last_ancestral;
int	nodes, ttlbr, maxst, last_nodes;
double	xx[100], yy[100];
double	*alphas, *betas;	/* continuous change, pulsed change, geographic change					*/
double	*last_alphas, *last_betas;	/* continuous change, pulsed change, geographic change					*/
double	*lnl_line_stasis, *last_lnl_line_stasis;	/* likelihood of stasis within observed lineages						*/
double	**ranges;	/* ranges are observed; divergences are part of the tree				*/
double	**divergences, **last_divergences;
double	***statelikes, ***last_statelikes;				/* conditonal likelihoods of each state for each character at each node	*/
double	lambda=0.585f, theta=0.1;
double	ma_otu=0.0f, ma_htu=0.0f;
double	lnL_tree, lnL_stasis, last_lnL_tree, last_lnL_stasis, last_lnL;
double	met_hast;
//double	**Gibbs;
//int		nst[15];


cout << "\n";
    cout << "Hello world ... CPP Version!\n";

        //x=pow(e,(-1*2*(0.01*1)+(0.05*1.05)));
    //y=exp(-1*2*(0.01*1)+(0.05*1.05));

    inapcode='-';

    //strcpy(input1,"/Users/peterjwagner/Documents/Projects/PaleoMCMC/DummyMatrix.txt");
    strcpy(input1,"DummyMatrix.txt");
    info=getmatrixinfo(input1);
    notu=info[1];
    
    ochars=info[0];		/* original number of characters	*/
    /* reach matrix	*/
    omatrix=readclmatrix(input1,notu,ochars,inapcode,UNKNOWN,INAP);	/* original matrix						*/
    nchars=countuniquechars(omatrix,notu,ochars);					/* unique character strings				*/
    novels=countcharreplicates(omatrix,notu,ochars,nchars);			/* number of characters sharing string	*/
    chmatrix=reducecladisticmatrix(omatrix,notu,ochars);			/* reduced character matrix				*/

    ctype=ivector(nchars);
    nstates=numberstates(chmatrix, notu, nchars, UNKNOWN, INAP);
    maxst=1+maxiarray(nstates, nchars);
    charswstates=lvector(maxst);
    a=0;
    
    nodes=notu-1;
    for (ch=0; ch<nchars; ++ch)	{
	if (nstates[ch]<2)	nstates[ch]=2;
	++charswstates[nstates[ch]];
	if (charswstates[nstates[ch]]>a)	
            a=charswstates[nstates[ch]];
	}
    ++a;
    free_lvector(charswstates);

    charswstates=lvector(maxst);
    statechars=lmatrix(maxst,a);
    for (ch=0; ch<nchars; ++ch)	{
	statechars[nstates[ch]][charswstates[nstates[ch]]]=ch;
	++charswstates[nstates[ch]];
	}

    ancestral=ulvector(notu+nodes);
    for (c=notu; c<(notu+nodes); ++c)	
        ancestral[c]=2;	/* nodes have 2 descendants	*/
    minbr=ulvector((2*notu)-1);
    clearulvector(minbr,(2*notu)-1,1);	/* set minimum branchings to 1 initially	*/
    //for (ch=0; ch<nchars; ++ch)	nst[ch]=nstates[ch];

    ttlbr=(2*notu)-1;
    //strcpy(input2,"/Users/peterjwagner/Documents/Projects/PaleoMCMC/DummyRanges.txt");
    strcpy(input2,"DummyRanges.txt");
	
    ranges=readdmatrix(input2, notu, 2);
    divergences=dmatrix(ttlbr,2);
    tree=lmatrix(ttlbr,notu+1);
    ape=lvector(ttlbr);
    //for (i=0; i<11; ++i)	tree[i][0]=i;
    apeT[0]=notu+1;
    apeT[1]=notu+1;
    apeT[2]=notu+3;
    apeT[3]=notu+3;
    apeT[4]=notu+4;
    apeT[5]=notu+5;
    apeT[6]=notu+5;
    apeT[notu+0]=-1;		// A
    apeT[notu+1]=notu;	// B
    apeT[notu+2]=notu;	// C
    apeT[notu+3]=notu+2;	// D
    apeT[notu+4]=notu+2;	// E
    apeT[notu+5]=notu+4;	// F
    ttlbr=nodes+notu;

    ape=apeT;
    convertApeToNexus(ape,tree,notu,nodes);
    datecladerealape_addexp(ape, ranges, divergences, nodes, notu, 2*lambda, ancestral);

    ma_htu=ma_otu=0.0f;
    for (s=0; s<notu; ++s)	{
	ma_otu+=(xx[s]=ranges[s][1]-ranges[s][0]);
	ma_htu+=(yy[s]=divergences[s][1]-divergences[s][0]);
	}
    for (c=(notu+1); c<(notu+nodes); ++c)	{
	ma_htu+=(yy[c]=divergences[c][1]-divergences[c][0]);
	}
    ma_htu=sum_divergence_times(divergences,notu,nodes);

    alphas=dvector(rates=2);
    for (i=0; i<rates; ++i)		
        alphas[i]=(18/((double) ochars))/(ma_htu+ma_otu);	/* 18 is changes in DummyMatrix.txt	*/
    betas=dvector(rates=2);
    for (i=0; i<rates; ++i)		
        betas[i]=(0/((double) ochars))/(ma_htu+ma_otu);	/* let's assume no pulsed change at first	*/
    statelikes=initialize_conditional_likelihoods(chmatrix, nstates, notu, (notu-1), nchars, maxst, UNKNOWN, INAP);
    //statelikes=dcube(ttlbr,nchars,maxst);	/* likelihoods of each state at each node	*/

    lnl_line_stasis=dvector(notu);
    //stasis_likelihood_simple(notu, chmatrix, ranges, lnl_line_stasis, statechars, charswstates, novels, int maxst, double *alphas, int rates, int UNKNOWN, int INAP);
    stasis_likelihood_simple(notu, chmatrix, ranges, lnl_line_stasis, statechars, charswstates, novels, maxst, alphas, rates, UNKNOWN, INAP);
    lnL_stasis=sumdvector(lnl_line_stasis,notu);

    MkSimpleWholeTree(ape, notu, nodes, divergences, statechars, charswstates, maxst, alphas, betas, rates, statelikes, lambda, ancestral, minbr);
    lnL_tree=get_tree_likelihood(statelikes[notu], novels, nstates, nchars);

    /* START HERE 2014-09-17	*/
    t=0;
    printf("Tree #%d: lnL=%4.3f; · divergences = %4.1f; ",++t,(lnL_tree+lnL_stasis), ma_htu);
    write_nexus(tree,0,notu);

    //outfile1=fopen("/Users/peterjwagner/Documents/Projects/PaleoMCMC/Ape_Test_Autumn.xls","w");
    outfile1=fopen("Ape_Test_Autumn.xls","w");
    fprintf(outfile1,"Swap");
    for (i=0; i<ttlbr; ++i)	{
	if (i<10)	
            fprintf(outfile1,"\tSp0%d",i);
	else		
            fprintf(outfile1,"\tSp%d",i);
	}
    fprintf(outfile1,"\n");
    for (i=0; i<ttlbr; ++i)	
        fprintf(outfile1,"\t%d",ape[i]);
    fprintf(outfile1,"\n");
    fclose(outfile1);

    last_lnL=lnL_tree+lnL_stasis;
    last_lnL_tree=lnL_tree;
    last_lnL_stasis=lnL_stasis;
    last_divergences=dmatrix((2*notu)-1,2);
    equaldmatrix(last_divergences,divergences,(2*notu)-1,2);
    last_ape=lvector((2*notu)-1);
    equallvector(last_ape,ape,(2*notu)-1);
    last_statelikes=dcube((2*notu)-1,nchars,maxst);
    equaldcube(last_statelikes,statelikes,(2*notu)-1,nchars,maxst);
    last_lnl_line_stasis=dvector(notu);
    equaldvector(last_lnl_line_stasis,lnl_line_stasis,notu);
    last_alphas=dvector(rates);
    equaldvector(last_alphas,alphas,rates);
    last_betas=dvector(rates);
    equaldvector(last_betas,betas,rates);
    last_ancestral=ulvector((2*notu)-1);
    equalulvector(last_ancestral,ancestral,(2*notu)-1);
    last_minbr=ulvector((2*notu)-1);
    equalulvector(last_minbr,minbr,(2*notu)-1);
    last_nodes=nodes;

    
    for (r=0; r<runs; ++r)		{
	/* here, we either:
		1. alter the cladogram;
			adjust this to undo polytomies sometimes
		2. alter a branch length
			NOTE: this could create AD statements, or break them up.
				this could also create polytomies
		3. alter the anagenetic rate of change
		4. alter the pulsed rate of change
		5. collapse branches into polytomies
	************************************************************/
	y=(((unsigned int) rand())%RAND_MAX)/((double) RAND_MAX);
	if (r==0)	y=0.985;
	if (r==1)	y=0.0f;
	if (y<0.3)			SWITCH=SW_TREE;
	else if (y<0.9)		SWITCH=SW_BRANCH;
	else if (y<0.95)	SWITCH=SW_ALPHA;
	else if (polytomies=='y' && y<0.99)	SWITCH=SW_POLYT;
//	void swap_and_recalculate(long *ape, int notu, int nodes, long **statechars, long *charswstates, int maxst, double *alphas, double *betas, int rates, double ***statelikes, double lambda, double **ranges, double **divergences)
	if (SWITCH==SW_TREE || (SWITCH==SW_BRANCH || SWITCH==SW_POLYT))   {
		if (SWITCH==SW_TREE)
			swap_and_recalculate(ape, notu, nodes, statechars, charswstates, maxst, alphas, betas, rates, statelikes, lambda, ranges, divergences, ancestral, minbr);
		else if (SWITCH==SW_BRANCH)
			change_branch_length(ape, notu, nodes, statechars, charswstates, nchars, maxst, alphas, betas, rates, statelikes, lambda, ranges, divergences, ancestral, minbr);
		else if (SWITCH==SW_POLYT)
			collapse_node_and_recalculate(ape, notu, &nodes, statechars, charswstates, nchars, maxst, alphas, betas, rates, statelikes, lambda, ranges, divergences, ancestral, minbr);
		ma_htu=sum_divergence_times(divergences,notu,nodes);
		lnL_tree=get_tree_likelihood(statelikes[notu], novels, nstates, nchars);
		convertApeToNexus(ape,tree,notu,nodes);
        }
//		change_branch_length(ape, notu, nodes, statechars, charswstates, maxst, alphas, betas, rates, statelikes, lambda, ranges, divergences, ancestral, minbr);
	else if (SWITCH==SW_ALPHA || SWITCH==SW_BETA)   {
		if (SWITCH==SW_ALPHA)	{
			change_alpha(alphas,rates);
			stasis_likelihood_simple(notu, chmatrix, ranges, lnl_line_stasis, statechars, charswstates, novels, maxst, alphas, rates, UNKNOWN, INAP);
			lnL_stasis=sumdvector(lnl_line_stasis,notu);
			}
		else if (SWITCH==SW_BETA)
			change_alpha(betas,rates);
        /* recalculate P[stasis | alphas]   */
		MkSimpleWholeTree(ape, notu, nodes, divergences, statechars, charswstates, maxst, alphas, betas, rates, statelikes, lambda, ancestral, minbr);
		lnL_tree=get_tree_likelihood(statelikes[notu], novels, nstates, nchars);
		}
//	swap_and_recalculate(ape, notu, nodes, divergences, statechars, charswstates, maxst, alphas, rates, statelikes, lambda, BIF, INAP, UNKNOWN);
  	printf("Tree #%d: lnL=%4.3f; Sum divergences = %4.1f; ",++t,(lnL_tree+lnL_stasis), ma_htu);
	write_nexus(tree,0,notu);
	met_hast=exp(last_lnL-(lnL_tree+lnL_stasis));
	x=1.0f;
	if (met_hast<1)	x=(((unsigned int) rand())%RAND_MAX)/((double) RAND_MAX);
	if (x<=met_hast)	{
		/* these always will change */
        equaldcube(last_statelikes,statelikes,(2*notu)-1,nchars,maxst);
		last_lnL=lnL_tree+lnL_stasis;
		last_lnL_tree=lnL_tree;
		if (SWITCH==SW_TREE)    equallvector(last_ape,ape,(2*notu)-1);
		if (SWITCH==SW_TREE || (SWITCH==SW_BRANCH || SWITCH==SW_POLYT))	{
            equaldmatrix(last_divergences,divergences,(2*notu)-1,2);
			equalulvector(last_ancestral, ancestral, (2*notu)-1);
			equalulvector(last_minbr,minbr,(2*notu)-1);
			}
		if (SWITCH==SW_ALPHA)   {
            last_lnL_stasis=lnL_stasis;
            equaldvector(last_alphas,alphas,rates);
			}
		if (SWITCH==SW_BETA)	equaldvector(last_betas,betas,rates);
		if (SWITCH==SW_POLYT)	last_nodes=nodes;
		}	/* keepers	*/
	else	{
		equaldcube(statelikes,last_statelikes,(2*notu)-1,nchars,maxst);
		if (SWITCH==SW_TREE)	equallvector(ape,last_ape,(2*notu)-1);
		if (SWITCH==SW_TREE || (SWITCH==SW_BRANCH || SWITCH==SW_POLYT))	{
            equaldmatrix(divergences,last_divergences,(2*notu)-1,2);
			equalulvector(ancestral, last_ancestral, (2*notu)-1);
			equalulvector(minbr,last_minbr,(2*notu)-1);
			}
		if (SWITCH==SW_ALPHA)   {
            lnL_stasis=last_lnL_stasis;
            equaldvector(alphas,last_alphas,rates);
            }
		if (SWITCH==SW_BETA)	equaldvector(betas,last_betas,rates);
		if (SWITCH==SW_POLYT)	nodes=last_nodes;
		}	/* go back to prior parameters	*/
		
	//outfile1=fopen("/Users/peterjwagner/Documents/Projects/PaleoMCMC/Ape_Test_Autumn_II.xls","a");
        outfile1=fopen("Ape_Test_Autumn_II.xls","a");
	fprintf(outfile1,"%d",a);
	for (i=0; i<ttlbr; ++i)	
            fprintf(outfile1,"\t%d",apeT[i]);
	fprintf(outfile1,"\n");
	fclose(outfile1);
    }

  
    cout << "Done with processing...\n";
    

    
    return 0;
}

