/* This is a Multi-Objective GA program.
**********************************************************************
*  This program is the implementation of the NSGA-2 proposed by      *
*                                                                    *
*  Prof. Kalyanmoy Deb and his students .                            *
*                                                                    *
*  copyright Kalyanmoy Deb
**********************************************************************

18.08.2003: The keepaliven.h file is modified to have normalized
            crowding distance calculation. The previous version of 
            the code did not have this feature. This way, maintaining
            a good distribution of solutions in problems having quite
            a different range of objective functions were difficult.
            Hopefully, with this modification, such difficulties will
            not appear. --  K. Deb
18.08.2003: Also the dfit.h file is deleted. It was not needed any way.

The user have to give the input manualy or through a data file.

The user needs to enter objective functions in func-con.h
The code can also take care of the constraints. Enter the constraints
in the space provided in the func-con.h file.
Constraints must be of the following type:
g(x) >= 0.0
Also normalize all constraints (see the example problem in func-con.h)

If your program asks you to increase the values of some parameters in the
program come to main program and accordingly changed the values which are
defined against #define ...

The program generates few output files. These are described as
1.output.out
*           This file has the detailed record for all the variables,
*           the fitness values, constraint values, overall constraint
            violation (penalty)  and their ranks for all the members
*           of old population in the left hand side of the |**|
*           and of new population in the right hand side.

2.all_fitness.out
*         This file prints the record of all the fitness values for
*         different individual of new popultion created at all
*         generations.

3.g_rank_record.out
*        This file maintains the record of individuals in global pop-
*        -ulation at different ranks for all the generations.

4.ranks.out
*         This file prints the number of individual at different ranks
*          in old and new population and finds rank ratios

5.final_fitness.out
*                 This file has the fitness value of all feasible and
                  non-dominated individuals at the final generation

6.final_var.out
*                 This file has the all the variables of the feasible
                  and non-dominated individuals at the final generation.
                  The i-th solutions here corresponds to the i-th solution
                  in the final_fitness.out file. 

7.plot.out        This file contains gnuplot-based file for plotting
                  the non-dominated feasible solutions obtained by the code.
*************************************************************************
*         This is recommended to delete or rename all the *.out files
*         obtained from the previous runs as some files are opened in
*         append mode so they give false resemblence of data if the
*         user is not careful

Compilation procedure:  gcc nsga2.c -lm
Run ./a.out with or without an input file

Input data files: Three files are included, but at one time one is needed
depending on the type of variables used:
inp-r (template file input-real)  : All variables are real-coded
inp-b (template file input-binary): All variables are binary-coded
inp-rb(template file input-rl+bin): Some variables are real and some are binary  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define square(x) ((x)*(x))
#define maxpop   500  /*Max population */
#define maxchrom 1000  /*Max chromosome length*/
#define maxvar    30  /*Max no. of variables*/
#define maxfun    10  /*Max no. of functions */
#define maxcons   20  /*Max no. of Constraints*/

int gener,       /*No of generations*/
  nvar,nchrom,          /*No of variables(real) and chromosome(binary-codes)*/
  ncons,         /*No of Constraints*/
  vlen[maxvar],  /*Array to store no of bits for each variable*/
  nmut,          /* No of Mutations */
  ncross,        /*No of crossovers*/
  ans;
float seed,      /*Random Seed*/
  pcross,        /*Cross-over Probability*/
  pmut_b, pmut_r,          /*Mutation Probability*/
  lim_b[maxvar][2], lim_r[maxvar][2];/*Limits of variable in array*/
float di,        /*Distribution Index for the Cross-over*/
  dim,           /*Distribution Index for the Mutation*/
  delta_fit,     /* variables required forfitness for fitness sharing */
  min_fit,
  front_ratio;
int optype,      /*Cross-over type*/
  nfunc,         /*No of functions*/
  sharespace;    /*Sharing space (either parameter or fitness)*/
  
float z[maxfun];	/*Goal Point*/

double coef[maxvar]; /*Variable used for decoding*/

static int popsize,  /*Population Size*/
  chrom;             /*Chromosome size*/

typedef struct       /*individual properties*/
{
  int genes[maxchrom], /*bianry chromosome*/
    rank,              /*Rank of the individual*/
    flag;              /*Flag for ranking*/
  float xreal[maxvar], /*list of real variables*/
    xbin[maxvar];      /*list of decoded value of the chromosome */
  float fitness[maxfun],/*Fitness values(cost) */
    constr[maxcons],     /*Constraints values*/
    //cub_len,             /*crowding distance of the individual*/
    error;              /* overall constraint violation for the individual*/
  float Position[maxvar],
    g;
}individual;        /*Structure defining individual*/

int T = 0;

typedef struct
{
	double lambda[maxfun];
	int Neighbors[maxpop],
		*NB_ptr;
}subproblems;

typedef struct
{
  int maxrank;            /*Maximum rank present in the population*/
  int IsDominated;
  float rankrat[maxpop];  /*Rank Ratio*/
  int rankno[maxpop];     /*Individual at different ranks*/
  individual ind[maxpop], /*Different Individuals*/
    *ind_ptr;
  subproblems sp[maxpop],
	*sp_ptr;
}population;             /*Popuation Structure*/


#include "random.h"       /*Random Number Generator*/

#include "input.h"        /*File Takes Input from user*/

//#include "FileInput.h"    /*File Takes Input from file*/

#include "realinit.h"     /*Random Initialization of the populaiton*/

#include "init.h"         /*Random Initialization of the population*/

//#include "CreateSubProblems.h"

#include "decode.h"       /*File decoding the binary dtrings*/

#include "ranking.h"      /*File Creating the Pareto Fronts*/

#include "rancon.h"       /*File Creating the Pareto Fronts when
			    Constraints are specified*/

#include "func-con.h"     /*File Having the Function*/

#include "DecomposedCost.h"

#include "DetermineDomination.h"

#include "Dominates.h"

#include "select.h"       /*File for Tournament Selection*/

#include "crossover.h"    /*Binary Cross-over*/

#include "uniformxr.h"    /*Uniform Cross-over*/

#include "realcross2.h"   /*Real Cross-over*/

#include "mut.h"          /*Binary Mutation*/

#include "realmut1.h"     /*Real Mutation*/

#include "keepaliven.h"   /*File For Elitism and Sharing Scheme*/

#include "report.h"       /*Printing the report*/

population oldpop,
  newpop,
  matepop,
  EP,
  *old_pop_ptr,
  *new_pop_ptr,
  *mate_pop_ptr,
  *EP_ptr;
/*Defining the population Structures*/

/*subproblems sp,
  *c_sp_ptr;*/

double norm(double x[])
{
	double m_sum = 0.0;
	for (int i=0; i < nfunc; i++)
	{
		m_sum  += pow(x[i], 2);
	}
	return sqrt(m_sum);
}

double CalculateDistance(double x1, double y1, double x2, double y2)
{
	double diffx = x1 - x2;
	double diffy = y1 - y2;
	double diffx_sqr = pow(diffx,2);
	double diffy_sqr = pow(diffy,2);
	double distance = sqrt(diffx_sqr + diffy_sqr);
	
	return distance;
}

void CreateSubProblems()
{
	double SO[maxpop][maxpop];
	/*subproblems *sp_tmp_ptr;
	sp_tmp_ptr = (subproblems*) malloc(sizeof(subproblems));
	double lamb[nfunc];
	double A[popsize], B[popsize], C[popsize];
	double D[maxpop][maxpop];
	subproblems *tmp_ptr;
	tmp_ptr = sp_tmp_ptr;
	subproblems *tmp1_ptr;
	tmp1_ptr = sp_tmp_ptr;
	for(int i = 0; i < popsize; i++)
	{
		for(int j = 0; j < nfunc; j++)
		{
			srand( time(NULL) );
			/* 產生 [0, 1) 的浮點數亂數 
			lamb[j] = (double) rand() / (RAND_MAX + 1.0);
		}
		double n = norm(lamb);
		for(int k = 0; k < nfunc; k++)
		{
			lamb[k] = lamb[k] / n;
		}
		for(int l = 0; l < nfunc; l++)
		{
			sp_tmp_ptr->lambda[l] = lamb[l];
		}
		//sp_ptr = sp_ptr->next;
	}
	//sp_ptr->next = 0;
	sp_tmp_ptr = tmp_ptr;
	for(int i = 0; i < popsize; i++)
	{
		for(int a = 0; a < nfunc; a++)
		{
			A[a] = sp_tmp_ptr->lambda[a];
		}
		for(int j = 0; j < popsize; j++)
		{
			for(int b = 0; b < nfunc; b++)
			{
				B[b] = sp_tmp_ptr->lambda[b];
			}
			D[i][j] = CalculateDistance(A[0], A[1], B[0], B[1]);
			//tmp_ptr = tmp_ptr->next;
		}
		//sp_ptr = sp_ptr->next;
	}
	sp_tmp_ptr = tmp1_ptr;
	for(int i = 0; i < popsize; i++)
	{
		double tmp[popsize];
		for(int c = 0; c < popsize; c++)
		{
			C[c] = D[i][c];
			tmp[c] = C[c];
		}
		for(int x = 1; x < popsize; x++)
		{
			if (tmp[x] > tmp[i])                //Comparing other array elements
			{
				int tmp_1 = tmp[i];         //Using temporary variable for storing last value
				tmp[i] = tmp[x];            //replacing value
				tmp[x] = tmp_1;             //storing last value
			}
		}
		for(int k = 0; k < popsize; k++)
		{
			SO[i][k] = tmp[k];
		}
		for(int j = 0; j < T; j++)
		{
			sp_tmp_ptr->Neighbors[j] = SO[i][j];
		}
		//sp_ptr = sp_ptr->next;
	}*/
	return;
}
  
main(int argc, char *argv[])
{

	/*Some Local variables to this Problem (Counters And some other pointers*/

	int i,j,l,f,maxrank1;
	float *ptr,tot,*fitn_ptr;
	FILE
		*fp,  
		*rep_ptr,
		*gen_ptr,
		*rep2_ptr,
		*end_ptr,
		*g_var,
		*lastit;
	/*File Pointers*/
	
	/*subproblems *sp1_ptr;
	sp1_ptr = (subproblems*) malloc(sizeof(subproblems));*/
	
	CreateSubProblems();

	rep_ptr = fopen("output.out","w");
	gen_ptr =fopen("all_fitness.out","w");
	rep2_ptr = fopen("ranks.out","w");
	end_ptr = fopen("final_fitness.out","w");
	g_var = fopen("final_var.out","w");
	lastit = fopen("plot.out","w");
	/*Opening the files*/
	
	old_pop_ptr = (population*) malloc(sizeof(population));

	old_pop_ptr = &(oldpop);

	nmut = 0;
	ncross = 0;

	//if(argc == 1)
	//{
		/*Get the input from the file input.h*/
		input(rep_ptr);
	/*}
	else if((fp = fopen(argv[1], "r")) == NULL)
	{
		printf("open_file_error");
		exit(1);
	}
	else
	{
		FileInput(fp, rep_ptr);
	}*/

	fprintf(rep_ptr,"Results in a file\n");
	fprintf(end_ptr,"# Last generation population (Feasible and non-dominated)\n");
	fprintf(end_ptr,"# Fitness_vector (first %d)  Constraint_violation (next %d)  Overall_penalty\n",nfunc,ncons);
	fprintf(g_var,"#Feasible Variable_vectors for non-dominated solutions at last generation\n");
	fprintf(g_var,"# Real (first %d)  Binary (next %d)\n",nvar,nchrom);
	fprintf(lastit,"# Feasible and Non-dominated Objective Vector\n");

	if(ceil(0.15*popsize) > 2)
		T = ceil(0.15*popsize);
	else
		T = 2;
	if(T < 2)
		T = 2;
	if(T > 15)
		T = 15;
	//T = min(max(T, 2), 15);	/*Number of Neighbors*/
	/*Initialize the random no generator*/
	warmup_random(seed);

	/*Binary Initializaton*/
	if (nchrom > 0)
	{
		init(old_pop_ptr);  
	}
	if (nvar > 0)
	{
		realinit(old_pop_ptr);
	}
	
	EP_ptr = (population*) malloc(sizeof(population));
	EP_ptr->ind_ptr = &(EP_ptr->ind[0]);
	
	for(int x = 0; x < popsize; x++)
	{
		for(int y = 0; y < chrom; y++)
		{
			EP_ptr->ind_ptr->genes[y] = 0;
		}
		EP_ptr->ind_ptr = &(EP_ptr->ind[x+1]);
		for(int z = 0; z < nvar; z++)
		{
			EP_ptr->ind[x].xreal[z] = 0;
		}
	}

	old_pop_ptr = &(oldpop);
	old_pop_ptr->sp_ptr = &(old_pop_ptr->sp[0]);
	//c_sp_ptr = &(sp);
	
	// Create Sub-problems
	//CreateSubProblems(old_pop_ptr->sp_ptr);

	// decode binary strings
	decode(old_pop_ptr); 

	old_pop_ptr = &(oldpop);
	new_pop_ptr = &(newpop);

	for(j = 0;j < popsize;j++)
	{
		/*Initializing the Rank array having different individuals
		at a particular  rank to zero*/
		old_pop_ptr->rankno[j] = 0;
		new_pop_ptr->rankno[j] = 0;
	}
	
	for(int k = 0; k < nfunc; k++)
	{
		z[k] = 100;
	}

	old_pop_ptr = &(oldpop);
	

	func(old_pop_ptr); 
	/*Function Calculaiton*/
	
	old_pop_ptr = &(oldpop);
	old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
	
	for(int a = 0; a < popsize; a++)
	{
		old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[a]);
		fitn_ptr = &(old_pop_ptr->ind_ptr->fitness[0]);
		for(int b = 0; b < nfunc; b++)
		{
			fitn_ptr = &(old_pop_ptr->ind_ptr->fitness[b]);
			if(z[b] > *fitn_ptr)
			{
				z[b] = *fitn_ptr;
			}
		}
		//old_pop_ptr = old_pop_ptr->next;
	}
	
	old_pop_ptr = &(oldpop);
	old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
	old_pop_ptr->sp_ptr = &(old_pop_ptr->sp[0]);
	
	for(int b = 0; b < popsize; b++)
	{
		old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[b]);
		old_pop_ptr->sp_ptr = &(old_pop_ptr->sp[b]);
		old_pop_ptr->ind_ptr->g = DecomposedCost(old_pop_ptr->ind_ptr, z, old_pop_ptr->sp_ptr->lambda);
		//c_sp_ptr = c_sp_ptr->next;
	}
	
	DetermineDomination(old_pop_ptr);
	
	//EP_ptr = &(EP);
	old_pop_ptr = &(oldpop);
	
	for(int c = 0; c < popsize; c++)
	{
		if(old_pop_ptr->IsDominated == 1)
		{
			EP_ptr->maxrank = old_pop_ptr->maxrank;
			EP_ptr->IsDominated = old_pop_ptr->IsDominated;
			EP_ptr->rankrat[c] = old_pop_ptr->rankrat[c];
			EP_ptr->rankno[c] = old_pop_ptr->rankno[c];
			EP_ptr->ind[c] = old_pop_ptr->ind[c];
			EP_ptr->sp[c] = old_pop_ptr->sp[c];
		}
		//old_pop_ptr = old_pop_ptr->next;
	}
	
	/*Node* start;
	Node* ptr = (Node*) malloc(sizeof(Node));
	start = ptr;
	Node* newNode = (Node*) malloc(sizeof(Node));
	ptr=>next = newNode;
	ptr->next = (Node*)malloc(sizeof(Node));
	ptr->next->value = 10;*/

	fprintf(rep_ptr,"----------------------------------------------------\n");
	fprintf(rep_ptr,"Statistics at Generation 0 ->\n");
	fprintf(rep_ptr,"--------------------------------------------------\n");
	
	//population p1, p2;
	population *p1_ptr, *p2_ptr;
	
	//individual tmp_ind;

	/********************************************************************/
	/*----------------------GENERATION STARTS HERE----------------------*/
	for (i = 0;i < gener;i++)
	{
		printf("Generation = %d\n",i+1);
		old_pop_ptr = &(oldpop);
		//p1_tmp_ptr = &(p1);
		//p2_tmp_ptr == &(p2);
		mate_pop_ptr = &(matepop);
		old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
		old_pop_ptr->sp_ptr = &(old_pop_ptr->sp[0]);
		fprintf(rep_ptr,"Population at generation no. -->%d\n",i+1);
		fprintf(gen_ptr,"#Generation No. -->%d\n",i+1);
		fprintf(gen_ptr,"#Variable_vector  Fitness_vector Constraint_violation Overall_penalty\n");

		/*for(int j = 0; j < popsize; j++)
		{		
			int rnd1 = rnd(0,T);
			int rnd2 = rnd(0,T);

			int j1 = 0;
			int j2 = 0;
		
			p1_ptr = old_pop_ptr->ind_ptr;
			p2_ptr = old_pop_ptr->ind_ptr;
			
			j1 = old_pop_ptr->sp_ptr->NB_ptr = &(old_pop_ptr->sp_ptr->Neighbors[rnd1]);
			/*for(int x = 0; x <= j1; x++)
			{
				p1_ptr = p1_ptr->next;
			}
			p1_ptr = &(p1_ptr->ind[j1]);
			//p1_tmp_ptr = p1_ptr;
			//p1_tmp_ptr = p1_tmp_ptr->next;
			
			//mate_pop_ptr = p1_ptr;
			//mate_pop_ptr = mate_pop_ptr->next;
			
			j2 = old_pop_ptr->sp_ptr->NB_ptr = &(old_pop_ptr->sp_ptr->Neighbors[rnd2]);
			/*for(int x = 0; x <= j2; x++)
			{
				p2_ptr = p2_ptr->next;
			}
			p2_ptr = &(p2_ptr->ind[j2]);
			//p2_tmp_ptr = p2_ptr;
			//p2_tmp_ptr = p2_tmp_ptr->next;
			
			//mate_pop_ptr = p2_ptr;
			//mate_pop_ptr = mate_pop_ptr->next;
			
			//c_sp_ptr = c_sp_ptr->next;
		}*/

		/*--------SELECT----------------*/
		nselect(old_pop_ptr ,mate_pop_ptr );

		new_pop_ptr = &(newpop);
		mate_pop_ptr = &(matepop);

		/*CROSSOVER----------------------------*/      
		if (nchrom > 0) 
		{

			if(optype == 1)
			{
			  crossover(new_pop_ptr ,mate_pop_ptr );
			  /*Binary Cross-over*/
			}

			if(optype == 2)
			{
			  unicross(new_pop_ptr ,mate_pop_ptr );
			  /*Binary Uniform Cross-over*/
			}
		}
		if (nvar > 0)
		{			
			realcross(new_pop_ptr ,mate_pop_ptr );
			/*Real Cross-over*/
		}


		/*------MUTATION-------------------*/
		/*new_pop_ptr = &(newpop);

		if (nchrom > 0)
		{
			mutate(new_pop_ptr );
			Binary Mutation 
		}

		if (nvar > 0)
		{
			real_mutate(new_pop_ptr );*/
			/*Real Mutation
		}*/

		new_pop_ptr = &(newpop);

		/*-------DECODING----------*/
		if(nchrom > 0)
		{
			decode(new_pop_ptr );
			/*Decoding for binary strings*/
		}

		/*----------FUNCTION EVALUATION-----------*/
		new_pop_ptr = &(newpop);
		func(new_pop_ptr );
	
		new_pop_ptr = &(newpop);
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);
		
		for(int a = 0; a < popsize; a++)
		{
			new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[a]);
			fitn_ptr = &(new_pop_ptr->ind_ptr->fitness[0]);
			for(int b = 0; b < nfunc; b++)
			{
				fitn_ptr = &(new_pop_ptr->ind_ptr->fitness[b]);
				if(z[b] > *fitn_ptr)
				{
					z[b] = *fitn_ptr;
				}
			}
			//new_pop_ptr = new_pop_ptr->next;
		}

		/*-------------------SELECTION KEEPING FRONTS ALIVE--------------*/
		old_pop_ptr = &(oldpop);
		old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[0]);
		old_pop_ptr->sp_ptr = &(old_pop_ptr->sp[0]);
		new_pop_ptr = &(newpop);
		new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);
		new_pop_ptr->sp_ptr = &(new_pop_ptr->sp[0]);
		mate_pop_ptr = &(matepop);
		int *old_ptr, *new_ptr;
		float *old_ptr_r, *new_ptr_r;
		//c_sp_ptr->NB_ptr = &(c_sp_ptr->Neighbors[0]);
		
		for(int a = 0; a < popsize; a++)
		{
			for(int b = 0; b < T; b++)
			{
				old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[b]);
				old_ptr = &(old_pop_ptr->ind_ptr->genes[0]);
				old_ptr_r = &(old_pop_ptr->ind_ptr->xreal[0]);
				old_pop_ptr->ind_ptr->g = DecomposedCost(old_pop_ptr->ind_ptr, z, old_pop_ptr->sp_ptr->lambda);
				new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[b]);
				new_ptr = &(new_pop_ptr->ind_ptr->genes[0]);
				new_ptr_r = &(new_pop_ptr->ind_ptr->xreal[0]);
				new_pop_ptr->ind_ptr->g = DecomposedCost(new_pop_ptr->ind_ptr, z, new_pop_ptr->sp_ptr->lambda);
				//c_sp_ptr->NB_ptr = &(c_sp_ptr->Neighbors[b]);
				if(new_pop_ptr->ind_ptr->g <= old_pop_ptr->ind_ptr->g)
				{
					for(i = 0;i < chrom;i++)
						*old_ptr++=*new_ptr++;
					for(i = 0;i < nvar;i++)
						*old_ptr_r++=*new_ptr_r++;
				}
				/*c_sp_ptr = c_sp_ptr->next;
				old_pop_ptr = old_pop_ptr->next;
				new_pop_ptr = new_pop_ptr->next;*/
			}
		}
		
		old_pop_ptr = &(oldpop);
		
		DetermineDomination(old_pop_ptr);
		
		//EP_ptr = &(EP);
		old_pop_ptr = &(oldpop);
		
		for(int c = 0; c < popsize; c++)
		{
			if(old_pop_ptr->IsDominated == 1)
			{
				EP_ptr->maxrank = old_pop_ptr->maxrank;
				EP_ptr->IsDominated = old_pop_ptr->IsDominated;
				EP_ptr->rankrat[c] = old_pop_ptr->rankrat[c];
				EP_ptr->rankno[c] = old_pop_ptr->rankno[c];
				EP_ptr->ind[c] = old_pop_ptr->ind[c];
				EP_ptr->sp[c] = old_pop_ptr->sp[c];
			}
			//old_pop_ptr = old_pop_ptr->next;
		}

		/*Elitism And Sharing Implemented*/
		//keepalive(old_pop_ptr ,new_pop_ptr ,mate_pop_ptr,i+1);      

		//mate_pop_ptr = &(matepop);
		EP_ptr = &(EP);
		if(nchrom > 0)
		{
			//decode(mate_pop_ptr);
			decode(EP_ptr);
		}

		//mate_pop_ptr = &(matepop);
		EP_ptr = &(EP);
		/*------------------REPORT PRINTING--------------------------------*/  
		//report(i ,old_pop_ptr ,mate_pop_ptr ,rep_ptr ,gen_ptr, lastit );
		report(i ,old_pop_ptr ,EP_ptr ,rep_ptr ,gen_ptr, lastit);

		/*==================================================================*/

		/*----------------Rank Ratio Calculation------------------------*/
		new_pop_ptr = &(matepop);
		old_pop_ptr = &(oldpop);

		/*Finding the greater maxrank among the two populations*/

		if(old_pop_ptr->maxrank > new_pop_ptr->maxrank)
		{
			maxrank1 = old_pop_ptr->maxrank;
		}
		else 
		{
			maxrank1 = new_pop_ptr->maxrank;
		}

		fprintf(rep2_ptr,"--------RANK AT GENERATION %d--------------\n",i+1);
		fprintf(rep2_ptr,"Rank old ranks   new ranks     rankratio\n");

		for(j = 0;j < maxrank1 ; j++)
		{ 
			/*Sum of the no of individuals at any rank in old population 
			and the new populaion*/

			tot = (old_pop_ptr->rankno[j])+ (new_pop_ptr->rankno[j]);

			/*Finding the rank ratio for new population at this rank*/

			new_pop_ptr->rankrat[j] = (new_pop_ptr->rankno[j])/tot;

			/*Printing this rank ratio to a file called ranks.dat*/

			fprintf(rep2_ptr," %d\t  %d\t\t %d\t %f\n",j+1,old_pop_ptr->rankno[j],new_pop_ptr->rankno[j],new_pop_ptr->rankrat[j]);

		}

		fprintf(rep2_ptr,"-----------------Rank Ratio-------------------\n");
		/*==================================================================*/

		/*=======Copying the new population to old population======*/

		old_pop_ptr = &(oldpop);
		new_pop_ptr = &(matepop);

		for(j = 0;j < popsize;j++)
		{
			old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
			new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);
			if(nchrom > 0)
			{
				/*For Binary GA copying of the chromosome*/

				for(l = 0;l < chrom;l++)
				{
					old_pop_ptr->ind_ptr->genes[l]=new_pop_ptr->ind_ptr->genes[l];
				}

				for(l = 0;l < nchrom;l++)
				{
					old_pop_ptr->ind_ptr->xbin[l] = new_pop_ptr->ind_ptr->xbin[l];
				}
			}
			if(nvar > 0)
			{
				/*For Real Coded GA copying of the chromosomes*/
				for(l = 0;l < nvar;l++)
				{
					old_pop_ptr->ind_ptr->xreal[l] = new_pop_ptr->ind_ptr->xreal[l];
				}
			}

			/*Copying the fitness vector */	  
			for(l = 0 ; l < nfunc ;l++)
			{
				old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];
			}

			/*Copying the dummy fitness*/
			//old_pop_ptr->ind_ptr->cub_len = new_pop_ptr->ind_ptr->cub_len;

			/*Copying the rank of the individuals*/
			old_pop_ptr->ind_ptr->rank = new_pop_ptr->ind_ptr->rank;

			/*Copying the error and constraints of the individual*/

			old_pop_ptr->ind_ptr->error = new_pop_ptr->ind_ptr->error;
			for(l = 0;l < ncons;l++)
			{
				old_pop_ptr->ind_ptr->constr[l] = new_pop_ptr->ind_ptr->constr[l];
			}

			/*Copying the flag of the individuals*/
			old_pop_ptr->ind_ptr->flag = new_pop_ptr->ind_ptr->flag;
		}   // end of j

		maxrank1 = new_pop_ptr->maxrank ;

		/*Copying the array having the record of the individual 
		at different ranks */
		for(l = 0;l < popsize;l++)
		{
			old_pop_ptr->rankno[l] = new_pop_ptr->rankno[l];
		}

		/*Copying the maxrank */
		old_pop_ptr->maxrank = new_pop_ptr->maxrank;

		/*Printing the fitness record for last generation in a file last*/
		if(i == gener-1)
		{
			// for the last generation 
			old_pop_ptr = &(matepop);
			for(f = 0;f < popsize ; f++) // for printing
			{
				old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[f]);

				if ((old_pop_ptr->ind_ptr->error <= 0.0) && (old_pop_ptr->ind_ptr->rank == 1))  // for all feasible solutions and non-dominated solutions
				{
					for(l = 0;l < nfunc;l++)
					{
						fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->fitness[l]);
					}
					for(l = 0;l < ncons;l++)
					{
						fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->constr[l]);
					}
					if (ncons > 0)
					{
						fprintf(end_ptr,"%f\t",old_pop_ptr->ind_ptr->error);
						fprintf(end_ptr,"\n");
					}

					if (nvar > 0)
					{
						for(l = 0;l < nvar ;l++)
						{
							fprintf(g_var,"%f\t",old_pop_ptr->ind_ptr->xreal[l]);
						}
						fprintf(g_var,"  ");
					}

					if(nchrom > 0)
					{
						for(l = 0;l < nchrom;l++)
						{
							fprintf(g_var,"%f\t",old_pop_ptr->ind_ptr->xbin[l]);
						}
					}
					fprintf(g_var,"\n");
				}  // feasibility check
			} // end of f (printing)

		} // for the last generation
	}  // end of i 

	/*                   Generation Loop Ends                                */
	/************************************************************************/

	fprintf(rep_ptr,"NO. OF CROSSOVER = %d\n",ncross);
	fprintf(rep_ptr,"NO. OF MUTATION = %d\n",nmut);
	fprintf(rep_ptr,"------------------------------------------------------------\n");
	fprintf(rep_ptr,"---------------------------------Thanks---------------------\n");
	fprintf(rep_ptr,"-------------------------------------------------------------\n");
	printf("NOW YOU CAN LOOK IN THE FILE OUTPUT2.DAT\n");

	/*Closing the files*/
	fclose(rep_ptr);
	fclose(gen_ptr);
	fclose(rep2_ptr);
	fclose(end_ptr);
	fclose(g_var);
	fclose(lastit);
	
}





