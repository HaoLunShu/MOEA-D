/*This also demarkates the different Pareto Fronts*/

void rankcon(population *pop_ptr);

int indcmp3(float *ptr1,float *ptr2);

void rankcon(population *pop_ptr)
{
  int i,j,k,       /*counters*/
    rnk,           /*rank*/
    val,           /*value obtained after comparing two individuals*/
    nondom,        /*no of non dominated members*/
    maxrank1,      /*Max rank of the population*/
    rankarr[maxpop], /*Array storing the individual number at a rank*/
    q;


  float *ptr1,*ptr2,*err_ptr1,*err_ptr2;

  /*------------------------------* RANKING *------------------------------*/

  /*Initializing the ranks to zero*/
  rnk = 0 ;

  nondom = 0 ;
  maxrank1 = 0;
  /*min_fit is initialize to start distributing the dummy fitness = 
    popsize to the rank one individuals and keeping the record such 
    that the minimum fitness of the better rank individual is always 
    greater than max fitness of the relatively worse rank*/
  
  min_fit = popsize;


  /*Difference in the fitness of minimum dummy fitness of better rank 
    and max fitness of the next ranked individuals*/
  
  delta_fit = 0.1 * popsize;

  /*Initializing all the flags to 2*/

  for( j = 0 ;j < popsize; j++)
    {
      pop_ptr->ind[j].flag = 2;
    }

  q = 0;

  for(k =  0;k < popsize;k++,q=0)
    {
      for(j = 0;j < popsize;j++)
	{
	  if (pop_ptr->ind[j].flag != 1)break;
	  /*Break if all the individuals are assigned a rank*/
	}
      if(j == popsize)break;

      rnk = rnk + 1;

      for( j = 0 ;j < popsize; j++)
	{
	  if(pop_ptr->ind[j].flag == 0) pop_ptr->ind[j].flag = 2;
	  /*Set the flag of dominated individuals to 2*/
	}

      for(i = 0;i < popsize ; i++)
	{
	  /*Select an individual which rank to be assigned*/


	  if(pop_ptr->ind[i].flag != 1 && pop_ptr->ind[i].flag != 0) 
	    {
	      ptr1 = &(pop_ptr->ind[i].fitness[0]);
	      err_ptr1 = &(pop_ptr->ind[i].error);

	      for(j = 0;j < popsize ; j++)
		{
		  
		  /*Select the other individual which has not got a rank*/
		  if( i!= j)
		    {
		      if(pop_ptr->ind[j].flag != 1)
			{
			  ptr2 = &(pop_ptr->ind[j].fitness[0]);
			  err_ptr2 = &(pop_ptr->ind[j].error);
			  
			  if(*err_ptr1 < 1.0e-6 && *err_ptr2 > 1.0e-6)
			    {
			      /*first ind is feasible second individaul 
				is infeasible*/
			      pop_ptr->ind[j].flag = 0;
			    }
			  else
			    {
			      if(*err_ptr1 > 1.0e-6 && *err_ptr2 < 1.0e-6)
				{
				  /*first individual is infeasible and 
				    second is feasible*/
				  pop_ptr->ind[i].flag = 0;
				  break;
				}
			      else
				{
				  /*both are feasible or both are infeasible*/
				  if(*err_ptr1 > *err_ptr2)
				    {
				      pop_ptr->ind[i].flag = 0;
				      /*first individual is more infeasible*/
				      break;
				    }
				  else
				    {
				      if(*err_ptr1 < *err_ptr2)
					{
					  pop_ptr->ind[j].flag = 0;
					  /*second individual is more 
					    infeasible*/
					}
				      else
					{
					  /*Compare the two individuals for 
					    fitness*/
					  val = indcmp3(ptr1,ptr2);
					  
					  /*VAL = 2 for dominated individual 
					    which rank to be given*/

					  /*VAL = 1 for dominating individual 
					    which rank to be given*/

					  /*VAL = 3 for non comparable 
					    individuals*/

					  if( val == 2)
					    { 
					      pop_ptr->ind[i].flag = 0;
					      /* individual 1 is dominated */
					      break;
					    }
					  
					  if(val == 1)
					    {
					      pop_ptr->ind[j].flag = 0;
					      /* individual 2 is dominated */
					    }
					  
					  if(val == 3)
					    {
					      nondom++;
					      /* individual 1 & 2 are 
						 non dominated */
					      if(pop_ptr->ind[j].flag != 0)
						pop_ptr->ind[j].flag = 3;
					    }
					  
					}   /*if loop ends*/
				    }       /* i != j loop ends*/
				}
			    }
			}
		    }
		}        /*loop over j ends*/
	      if( j == popsize)
		{
		  /*Assign the rank and set the flag*/
		  pop_ptr->ind[i].rank = rnk;
		  pop_ptr->ind[i].flag = 1;
		  rankarr[q] = i;
		  q++;
		}
	    }       /*Loop over flag check ends*/
	}           /*Loop over i ends */
      pop_ptr->rankno[rnk-1] = q ;
    }
  maxrank1 = rnk;
  
  
  /*     Find Max Rank of the population    */
  for(i = 0;i < popsize;i++)
    {
      rnk = pop_ptr->ind[i].rank;

      if(rnk > maxrank1)maxrank1 = rnk;

    }

  pop_ptr->maxrank = maxrank1;

  return;
}



/*Routine Comparing the two individuals*/

int indcmp3(float *ptr1,float *ptr2)
{
  float fit1[maxfun],fit2[maxfun];
  int i,value,m,n;
  for(i = 0;i < nfunc ;i++)
    {
      fit1[i] = *ptr1++;
      fit2[i] = *ptr2++;
    }
  m = 0;
  n = 0;
  while(m < nfunc && fit1[m] <= fit2[m]) 
    {
      if(fit1[m]== fit2[m]) n++;
      m++;
    }
  if(m == nfunc) 
    {
      if(n == nfunc) value = 3;
      else value = 1;             /*value = 1 for dominationg*/
    }
  else 
    {
      m = 0;
      n = 0;
      while(m < nfunc && fit1[m] >= fit2[m])
	{
	  if(fit1[m]== fit2[m]) n++;
	  m++;
	}
      if(m == nfunc)
	{
	  if(n != nfunc)
	    value = 2;                       /*value =  2 for dominated */
	  else value =3;
	}
      else value = 3;                   /*value = 3 for incomparable*/
    }
  
  return value;
}



