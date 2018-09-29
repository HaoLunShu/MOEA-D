/*This program subroutine is used to print the report*/

void report(int t,population *pop1_ptr, population *pop2_ptr,FILE *rep_ptr,FILE *gen_ptr, FILE *lastit);

void report(int t,population *pop1_ptr,population *pop2_ptr,FILE *rep_ptr,FILE *gen_ptr, FILE *lastit )
{
  int i,j,*be,*rptr,*rptr1; 

  float *ptr1,*ptr,*fptr,*fptr1, *ptr1_b, *ptr2_b;

  float *ptr2,*cons_ptr1,*cons_ptr2, *err2;
  
  int *gene_ptr1 ,*gene_ptr2 ;

  fprintf(rep_ptr,"\n\n---------------------------------------------------\n");
  fprintf(rep_ptr,"Generation No.     ->%d\n",t+1);
  fprintf(rep_ptr,"------------------------------------------------------\n");
  if(ncons == 0)
    fprintf(rep_ptr," variables (real %d binary %d)  fitness (%d)  rank cublen || variables  fitness rank cublen\n",nvar,nchrom,nfunc);
  else
    fprintf(rep_ptr," variables (real %d binary %d)  fitness (%d) constraint (%d) penalty rank cublen || variables  fitness constraint penalty rank cublen\n",nvar,nchrom,nfunc,ncons);
  


  for(i = 0;i < popsize;i++)
    {
      fprintf(rep_ptr,"\n------------------------------------------------\n");

      for(j = 0;j < nvar;j++)
	fprintf(rep_ptr,"%f ",pop1_ptr->ind[i].xreal[j]);

      for(j = 0;j < nchrom; j++)
	fprintf(rep_ptr,"%f ",pop1_ptr->ind[i].xbin[j]);
      if (t == gener-1)
	{
	  for(j = 0;j < nfunc;j++)
	    {   
	      if ((pop2_ptr->ind[i].error <= 0.0) && (pop2_ptr->ind[i].rank == 1))
		fprintf(lastit,"%f\t",pop2_ptr->ind[i].fitness[0]);
	    }
	  if ((pop2_ptr->ind[i].error <= 0.0) && (pop2_ptr->ind[i].rank == 1))
	    fprintf(lastit,"\n");
	}
      
      for(j = 0;j < nfunc;j++)
	fprintf(rep_ptr,"  %.4f",pop1_ptr->ind[i].fitness[j]);
      
      if(ncons != 0)
	{
	  for(j = 0;j < ncons;j++)
	    fprintf(rep_ptr,"  %.2e",pop1_ptr->ind[i].constr[j]);
	  fprintf(rep_ptr," %.2e",pop1_ptr->ind[i].error);
	}
      
      fprintf(rep_ptr," %d ",pop1_ptr->ind[i].rank);
      
      fprintf(rep_ptr,"|**|");

      for(j = 0;j < nvar;j++)
	{
	  fprintf(rep_ptr," %f ",pop2_ptr->ind[i].xreal[j]);
	  fprintf(gen_ptr," %f ",pop2_ptr->ind[i].xreal[j]);
	}
      for(j = 0;j < nchrom; j++)
	{
	  fprintf(rep_ptr,"%f ",pop2_ptr->ind[i].xbin[j]); 
	  fprintf(gen_ptr,"%f ",pop2_ptr->ind[i].xbin[j]);
	}
      for(j = 0;j < nfunc;j++)
	{	
	  fprintf(rep_ptr,"  %f",pop2_ptr->ind[i].fitness[j]);
	  fprintf(gen_ptr,"  %f",pop2_ptr->ind[i].fitness[j]);
	}
      fprintf(rep_ptr," %d ",pop2_ptr->ind[i].rank);
      
      if(ncons != 0)
	{
	  for(j = 0;j < ncons;j++)
	    {
	      fprintf(rep_ptr,"  %.2e",pop2_ptr->ind[i].constr[j]);
	      fprintf(gen_ptr,"  %.2e",pop2_ptr->ind[i].constr[j]);
	    }
	  fprintf(rep_ptr," %.2e",pop2_ptr->ind[i].error);
	  fprintf(gen_ptr," %.2e",pop2_ptr->ind[i].error);
	}
      
      fprintf(gen_ptr,"\n");
    }
  
  fprintf(rep_ptr,"\n--------------------------------------------------\n\n"); 
  fprintf(rep_ptr,"-------------------------------------------------------\n");
  fprintf(gen_ptr,"\n--------------------------------------------------\n\n");
  return;
}






