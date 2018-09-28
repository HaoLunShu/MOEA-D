/* This is the module used to formulate the mutation routine*/

void mutate(population *new_pop_ptr);             

void mutate(population *new_pop_ptr)
{
  int i,*ptr,j,r;
  float rand1,*rand_float_ptr;

  rand1=randomperc();
  
  for(j = 0;j < popsize;j++)
    {
      ptr= &(new_pop_ptr->ind[j].genes[0]);     
      
      /*Select bit */
      for (i = 0;i < chrom;i++)
	{
	  rand1 = randomperc();
	  
	  /*Check whether to do mutation or not*/
	  if(rand1 <= pmut_b)
	    {
	      if(*ptr == 0)
		*ptr =1;
	      else
		*ptr=0;
	      nmut++;
	    }
	  ptr++;
	}
    }
  return;
}
