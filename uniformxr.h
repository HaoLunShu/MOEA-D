/* This is the header file to do the uniform crossover */

void unicross(population *new_pop_ptr,population *mate_pop_ptr);

void unicross(population *new_pop_ptr, population *mate_pop_ptr)
{
  int i,j,r,*gene,y,n,*par1,*par2,*chld1,*chld2;
  float rnd;
  for(i = 0,y = 0,n = 0;i < popsize;i++)
    {
      for(j = 0;j < chrom;j++)
	{

	  /*Select a bit for doing cross-over*/	
	  chld1 = &(new_pop_ptr->ind[y].genes[j]);

	  chld2 = &(new_pop_ptr->ind[y+1].genes[j]);
	
	  par1 = &(mate_pop_ptr->ind[n].genes[j]);

	  par2 = &(mate_pop_ptr->ind[n+1].genes[j]);

	  rnd = randomperc();

	  /*Checking whether to do cross-over or not*/
	  if(rnd <= pcross)
		{
		  ncross++;
		  *chld1 = *par2;
		  *chld2 = *par2;
		}

	  else
	    {
	      *chld1 = *par1;
	      *chld2 = *par2;
	    }
	}

      y = y+2;
      n = n+1;

    }

  for(i = 0;i < popsize;i++)
    {
      gene = &(new_pop_ptr->ind[i].genes[0]);
      for(j = 0;j < chrom;j++)
	{
	  gene = &(new_pop_ptr->ind[i].genes[j]);
	}
    }
  return;
}


