/*This is the file for formulating the crossover process*/

void crossover(population *new_pop_ptr,population *mate_pop_ptr);

void crossover(population *new_pop_ptr,population *mate_pop_ptr)
{
  int i,j,k,l,m,n,y,mating_site,*par1,*par2,*chld1,*chld2,c;
  float rnd;
  int r;
  rnd=randomperc();  

  for (i = 0,y = 0,n = 0;i < popsize/2;i++)
    {
      chld1=&(new_pop_ptr->ind[n].genes[0]);
      n = n+1;
      
      chld2=&(new_pop_ptr->ind[n].genes[0]);
      n = n+1;
      
      par1 = &(mate_pop_ptr->ind[y].genes[0]);
      y = y+1;
      
      par2 = &(mate_pop_ptr->ind[y].genes[0]);  
      y = y+1;
      
      rnd = randomperc();
      if (rnd < pcross)
	{
	  ncross++;
	  rnd = randomperc();
	  c = floor(rnd*(chrom+10));
	  mating_site = c;
	  if(mating_site >= chrom)
	    {
	      mating_site = mating_site/2;
	    }
	  
	  for(k = 0;k < chrom;k++)
	    {
	      if(k > mating_site-1)
		{
		  *chld1++ = *par2++;
		  *chld2++ = *par1++;
		}
	      else
		{
		  *chld1++ = *par1++;
		  *chld2++ = *par2++;
		}
	    }
	}
    else 
      {
	for (k = 0;k < chrom;k++)
	  {
	    *chld1++ = *par1++;
	    *chld2++ = *par2++;
	  }
      }
    }
  return;
}

