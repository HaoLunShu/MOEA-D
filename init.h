/*This is the file which initializes the population*/
void init(population *pop_ptr);

void init(population *pop_ptr)
{
  int i,j,r;
  float d; 
  
  /*Loop Over the population size*/
  for (i = 0 ; i < popsize ; i++)
    { 
      
      /*Loop over the chromosome length*/
      for (j = 0;j < chrom;j++)
	{
	  /*Generate a Random No. if it is less than 0.5 it 
	    generates a 0 in the string otherwise 1*/
	  d = randomperc();
	  if(d >= 0.5)
	    {
	      pop_ptr->ind[i].genes[j] = 1;
	    }
	  else
	    {
	      pop_ptr->ind[i].genes[j] = 0;
	    } 
	}
    }
 return;
}






























































