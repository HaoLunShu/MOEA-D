//#include "Dominates.h"

void DetermineDomination(population *pop_ptr);

void DetermineDomination(population *pop_ptr)
{
	population *tmp_pop_ptr;
	tmp_pop_ptr = pop_ptr;
	for(int i = 0; i < popsize; i++)
	{
		pop_ptr->IsDominated = 0;
		//pop_ptr = pop_ptr->next;
	}
	pop_ptr = tmp_pop_ptr;
	for(int j = 0; j < popsize; j++)
	{
		//tmp_pop_ptr = tmp_pop_ptr->next;
		for(int k = 1; k < popsize; k++)
		{
			if(Dominates(pop_ptr, tmp_pop_ptr) == 1)
			{
				tmp_pop_ptr->IsDominated = 1;
			}
			else if(Dominates(tmp_pop_ptr, pop_ptr) == 0)
			{
				pop_ptr->IsDominated = 0;
			}
		}
		//pop_ptr = pop_ptr->next;
	}
}