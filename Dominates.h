//int Dominates(population *x, population *y);

int Dominates(population *x, population *y)
{
	float *x_fitn_ptr;
	float *y_fitn_ptr;
	int flag[nfunc];
	int ans = 0;
	for(int i = 0; i < nfunc; i++)
	{
		if(x->ind[i].fitness[i] < y->ind[i].fitness[i])
		{
			flag[i] = 1;
		}
		else if(x->ind[i].fitness[i] > y->ind[i].fitness[i])
		{
			flag[i] = 0;
		}
	}
	for(int j = 0; j < nfunc; j++)
	{
		if(flag[j] == 0)
		{
			ans = 0;
		}
		else if(flag[j] == 1)
		{
			ans = 1;
		}
	}
	
	return ans;
}