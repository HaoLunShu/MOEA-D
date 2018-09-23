//int Dominates(population *x, population *y);

int Dominates(population *x, population *y)
{
	float *x_fitn_ptr;
	float *y_fitn_ptr;
	int flag[nfunc];
	int ans = 0;
	x->ind_ptr = &(x->ind[0]);
	y->ind_ptr = &(y->ind[0]);
	x_fitn_ptr = &(x->ind_ptr->fitness[0]);
	y_fitn_ptr = &(y->ind_ptr->fitness[0]);
	for(int i = 0; i < nfunc; i++)
	{
		x_fitn_ptr = &(x->ind_ptr->fitness[i]);
		y_fitn_ptr = &(y->ind_ptr->fitness[i]);
		if(*x_fitn_ptr < *y_fitn_ptr)
		{
			flag[i] = 1;
		}
		else if(*x_fitn_ptr > *y_fitn_ptr)
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