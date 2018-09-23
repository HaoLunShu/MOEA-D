//float DecomposedCost(individual *ind, float[] z, double[] lambda);

float DecomposedCost(individual *ind, float z[], double lambda[])
{
	float fx[maxfun];
	float *cost_ptr;
	double fabs_fx_z[maxfun];
	cost_ptr = &(ind->fitness[0]);
	for(int i = 0; i < nfunc; i++)
	{
		cost_ptr = &(ind->fitness[i]);
		fx[i] = *cost_ptr;
	}
	for(int j = 0; j < nfunc; j++)
	{
		fabs_fx_z[j] = lambda[j] * fabs(fx[j] - z[j]);
	}
	float maximum = fabs_fx_z[0];
	for(int k = 0; k < nfunc; k++)
	{
		if(fabs_fx_z[k] > maximum)
			maximum = fabs_fx_z[k];
	}
	return maximum;
}