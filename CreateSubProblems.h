#include <time.h>
#include <math.h>

void CreateSubProblems(subproblems *sp_tmp_ptr);

double norm(double x[]);

double CalculateDistance(double x1, double y1, double x2, double y2);

//double[] sort(double[] a[]);

void CreateSubProblems(subproblems sp_tmp[])
{
	double lamb[nfunc];
	double A[popsize], B[popsize], C[popsize];
	double D[popsize][popsize];
	double SO[popsize][popsize];
	for(int i = 0; i < popsize; i++)
	{
		for(int j = 0; j < nfunc; j++)
		{
			srand( time(NULL) );
			/* 產生 [0, 1) 的浮點數亂數 */
			lamb[j] = (double) rand() / (RAND_MAX + 1.0);
		}
		double n = norm(lamb);
		for(int k = 0; k < nfunc; k++)
		{
			lamb[k] = lamb[k] / n;
		}
		for(int l = 0; l < nfunc; l++)
		{
			sp_tmp[i].lambda[l] = lamb[l];
		}
	}
	for(int i = 0; i < popsize; i++)
	{
		for(int a = 0; a < nfunc; a++)
		{
			A[a] = sp_tmp[i].lambda[a];
		}
		for(int j = 0; j < popsize; j++)
		{
			for(int b = 0; b < nfunc; b++)
			{
				B[b] = sp_tmp[i].lambda[b];
			}
			D[i][j] = CalculateDistance(A[0], A[1], B[0], B[1]);
		}
	}
	for(int i = 0; i < popsize; i++)
	{
		double tmp[popsize];
		for(int c = 0; c < popsize; c++)
		{
			C[c] = D[i][c];
			tmp[c] = C[c];
		}
		for(int x = 1; x < popsize; x++)
		{
			if (tmp[x] > tmp[i])                //Comparing other array elements
			{
				int tmp_1 = tmp[i];         //Using temporary variable for storing last value
				tmp[i] = tmp[x];            //replacing value
				tmp[x] = tmp_1;             //storing last value
			}
		}
		for(int k = 0; k < popsize; k++)
		{
			SO[i][k] = tmp[k];
		}
		for(int j = 0; j < T; j++)
		{
			sp_tmp[i].Neighbors[j] = SO[i][j];
		}
		//sp_ptr = sp_ptr->next;
	}
	return;
}

double norm(double x[])
{
	double m_sum = 0.0;
	for (int i=0; i < nfunc; i++)
	{
		m_sum  += pow(x[i], 2);
	}
	return sqrt(m_sum);
}

double CalculateDistance(double x1, double y1, double x2, double y2)
{
	double diffx = x1 - x2;
	double diffy = y1 - y2;
	double diffx_sqr = pow(diffx,2);
	double diffy_sqr = pow(diffy,2);
	double distance = sqrt(diffx_sqr + diffy_sqr);
	
	return distance;
}