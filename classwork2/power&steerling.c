#include <stdio.h>
#include <math.h>

double npower(double num, unsigned long long pow)
{
	if (pow == 0)
		return 1;
	if (pow == 1)
		return num;
	if (pow == 2)
		return num * num;
	
	if (pow & 1) 
		return num * npower(num, pow - 1);
	return npower(npower(num, pow >> 1), 2);
}

double betternpower(double num, unsigned long long pow)
{
	double res = 1;
	while (pow) {
		if (pow & 1) {
			res *= num;
			--pow;
		}
		else {
			num *= num;
			pow >>= 1;
		}
	}
	return res;
}

double steerling(unsigned long long n)
{
	return 	sqrtf(2 * M_PI * (double)n) *
			betternpower((double)n / M_E, n);
}

unsigned long long fact(unsigned long long n)
{
	unsigned long long retval = 1;
	for (unsigned long long i = 1; i < n + 1; ++i)
		retval *= i;
	return retval;
}

int main() 
{
	for (unsigned long long n = 1; n < 20; ++n) {
		double s = steerling(n);
		unsigned long long f = fact(n);
		double d = (double)f - s;
		printf("%llu - %lg = %lg\n", f, s, d);
		if (fabs(d) < (double)f * 0.01)
			printf("-----%llu-----\n", n);
	}
}
