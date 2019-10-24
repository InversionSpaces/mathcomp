#include <stdio.h>
#include <stdlib.h>

double And(double a, double b)
{
	return a * b;
}

double Or(double a, double b)
{
	return a + b - a * b;
}

double Not(double a)
{
	return 1 - a;
}


int main() {
	stdout = fopen("output.txt", "wt");
	stdin = fopen("input.txt", "rt");

	double a, b;
	char op;
	int count = scanf("%lg%c%lg", &a, &op, &b);
	if(!count) scanf("%c%lg", &op, &a);
	double ans = 0;
	if(op == '+')
		ans = Or(a, b);
	else if(op == '*')
		ans = And(a, b);
	else
		ans = Not(a);
	printf("%.2lf", ans);
}
