#include <iostream>
#include <stdio.h>


using namespace std;

double root2(size_t deep) 
{
	double value = 1.;
	while (deep--) value = 1. / (2. + value);
	return (value + 1);
}

int main() 
{
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

	size_t d = 0;
	cin >> d;
	cout << root2(d);
}
