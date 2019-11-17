#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

double secr(long long n, long long depth)
{
	if (n == depth) return 1.;
	return 2. + 3. / secr(n, depth + 1);
}

double seci(long long n)
{
	double r = 1.;
	for (int i = 0; i < n; ++i)
		r = 2. + 3. / r;
	
	return r;
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	long long n = 0;
	in >> n;
	
	out << setprecision(3) << fixed;
	out << secr(n, 0) << endl;
	
	in.close();
	out.close();
}
