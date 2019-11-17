#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

typedef unsigned long long ull;

ull _rec(ull n)
{
	return ull(pow(9, n)) + ull(pow(3, n)) + 1;
}

ull rec(ull n)
{
	if (n == 0) return 3;
	if (n == 1) return 13;
	if (n == 2) return 91;
	return 13UL * rec(n - 1) - 39UL * rec(n - 2) + 27UL * rec(n - 3);
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	ull n = 0;
	in >> n;
	out << rec(n) << endl;
	
	in.close();
	out.close();
}
