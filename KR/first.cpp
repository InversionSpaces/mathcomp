#include <iostream>
#include <fstream>

using namespace std;

typedef unsigned long long ull;

ull nim(ull n)
{
	if (n <= 100UL) return n - 1;
	if (n % 100UL == 0) return 98;
	return n % 100UL - 1;
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	ull n = 0;
	in >> n;
	cout << n;
	out << nim(n);
	
	in.close();
	out.close();
}
