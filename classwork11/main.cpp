#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

typedef unsigned long long set_t;

inline void set(set_t& s, size_t i)
{
	s |= 1ULL << i;
}

inline void reset(set_t& s, size_t i)
{
	s &= !(1ULL << i);
}

inline int get(const set_t& s, size_t i)
{
	return !(!(s & (1ULL << i)));
}

inline int next(set_t& s, const int count)
{
	s++;
	return s < (1 << count);
}

double count(const set_t& s, const double vals[], int n)
{
	double retval = 0;
	
	for (int i = 0; i < n; ++i)
		if (get(s, i)) retval += vals[i];
		
	return retval;
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	int n = 0;
	double w = 0;
	in >> n >> w;
	
	double vals[2][n];
	for (int i = 0; i < n; ++i)
		in >> vals[0][i] >> vals[1][i];
	
	double best_s = 0;
	double best_w = 0;
	
	set_t s = 0;
	do {
		double tmp_w = count(s, vals[0], n);
		if (tmp_w > w) continue;
		
		double tmp_s = count(s, vals[1], n);
		if (tmp_s > best_s) {
			best_s = tmp_s;
			best_w = tmp_w;
		}
	} while (next(s, n));
	
	out << fixed << setprecision(2) << best_w << " " << best_s << endl;
	
	in.close();
	out.close();
}
