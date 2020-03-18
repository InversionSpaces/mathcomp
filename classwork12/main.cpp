#include <fstream>
#include <iomanip>

using namespace std;

int main() {
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	int n = 0;
	in >> n;
	
	double res = 1;
	for (int i = 1; i < n; ++i) {
		res *= (365 - i) / 365.;
	}
	
	out << setprecision(3) << fixed << 1 - res;
	
	in.close();
	out.close();
}
