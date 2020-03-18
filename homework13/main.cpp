#include <iostream>
#include <fstream>
#include <cinttypes>
#include <iomanip>

using namespace std;

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	// sign = 1 so number is negative.
	// to get smallest module we need
	// to get smallest fraction part
	// and get smallest exponential part.
	// so we make number unnormalized.
	// exp = 0, so exp = -126, smallest possible
	// to get nonzero number, we set fraction part
	// to one, so frac = 0.0...01, smallest possible
	uint32_t num = (1UL << 31) + 1UL;
	
	out 	<< fixed << setprecision(80)
			<< *reinterpret_cast<float*>(&num) << endl;
	
	out	<< hex << num << endl;
	
	in.close();
	out.close();
}
