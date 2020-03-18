#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <bitset>
#include <cinttypes>

using namespace std;

int main()
{
	//cout << sizeof(float) << endl;
	
	//ifstream in("input.txt");
	//ofstream out("output.txt");
	
	float num = 0;
	cin >> num;
	
	bitset<sizeof(uint32_t) * 8> bin1(
		*reinterpret_cast<uint32_t*>(&num)
	);
	
	cout << bin1 << endl;
	
	*reinterpret_cast<uint32_t*>(&num) += 100UL << 23;
	
	cout << fixed << setprecision(3) << num << endl;
	
	bitset<sizeof(uint32_t) * 8> bin2(
		*reinterpret_cast<uint32_t*>(&num)
	);
	
	cout << bin2 << endl;
	
	//in.close();
	//out.close();
}
