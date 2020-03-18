#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <fstream>

using namespace std;

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	int n = 0;
	in >> n;
	
	vector<pair<string, char> > codes(n);
	for (auto &i: codes)
		in >> i.first >> i.second;
	
	string code;
	in >> code;
	
	int pos = 0;
	while (pos < code.size()) {
		for (const auto& c: codes) {
			if (code.compare(pos, c.first.size(), c.first) == 0) {
				out << c.second;
				pos += c.first.size();
				continue;
			}
		}
	}
	
	in.close();
	out.close();
}
