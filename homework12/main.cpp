#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cstring>

using namespace std;

int main() {
	FILE* in = fopen("input.txt", "r");
	FILE* out = fopen("output.txt", "w");
	
	const int size = 1 << (sizeof(char) * CHAR_BIT);
	
	char* str = 0;
	size_t n = 0;
	
	getline(&str, &n, in);
	
	int counter = 0;
	int count[size] = {};
	
	n = strlen(str); // -1 for \n from getline // UPD NO MORE
	
	for (int i = 0; i < n; ++i) {
		if (!count[str[i]]) counter++;
		count[str[i]]++;
	}
	
	double probability[counter];
	for (int i = 0, k = 0; i < size; ++i) {
		if (count[i]) {
			probability[k++] = double(count[i]) / double(n);
		}
	}
	
	double entropy = 0.;
	for (int i = 0; i < counter; i++) {
		entropy -= probability[i] * log2l(probability[i]);
	}
	
	// I don't understand it
	fprintf(out, "%0.2lf\n", entropy > 1e-9 ? 8.f / entropy : 8.f); 
	
	free(str);
	
	fclose(out);
	fclose(in);
}
