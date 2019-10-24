#include <iostream>
#include <cstdio>

using namespace std;

typedef unsigned long long ull;

ull C_combinations(ull n, ull k) {
	if (k == 0) 
		return 1;
		
	if (n == 0)
		return 0;
		
	if (n == k)
		return 1;
		
	if (n < k)
		return 0;
	
	if (2 * k > n) 
		return C_combinations(n, n - k);
	
	ull dynamic[k + 1] = {0}, temporary[k + 1];
	dynamic[0] = temporary[0] = 1;
	
	ull *dynp = dynamic;
	ull *tmpp = temporary;
	
	while (n--) {
		for (ull i = 1; i < k + 1; ++i)
			tmpp[i] = dynp[i] + dynp[i - 1];
		swap(tmpp, dynp);
	}
	
	return dynamic[k];
}

ull A_combinations(ull n, ull k) {
	if (n < k)
		return 0;
	
	ull retval = 1;
	
	for (ull i = n - k + 1; i < n + 1; ++i)
		retval *= i;
	
	return retval;
}

int main() {
	freopen("input.txt", "r", stdin);
	
	char choice = 0;
	ull n = 0;
	ull k = 0;
	
	cin >> choice >> n >> k;
	cout << choice << " " << n << " " << k;
	
	freopen("output.txt", "w", stdout);
	
	if (choice == 'A')
		cout << A_combinations(n, k);
	else
		cout << C_combinations(n, k);
}
