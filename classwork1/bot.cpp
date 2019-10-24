#include <stdio.h>
#include <stdlib.h>

int move(int n) {
	if (n == 1) return 1;
	if (n <= 4) return n - 1;
	
	n -= 5;
	n = n % 4;
	if (n == 0) return rand() % 4;
	return n;
}

int main() {
	int n;
	printf("Initial:");
	scanf("%d", &n);
	
	while (1) {
		int tmp;
		printf("n: %d. Your turn:", n);
		scanf("%d", &tmp);
		if (tmp < 1 || tmp > 3) {
			printf("Not by rules\n");
			return 1;
		}
		n -= tmp;
		if (n == 0) {
			printf("I win\n");
			return 0;
		}
		n -= move(n);
		if (n == 0) {
			printf("You win\n");
			return 0;
		}
	}
}
