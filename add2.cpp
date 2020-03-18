#include <stdio.h>

int main() {
	FILE* in = fopen("test.txt", "rt");
	unsigned long long a, b;
	fscanf(in, "%llu %llu", &a, &b);
	fclose(in);

	FILE* out = fopen("test.out", "wt");
	fprintf(out, "%llu", a + b);
	fclose(out);
}
