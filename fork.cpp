#include <unistd.h>
#include <sys/types.h>

int main() {
	while(1)
		fork();
}
