#include <vector>

using namespace std;

template<typename T>
class Test
{
	/*
	friend void foo(Test<T> bar);
	*/
	
	template<typename S, typename U>
	friend void baz(Test<S> foo, Test<U> bee);
};

template<typename T, typename U>
void baz(Test<T> foo, Test<U> bee)
{
	//test
}

int main() {
	Test<int> a;
	Test<int> b;
	
	baz(a, b);
}
