#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>

//#include <stdio.h>

using namespace std;

template<typename T>
class vec
{
private:
	T* data;
	int size;
public:
	vec(int n)
	{
		size = n;
		data = new T[size];
	}
	
	T operator*(const vec<T>& other)
	{
		T retval = 0;
		for (int i = 0; i < size; ++i) 
			retval += other.data[i] * data[i];
		return retval;
	}
	
	float len()
	{
		return sqrtf((*this) * (*this));
	}
	
	void normalize()
	{
		float l = len();
		if (l == 0) return;
		for (int i = 0; i < size; ++i)
			data[i] /= l;
	}
	
	template<typename U>
	friend istream& operator>>(istream& in, vec<U>& v);
	
	template<typename U>
	friend ostream& operator<<(ostream& out, vec<U>& v);
	
	~vec()
	{
		delete[] data;
	}
};

template<typename T>
float angle(vec<T>& a, vec<T>& b)
{
	a.normalize();
	b.normalize();
	
	return acos(a * b);
}

template<typename T>
istream& operator>>(istream& in, vec<T>& v)
{
	for (int i = 0; i < v.size; ++i)
		in >> v.data[i];
	return in;
}

template<typename T>
ostream& operator<<(ostream& out, vec<T>& v)
{
	for (int i = 0; i < v.size; ++i)
		out << v.data[i] << " ";
	return out;
}

int main()
{	
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	int n = 0;
	in >> n;
	
	cout << n << endl;
	
	if (n <= 0)
		return 0;
	
	vec<float> v1(n);
	vec<float> v2(n);
	in>>v1>>v2;
	
	cout<<v1<<endl<<v2<<endl;
	
	string c;
	in>>c;
	cout<<c<<endl;
	
	out<<setprecision(2)<<fixed;
	if (c == "dot")
		out<<(v1 * v2);
	if (c == "angle")
		out<<angle(v1, v2);
		
	in.close();
	out.close();
}
