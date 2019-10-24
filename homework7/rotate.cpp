#include <iostream>
#include <fstream>
#include <array>
#include <cmath>
#include <string>
#include <iomanip>

using namespace std;

typedef array<float, 3> vec3f;

float len(const vec3f& v)
{
	float l = 0;
	for (int i = 0; i < 3; ++i)
		l += v[i] * v[i];
	return sqrtf(l);
}

inline void set_len(vec3f& v, float l)
{
	float old = len(v);
	for (int i = 0; i < 3; ++i)
		v[i] *= l / old;
}

vec3f operator*(float a, const vec3f& v)
{
	vec3f retval = v;
	for (int i = 0; i < 3; ++i)
		retval[i] *= a;
	return retval;
}

inline void rotate(vec3f& v, float a, int d)
{
	float c = cos(a);
	float s = sin(a);
	
	float x = v[0];
	float y = v[1];
	float z = v[2];
	
	if (d == 0) {
		v[1] = c * y - s * z;
		v[2] = s * y + c * z;
	}
	if (d == 1) {
		v[0] = c * x + s * z;
		v[2] = c * z - s * x;
	}
	if (d == 2) {
		v[0] = c * x - s * y;
		v[1] = s * x + c * y;
	}
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	string s;
	
	vec3f v = {0, 1, 0};
	float l = 1;
	
	while (in >> s) {
		if (s == "len") {
			float tmp;
			in >> tmp;
			l += tmp;
		}
		if (s == "rot") {
			string d;
			float a;
			in >> d >> a;
			rotate(v, a, d == "X" ? 0 : d == "Y" ? 1 : 2);
		}
	}
	
	set_len(v, l);
	
	out << setprecision(1) << fixed;
	for (int i = 0; i < 3; ++i) {
		if (fabs(v[i]) < 0.01) {
			v[i] = 0.0;
		}
		out << v[i] << " ";
	}
	
	in.close();
	out.close();
}

