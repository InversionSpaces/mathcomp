#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <iomanip>

using namespace std;

struct vec2f
{
	float x, y;
};

struct vec3f
{
	float a, b, c;
};

float operator*(const vec2f& a, const vec2f& b)
{
	return a.x * b.x + a.y * b.y;
}

vec2f operator-(const vec2f& a, const vec2f& b)
{
	return {a.x - b.x, a.y - b.y};
}

float len2(const vec2f& a)
{
	return a * a;
}

typedef pair<vec2f, vec2f> entry;

void print(const vector<entry>& entries)
{
	for (int i = 0; i < entries.size(); ++i) {
		cout << entries[i].first.x << " " << entries[i].first.y << endl;
		cout << entries[i].second.x << " " << entries[i].second.y << endl;
		cout << "*" << endl;
	}
	cout << "-----" << endl;
}

inline vec2f func(const vec2f& data, const vec3f& koef)
{
	float x = 1 - (data.x / koef.a) * (data.x / koef.a) - 
					(data.y / koef.b) * (data.y / koef.b);
	float y =  koef.c * data.x / data.y;
	
	return {x, y};
}

inline vector<entry> read_file()
{
	ifstream in("input.txt");
	
	vector<entry> retval;
	
	float x, y, X, Y;
	while (in >> x >> y >> X >> Y) {
		vec2f in = {x, y};
		vec2f out = {X, Y};
		
		retval.push_back(make_pair(in, out));
	}
	
	return retval;
}

inline float count_rmse(vector<entry> entries)
{
	float retval = 0;
	for (int i = 0; i < entries.size(); ++i) {
		entry e = entries[i];
		retval += len2(e.first - e.second);
	}
	
	return sqrtf(retval);
}

inline void update_entries(vector<entry>& entries, 
	const vector<vec2f>& inputs, const vec3f& koefs)
{
	for (int i = 0; i < inputs.size(); ++i)
		entries[i].first = func(inputs[i], koefs);
}

inline vector<vec2f> get_inputs(const vector<entry>& entries)
{
	vector<vec2f> retval(entries.size());
	for (int i = 0; i < retval.size(); ++i)
		retval[i] = entries[i].first;
	
	return retval;
}

vec3f optimize(vector<entry> entries)
{
	vector<vec2f> inputs = get_inputs(entries);
	
	vec3f koefs = {8, 8, 0.8};
	update_entries(entries, inputs, koefs);
	float rmse = count_rmse(entries);
	
	for (float a = 8; a < 12; a += 0.1)
		for (float b = 8; b < 12; b += 0.1)
			for (float c = 0.8; c < 1.2; c += 0.1) {
				vec3f _koefs = {a, b, c};
				update_entries(entries, inputs, _koefs);
				float _rmse = count_rmse(entries);
				
				//cout << _rmse << endl;
				
				if (_rmse < rmse) {
					rmse = _rmse;
					koefs = _koefs;
				}
			}
			
	return koefs;
}

int main()
{
	vector<entry> entries = read_file();
	
	vec3f koefs = optimize(entries);
	ofstream out("output.txt");
	
	out << setprecision(1) << fixed;
	out << koefs.a << " " << koefs.b << " " << koefs.c << endl;
	
	out.close();
}

