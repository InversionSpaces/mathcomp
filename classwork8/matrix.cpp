#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>
#include <cinttypes>
#include <cmath>
#include <algorithm>

using namespace std;

inline uint32_t greater2power(uint32_t val)
{
	if (val == 0) return 1;
	
	val--;
	val |= val >> 1;
	val |= val >> 2;
	val |= val >> 4;
	val |= val >> 8;
	val |= val >> 16;
	val++;
	
	return val;
}

class matrix
{
private:
	vector<float> data;
public:
	const size_t sizex;
	const size_t sizey;
	
	matrix(const matrix& m) :
		sizex(m.sizex),
		sizey(m.sizey),
		data(m.data)
	{
	}

	matrix(size_t size) :
		sizex(size),
		sizey(size),
		data(size * size)
	{
	}
	
	matrix(size_t sizex, size_t sizey) : 
		sizex(sizex),
		sizey(sizey),
		data(sizex * sizey)
	{
	}
	
	void set(float val)
	{
		fill(data.begin(), data.end(), val);
	}
	
	void set(const matrix& other, size_t sx, size_t sy)
	{
		assert(other.sizex + sx <= sizex);
		assert(other.sizey + sy <= sizey);
		
		for (int i = 0; i < other.sizex; ++i)
			for (int j = 0; j < other.sizex; ++j)
				data[(i + sx) * sizey + j + sy] = other[i][j];
	}
	
	float* operator[](const size_t& i)
	{
		return data.data() + i * sizey;
	}
	
	const float* operator[](const size_t& i) const
	{
		return data.data() + i * sizey;
	}
	
	matrix submatrix(size_t is, size_t ie, size_t js, size_t je) const
	{
		assert(is < ie);
		assert(js < je);
		assert(ie <= sizex);
		assert(je <= sizey);
		
		size_t nx = ie - is;
		size_t ny = je - js;
		
		matrix retval(nx, ny);
		for (size_t i = 0; i < nx; ++i)
			for (size_t j = 0; j < ny; ++j)
				retval[i][j] = data[j + js + (i + is) * sizey];
		
		return retval;
	}
	
	const matrix& operator=(const matrix& other)
	{
		assert(sizex == other.sizex);
		assert(sizey == other.sizey);
		
		set(other, 0, 0);
		
		return *this;
	}
	
	friend matrix operator/(const matrix& b, const float& a);
	
	friend matrix operator*(const matrix& b, const float& a);
	friend matrix operator*(const float& a, const matrix& b);
	
	friend matrix operator*(const matrix& a, const matrix& b);
	friend matrix operator-(const matrix& a, const matrix& b);
	friend matrix operator+(const matrix& a, const matrix& b);
	
	matrix transposed() const
	{
		matrix retval(sizey, sizex);
		for (size_t i = 0; i < sizex; ++i)
			for (size_t j = 0; j < sizey; ++j)
				retval[j][i] = data[i * sizey + j];
		
		return retval;
	}
	
	friend istream& operator>>(istream& in, matrix& m);
	friend ostream& operator<<(ostream& out, const matrix& m);
	
	~matrix() = default;
};

matrix operator/(const matrix& b, const float& a)
{
	matrix retval(b.sizex, b.sizey);
	for (size_t i = 0; i < retval.sizex * retval.sizey; ++i)
		retval.data[i] = b.data[i] / a;
	
	return retval;
}

matrix operator*(const float& a, const matrix& b)
{
	matrix retval(b.sizex, b.sizey);
	for (size_t i = 0; i < retval.sizex * retval.sizey; ++i)
		retval.data[i] = a * b.data[i];
		
	return retval;
}

matrix operator*(const matrix& b, const float& a)
{
	matrix retval(b.sizex, b.sizey);
	for (size_t i = 0; i < retval.sizex * retval.sizey; ++i)
		retval.data[i] = a * b.data[i];
	
	return retval;
}

matrix operator*(const matrix& a, const matrix& b)
{
	assert(a.sizey == b.sizex);
	
	matrix retval(a.sizex, b.sizey);
	for (size_t i = 0; i < a.sizex; ++i)
		for (size_t j = 0; j < b.sizey; ++j) {
			float tmp = 0;
			for (size_t k = 0; k < b.sizex; ++k)
				tmp += a[i][k] * b[k][j];
			retval[i][j] = tmp;
		}
	
	return retval;
}

matrix operator+(const matrix& a, const matrix& b)
{
	assert(a.sizex == b.sizex);
	assert(a.sizey == b.sizey);
	
	matrix retval(a.sizex, a.sizey);
	for (size_t i = 0; i < a.sizex; ++i)
		for (size_t j = 0; j < a.sizey; ++j)
			retval[i][j] = a[i][j] + b[i][j];
	
	return retval;
}

matrix operator-(const matrix& a, const matrix& b)
{
	assert(a.sizex == b.sizex);
	assert(a.sizey == b.sizey);
	
	matrix retval(a.sizex, a.sizey);
	for (size_t i = 0; i < a.sizex; ++i)
		for (size_t j = 0; j < a.sizey; ++j)
			retval[i][j] = a[i][j] - b[i][j];
	
	return retval;
}

istream& operator>>(istream& in, matrix& m)
{
	for (size_t i = 0; i < m.sizex * m.sizey; ++i)
		in >> m.data[i];
		
	return in;
}

ostream& operator<<(ostream& out, const matrix& m)
{
	for (size_t i = 0; i < m.sizex; ++i) {
		for (size_t j = 0; j < m.sizey; ++j)
			out << m[i][j] << " ";
		out << endl;
	}
	
	return out;
}

typedef struct
{
	matrix L, U;
} LU;

LU LUfact(const matrix& A)
{
	assert(A.sizex == A.sizey);
	
	if (A.sizex == 1)
	{
		matrix L(1, 1);
		matrix U(1, 1);
		
		L[0][0] = 1;
		U[0][0] = A[0][0];
		
		return {L, U};
	}
	
	float a = A[0][0];
	matrix v = A.submatrix(1, A.sizex, 0, 1) / a;
	matrix w = A.submatrix(0, 1, 1, A.sizey);
	matrix B = A.submatrix(1, A.sizex, 1, A.sizey);
	
	LU lu = LUfact(B - v * w);
	
	matrix L(A.sizex, A.sizey);
	matrix U(A.sizex, A.sizey);
	
	L[0][0] = 1;
	U[0][0] = a;
	
	for (size_t i = 1; i < A.sizex; ++i) {
		L[0][i] = U[i][0] = 0;
		L[i][0] = v[i - 1][0];
		U[0][i] = w[0][i - 1]; 
	}
	
	for (size_t i = 1; i < A.sizex; ++i)
		for (size_t j = 1; j < A.sizey; ++j) {
			L[i][j] = lu.L[i - 1][j - 1];
			U[i][j] = lu.U[i - 1][j - 1];
		}
		
	return {L, U};
}

matrix __fast_mul(const matrix& _a, const matrix& _b)
{
	assert(_a.sizex == _a.sizey);
	assert(_b.sizex == _b.sizey);
	assert(_a.sizex == _b.sizex);
	
	size_t size = _a.sizex;
	
	if (size <= 2)
		return _a * _b;
	
	matrix a11 = _a.submatrix(0, size / 2, 0, size / 2);
	matrix a12 = _a.submatrix(0, size / 2, size / 2, size);
	matrix a21 = _a.submatrix(size / 2, size, 0, size / 2);
	matrix a22 = _a.submatrix(size / 2, size, size / 2, size);
	
	matrix b11 = _b.submatrix(0, size / 2, 0, size / 2);
	matrix b12 = _b.submatrix(0, size / 2, size / 2, size);
	matrix b21 = _b.submatrix(size / 2, size, 0, size / 2);
	matrix b22 = _b.submatrix(size / 2, size, size / 2, size);
	
	matrix p1 = __fast_mul(a11 + a22, b11 + b22);
	matrix p2 = __fast_mul(a21 + a22, b11);
	matrix p3 = __fast_mul(a11, b12 - b22);
	matrix p4 = __fast_mul(a22, b21 - b11);
	matrix p5 = __fast_mul(a11 + a12, b22);
	matrix p6 = __fast_mul(a21 - a11, b11 + b12);
	matrix p7 = __fast_mul(a12 - a22, b21 + b22);
	
	matrix c11 = p1 + p4 - p5 + p7;
	matrix c12 = p3 + p5;
	matrix c21 = p2 + p4;
	matrix c22 = p1 - p2 + p3 + p6;
	
	matrix retval(size, size);
	
	retval.set(c11, 0, 0);
	retval.set(c12, 0, size / 2);
	retval.set(c21, size / 2, 0);
	retval.set(c22, size / 2, size / 2);
	
	return retval;
}

matrix fast_mul(const matrix& a, const matrix& b)
{
	assert(a.sizey == b.sizex);
	
	int size = max(	
			max(greater2power(a.sizex), greater2power(a.sizey)),
			max(greater2power(b.sizex), greater2power(b.sizey))
		);
	
	matrix _a(size), _b(size);
	
	_a.set(0.f); 
	_b.set(0.f);
	
	_a.set(a, 0, 0);
	_b.set(b, 0, 0);
	
	return __fast_mul(_a, _b);
}

matrix upreverse(const matrix& a)
{
	assert(a.sizex == a.sizey);
	
	matrix retval(a.sizex, a.sizey);
	
	for (int i = 0; i < a.sizex; ++i)
		for (int j = 0; j < a.sizey; ++j) {
			if (i == j) {
				retval[i][j] = 1.f / a[i][j];
				continue;
			}
			
			if (j < i) {
				retval[i][j] = 0.f;
				continue;
			}
			
			float tmp = 0.f;
			for (int k = i; k < j; ++k)
				tmp += retval[i][k] * a[k][j];
			
			retval[i][j] = -tmp / a[j][j];
		}
		
	return retval;
}

matrix reverse(const matrix& a)
{
	assert(a.sizex == a.sizey);
	
	LU lu = LUfact(a);
	
	/*
	cout << "L:" << endl;
	cout << lu.L;
	cout << "U:" << endl;
	cout << lu.U;
	cout << "L * U:" << endl;
	cout << lu.L * lu.U;
	*/
	
	const matrix ur = upreverse(lu.U);
	
	/*
	cout << "UR:" << endl;
	cout << ur;
	cout << "UR * U:" << endl;
	cout << ur * lu.U;
	*/
	
	const matrix lt = lu.L.transposed();
	const matrix ltr = upreverse(lt);
	const matrix lr = ltr.transposed();
	
	/*
	cout << "LR:" << endl;
	cout << lr;
	cout << "LR * L:" << endl;
	cout << lr * lu.L;
	
	cout << "UR * LR:" << endl;
	cout << (ur * lr);
	*/
	
	return (ur * lr);
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	int n = 0;
	in >> n;
	
	matrix m1(n, n);
	
	in >> m1;
	//cout << m1;
	
	matrix t = reverse(m1);
	
	//cout << t;	
	
	out << setprecision(1) << fixed;
	out << t;
	//cout << t * m1;
	//cout << endl << t * m1;
	/*
	LU lu = LUfact(m1);
	
	cout << lu.L << endl << lu.U;
	*/
	in.close();
	out.close();
}
