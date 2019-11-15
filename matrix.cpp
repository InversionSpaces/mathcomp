#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>
#include <cinttypes>
#include <cmath>
#include <algorithm>
#include <type_traits>

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

template<typename T>
class matrix
{
private:
	vector<T> data;
public:
	const size_t sizex;
	const size_t sizey;
	
	matrix(const matrix<T>& m) :
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
	
	void set(T val)
	{
		fill(data.begin(), data.end(), val);
	}
	
	void set(const matrix<T>& other, size_t sx, size_t sy)
	{
		assert(other.sizex + sx <= sizex);
		assert(other.sizey + sy <= sizey);
		
		for (int i = 0; i < other.sizex; ++i)
			for (int j = 0; j < other.sizex; ++j)
				data[(i + sx) * sizey + j + sy] = other[i][j];
	}
	
	T* operator[](const size_t& i)
	{
		return data.data() + i * sizey;
	}
	
	const T* operator[](const size_t& i) const
	{
		return data.data() + i * sizey;
	}
	
	matrix<T> submatrix(size_t is, size_t ie, size_t js, size_t je) const
	{
		assert(is < ie);
		assert(js < je);
		assert(ie <= sizex);
		assert(je <= sizey);
		
		size_t nx = ie - is;
		size_t ny = je - js;
		
		matrix<T> retval(nx, ny);
		for (size_t i = 0; i < nx; ++i)
			for (size_t j = 0; j < ny; ++j)
				retval[i][j] = data[j + js + (i + is) * sizey];
		
		return retval;
	}
	
	const matrix<T>& operator=(const matrix<T>& other)
	{
		assert(sizex == other.sizex);
		assert(sizey == other.sizey);
		
		set(other, 0, 0);
		
		return *this;
	}
	
	template<typename U>
	friend matrix<typename std::common_type<T, U>::type> 
		operator/(const matrix<T>& b, const typename std::enable_if
			<std::is_arithmetic<U>::value, U>::type& a);
		
	template<typename U>
	friend matrix<typename std::common_type<T, U>::type> 
		operator*(const matrix<T>& b, typename std::enable_if
			<std::is_arithmetic<U>::value, U>::type& a);
		
	template<typename U>
	friend matrix<typename std::common_type<T, U>::type> 
		operator*(const typename std::enable_if
			<std::is_arithmetic<U>::value, U>::type& a, 
			const matrix<T>& b);
	
	template<typename U>
	friend matrix<typename std::common_type<T, U>::type> 
		operator*(const matrix<T>& a, const matrix<U>& b);
		
	template<typename U>
	friend matrix<typename std::common_type<T, U>::type> 
		operator-(const matrix<T>& a, const matrix<U>& b);
		
	template<typename U>
	friend matrix<typename std::common_type<T, U>::type> 
		operator+(const matrix<T>& a, const matrix<U>& b);
	
	matrix<T> transposed() const
	{
		matrix<T> retval(sizey, sizex);
		for (size_t i = 0; i < sizex; ++i)
			for (size_t j = 0; j < sizey; ++j)
				retval[j][i] = data[i * sizey + j];
		
		return retval;
	}
	
	template<typename U>
	friend istream& operator>>(istream& in, matrix<U>& m);
	
	template<typename U>
	friend ostream& operator<<(ostream& out, const matrix<U>& m);
	
	~matrix() = default;
};

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type> 
	operator/(const matrix<T>& b, const U& a)
{
	matrix<typename std::common_type<T, U>::type> 
		retval(b.sizex, b.sizey);
	for (size_t i = 0; i < retval.sizex * retval.sizey; ++i)
		retval.data[i] = b.data[i] / a;
	
	return retval;
}

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type>
	operator*(const U& a, const matrix<T>& b)
{
	matrix<typename std::common_type<T, U>::type>
		retval(b.sizex, b.sizey);
	for (size_t i = 0; i < retval.sizex * retval.sizey; ++i)
		retval.data[i] = a * b.data[i];
		
	return retval;
}

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type>
	operator*(const matrix<T>& b, const U& a)
{
	matrix<typename std::common_type<T, U>::type>
		retval(b.sizex, b.sizey);
	for (size_t i = 0; i < retval.sizex * retval.sizey; ++i)
		retval.data[i] = a * b.data[i];
	
	return retval;
}

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type> 
	operator*(const matrix<T>& a, const matrix<U>& b)
{
	assert(a.sizey == b.sizex);
	
	matrix<typename std::common_type<T, U>::type>  
		retval(a.sizex, b.sizey);
	for (size_t i = 0; i < a.sizex; ++i)
		for (size_t j = 0; j < b.sizey; ++j) {
			typename std::common_type<T, U>::type tmp = 0;
			for (size_t k = 0; k < b.sizex; ++k)
				tmp += a[i][k] * b[k][j];
			retval[i][j] = tmp;
		}
	
	return retval;
}

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type>
	operator+(const matrix<T>& a, const matrix<U>& b)
{
	assert(a.sizex == b.sizex);
	assert(a.sizey == b.sizey);
	
	matrix<typename std::common_type<T, U>::type> 
		retval(a.sizex, a.sizey);
	for (size_t i = 0; i < a.sizex; ++i)
		for (size_t j = 0; j < a.sizey; ++j)
			retval[i][j] = a[i][j] + b[i][j];
	
	return retval;
}

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type> 
	operator-(const matrix<T>& a, const matrix<U>& b)
{
	assert(a.sizex == b.sizex);
	assert(a.sizey == b.sizey);
	
	matrix<typename std::common_type<T, U>::type>
		retval(a.sizex, a.sizey);
	for (size_t i = 0; i < a.sizex; ++i)
		for (size_t j = 0; j < a.sizey; ++j)
			retval[i][j] = a[i][j] - b[i][j];
	
	return retval;
}

template<typename T>
istream& operator>>(istream& in, matrix<T>& m)
{
	for (size_t i = 0; i < m.sizex * m.sizey; ++i)
		in >> m.data[i];
		
	return in;
}

template<typename T>
ostream& operator<<(ostream& out, const matrix<T>& m)
{
	for (size_t i = 0; i < m.sizex; ++i) {
		for (size_t j = 0; j < m.sizey; ++j)
			out << m[i][j] << " ";
		out << endl;
	}
	
	return out;
}

template<typename T>
struct LU
{
	matrix<T> L, U;
};

template<typename T>
LU<T> LUfact(const matrix<T>& A)
{
	assert(A.sizex == A.sizey);
	
	if (A.sizex == 1)
	{
		matrix<T> L(1, 1);
		matrix<T> U(1, 1);
		
		L[0][0] = 1;
		U[0][0] = A[0][0];
		
		return {L, U};
	}
	
	T a = A[0][0];
	matrix<T> v = A.submatrix(1, A.sizex, 0, 1) / a;
	matrix<T> w = A.submatrix(0, 1, 1, A.sizey);
	matrix<T> B = A.submatrix(1, A.sizex, 1, A.sizey);
	
	LU<T> lu = LUfact(B - v * w);
	
	matrix<T> L(A.sizex, A.sizey);
	matrix<T> U(A.sizex, A.sizey);
	
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

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type> 
	__fast_mul(const matrix<T>& _a, const matrix<U>& _b)
{
	using C = typename std::common_type<T, U>::type;
	
	assert(_a.sizex == _a.sizey);
	assert(_b.sizex == _b.sizey);
	assert(_a.sizex == _b.sizex);
	
	size_t size = _a.sizex;
	
	if (size <= 2)
		return _a * _b;
	
	matrix<T> a11 = _a.submatrix(0, size / 2, 0, size / 2);
	matrix<T> a12 = _a.submatrix(0, size / 2, size / 2, size);
	matrix<T> a21 = _a.submatrix(size / 2, size, 0, size / 2);
	matrix<T> a22 = _a.submatrix(size / 2, size, size / 2, size);
	
	matrix<U> b11 = _b.submatrix(0, size / 2, 0, size / 2);
	matrix<U> b12 = _b.submatrix(0, size / 2, size / 2, size);
	matrix<U> b21 = _b.submatrix(size / 2, size, 0, size / 2);
	matrix<U> b22 = _b.submatrix(size / 2, size, size / 2, size);
	
	matrix<C> p1 = __fast_mul(a11 + a22, b11 + b22);
	matrix<C> p2 = __fast_mul(a21 + a22, b11);
	matrix<C> p3 = __fast_mul(a11, b12 - b22);
	matrix<C> p4 = __fast_mul(a22, b21 - b11);
	matrix<C> p5 = __fast_mul(a11 + a12, b22);
	matrix<C> p6 = __fast_mul(a21 - a11, b11 + b12);
	matrix<C> p7 = __fast_mul(a12 - a22, b21 + b22);
	
	matrix<C> c11 = p1 + p4 - p5 + p7;
	matrix<C> c12 = p3 + p5;
	matrix<C> c21 = p2 + p4;
	matrix<C> c22 = p1 - p2 + p3 + p6;
	
	matrix<C> retval(size, size);
	
	retval.set(c11, 0, 0);
	retval.set(c12, 0, size / 2);
	retval.set(c21, size / 2, 0);
	retval.set(c22, size / 2, size / 2);
	
	return retval;
}

template<typename T, typename U>
matrix<typename std::common_type<T, U>::type> 
	fast_mul(const matrix<T>& a, const matrix<U>& b)
{
	assert(a.sizey == b.sizex);
	
	int size = max(	
			max(greater2power(a.sizex), greater2power(a.sizey)),
			max(greater2power(b.sizex), greater2power(b.sizey))
		);
	
	matrix<T> _a(size);
	matrix<U> _b(size);
	
	_a.set(0); 
	_b.set(0);
	
	_a.set(a, 0, 0);
	_b.set(b, 0, 0);
	
	return __fast_mul(_a, _b);
}

template<typename T>
matrix<T> upreverse(const matrix<T>& a)
{
	assert(a.sizex == a.sizey);
	
	matrix<T> retval(a.sizex, a.sizey);
	
	for (int i = 0; i < a.sizex; ++i)
		for (int j = 0; j < a.sizey; ++j) {
			if (i == j) {
				retval[i][j] = 1. / a[i][j];
				continue;
			}
			
			if (j < i) {
				retval[i][j] = 0.;
				continue;
			}
			
			float tmp = 0.;
			for (int k = i; k < j; ++k)
				tmp += retval[i][k] * a[k][j];
			
			retval[i][j] = -tmp / a[j][j];
		}
		
	return retval;
}

template<typename T>
matrix<T> reverse(const matrix<T>& a)
{
	assert(a.sizex == a.sizey);
	
	LU<T> lu = LUfact(a);
	
	/*
	cout << "L:" << endl;
	cout << lu.L;
	cout << "U:" << endl;
	cout << lu.U;
	cout << "L * U:" << endl;
	cout << lu.L * lu.U;
	*/
	
	const matrix<T> ur = upreverse(lu.U);
	
	/*
	cout << "UR:" << endl;
	cout << ur;
	cout << "UR * U:" << endl;
	cout << ur * lu.U;
	*/
	
	const matrix<T> lt = lu.L.transposed();
	const matrix<T> ltr = upreverse(lt);
	const matrix<T> lr = ltr.transposed();
	
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

template<typename T>
matrix<T> adj_matrix(matrix<T> adj)
{
	assert(adj.sizex == adj.sizey);
	
	for (int i = 0; i < adj.sizex; ++i)
		for (int j = 0; j < adj.sizex; ++j)
			for (int k = 0; k < adj.sizex; ++k)
				adj[j][k] = adj[j][k] + adj[j][i] * adj[i][k];

	return adj;
}

void kirf()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	int n = 0;
	in >> n;
	
	matrix<int> adj(n);
	
	in >> adj;
	
	matrix<int> power(n);
	power.set(0);
	
	for (int i = 0; i < n; ++i) {
		int tmp = 0;
		
		for (int j = 0; j < n; ++j) {
			if (i == j) continue;
			
			tmp += lround(adj[i][j]);
		}
		
		power[i][i] = tmp;
		adj[i][i] = 0.f;
	}
	
	out << n << endl;
	//out << (power - adj);
	out << adj * power;
	
	in.close();
	out.close();
}

int main()
{
}
