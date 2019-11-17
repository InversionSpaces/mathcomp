#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <iomanip>

using namespace std;

typedef unsigned long long ull;

class BigInteger
{
private:
	vector<int> parts;
	int sign;
public:
	static const int ndigits = 2;
	static const int base = static_cast<int>(pow(10, ndigits)); 

	friend istream& operator>>(istream& in, BigInteger& bint);
	friend ostream& operator<<(ostream& out, const BigInteger& bint);
	
	friend int operator<(const BigInteger& a, const BigInteger& b);
	friend int operator>(const BigInteger& a, const BigInteger& b);
	friend int operator==(const BigInteger& a, const BigInteger& b);
	
	friend BigInteger operator*(const BigInteger& a, const BigInteger& b);
	friend BigInteger operator+(const BigInteger& a, const BigInteger& b);
	friend BigInteger operator-(const BigInteger& a, const BigInteger& b);
	
	friend BigInteger operator-(const BigInteger& a);
};

int ceil(int a, int b)
{
	return (a / b) + ((a % b) > 0);
}

int powi(int n, int p)
{
	if (p == 0) return 1;
	while (--p) n *= n;
	return n;
}

int operator==(const BigInteger& a, const BigInteger& b)
{
	if (a.parts.size() < b.parts.size()) return 0;
	if (a.parts.size() > b.parts.size()) return 0;
	
	int size = a.parts.size();
	
	for (int i = 0; i < size; i++)
		if (a.parts[i] != b.parts[i]) return 0;
	
	return 1;
}

int operator<(const BigInteger& a, const BigInteger& b)
{
	if (a.parts.size() < b.parts.size()) return 1;
	if (a.parts.size() > b.parts.size()) return 0;
	
	int size = a.parts.size();
	
	for (int i = size - 1; i >= 0; i++) {
		if (a.parts[i] < b.parts[i]) return 1;
		if (a.parts[i] > b.parts[i]) return 0;
	}
	
	return 0;
}

int operator>(const BigInteger& a, const BigInteger& b)
{
	return (b < a);
}

BigInteger operator*(const BigInteger& a, const BigInteger& b)
{		
	BigInteger retval;
	retval.sign = 0;
	retval.parts.push_back(0);
	
	BigInteger _a = a;
	BigInteger _b = b;
	
	cout << "A " << _a << endl;
	cout << "B " << _b << endl;
	
	for (int i = 0; i < _b.parts.size(); ++i) {
		int part = _b.parts[i];
		
		BigInteger mul;
		mul.sign = 0;
		
		int tmp = 0;
		for (int j = 0; j < _a.parts.size(); ++j) {
			int mpart = part * _a.parts[j] + tmp;
			tmp = mpart / mul.base;
			mul.parts.push_back(mpart % mul.base);
		}
		
		while (tmp) {
			mul.parts.push_back(tmp);
			tmp /= mul.base;
		}
		
		cout << mul << endl;
		
		for (int j = 0; j < mul.parts.size(); ++j) {
			if (j + i < retval.parts.size())
				retval.parts[j + i] += mul.parts[j];
			else
				retval.parts.push_back(mul.parts[j]);
		}
		
		tmp = 0;
		for (int j = 0; j < retval.parts.size(); ++j) {
			int rpart = retval.parts[j] + tmp;
			retval.parts[j] = rpart % retval.base;
			tmp = rpart / retval.base;
		}
		
		while (tmp) {
			retval.parts.push_back(tmp);
			tmp /= retval.base;
		}
	}
	
	retval.sign = a.sign ^ b.sign;
	
	return retval;
}

BigInteger operator-(const BigInteger& a)
{
	BigInteger retval = a;
	
	retval.sign = !a.sign;
	
	return retval;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b)
{
	BigInteger _a = a;
	BigInteger _b = b;
	
	if (_b.sign) {
		_b.sign = 0;
		
		return _a + _b;
	}
	
	if (_a.sign) {
		_a.sign = 0;
		
		return -(_a + _b);
	}
	
	BigInteger retval;
	retval.sign = 0;
	
	if (_a == _b) {
		retval.parts.push_back(0);
		
		return retval;
	}
	
	if (_a < _b) {
		return -(_b - _a);
	}
	
	for (int i = 0; i < _b.parts.size(); ++i)
		_a.parts[i] -= _b.parts[i];
	
	int len = 0;
	for (int i = 0; i < _a.parts.size(); ++i) {
		if (_a.parts[i] < 0) {
			_a.parts[i] += retval.base;
			_a.parts[i + 1] -= 1;
		}
		
		if (_a.parts[i]) len = i + 1;
	}
	
	for (int i = 0; i < len; ++i)
		retval.parts.push_back(_a.parts[i]);
	
	return retval;
}

BigInteger operator+(const BigInteger& a, const BigInteger& b)
{
	BigInteger _a = a;
	BigInteger _b = b;
	
	if (_a.sign && _b.sign) {
		_a.sign = _b.sign = 0;
		
		return -(_a + _b);
	}
	
	if (_a.sign) {
		_a.sign = 0;
		
		return _b - _a;
	}
	
	if (_b.sign) {
		_b.sign = 0;
		
		return _a - _b;
	}
	
	BigInteger retval;
	
	retval.sign = 0;
	
	int equal = min(a.parts.size(), b.parts.size());
	
	int tmp = 0;
	for (int i = 0; i < equal; ++i) {
		int part = (a.parts[i] + b.parts[i] + tmp);
		retval.parts.push_back(part % retval.base);
		tmp = part / retval.base;
	}
	
	for (int i = equal; i < a.parts.size(); ++i) {
		int part = (a.parts[i] + tmp);
		retval.parts.push_back(part % retval.base);
		tmp = part / retval.base;
	}
	
	for (int i = equal; i < b.parts.size(); ++i) {
		int part = (b.parts[i] + tmp);
		retval.parts.push_back(part % retval.base);
		tmp = part / retval.base;
	}
	
	if (tmp) retval.parts.push_back(tmp);
	
	return retval;
	
}

ostream& operator<<(ostream& out, const BigInteger& bint)
{
	if (bint.parts.size() == 0) {
		out << "<empty int>";
		
		return out;
	}
	
	if (bint.sign) out << '-';
	
	out << bint.parts[bint.parts.size() - 1];
	out << setfill('0') << setw(bint.ndigits);
	for (int i = bint.parts.size() - 2; i >= 0; --i)
		out << bint.parts[i];
	
	out << setfill(' ') << setw(0);
	
	return out;
}

istream& operator>>(istream& in, BigInteger& bint)
{
	string s;
	in >> s;
	
	bint.parts.clear();
	bint.sign = s[0] == '-';
	
	int size = s.size() - bint.sign;
	int pos = s.size();
	
	for (int i = 0; i < ceil(size, bint.ndigits); ++i) {
		int tmp = 0;
		for (int j = max(bint.sign, pos - bint.ndigits); j < pos; ++j) {
			int power = powi(10, pos - j - 1);
			tmp += int(s[j] - '0') * power;
		}
		pos -= bint.ndigits;
		bint.parts.push_back(tmp);
	}
	
	return in;
}

int main()
{
	ifstream in("input.txt");
	ofstream out("output.txt");
	
	BigInteger a;
	BigInteger b;
	
	cin >> a >> b;
	cout << (a * b) << endl;
	/*
	cout << (a + b) << endl;
	cout << (a - b) << endl;
	cout << (a > b) << endl;
	cout << (a < b) << endl;
	cout << (a == b) << endl;
	*/
	
	in.close();
	out.close();
}
