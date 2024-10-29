#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
//#include <omp.h>
#include <iomanip>
#include <random>
#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;

class FP32 {
	uint32_t getsign() const noexcept {
		return data >> 31;
	}
	int32_t getexp() const noexcept {
		return (data << 1) >> 24;
	}
	uint64_t getmantissa() const noexcept {
		return data & 0x007FFFFF;
	}
	bool isfpinf() const noexcept {
		return ((getexp() == 0xFF) && getmantissa() == 0);
	}
	bool isfpnan() const noexcept {
		return ((getexp() == 0xFF) && getmantissa() != 0);
	}
	uint64_t roundDiv64(uint64_t a, uint64_t qbits) const noexcept {
		uint64_t tmp = (a % (uint64_t(1) << qbits)) >= (uint64_t(1) << (qbits - 1));
		return (a >> qbits) + tmp;
	}
	int32_t roundDiv32(int32_t a, int32_t qbits) const noexcept {
		if (qbits >= 32) return 0;
		int32_t unit;
		if (a >= 0) unit = 1;
		else unit = -1;
		//                              '>' NOT '>=' mathematically
		int32_t tmp = (a % (1 << qbits)) >= (unit << (qbits - 1));
		return (a >> qbits) + tmp;
	}
	int32_t myabs(int32_t a) {
		if (a > 0) return -a;
		return a;
	}

public:
	float example;
	uint32_t data;
	FP32() = default;
	FP32(float f) noexcept {
		data = reinterpret_cast<uint32_t*>(&f)[0];
		example = f;
	}
	FP32(const FP32& bf) noexcept {
		example = bf.example;
		data = bf.data;
	}
	FP32& operator= (const FP32& bf) noexcept {
		example = bf.example;
		data = bf.data;
		return *this;
	}
	FP32& operator= (const uint32_t& bf) noexcept {
		data = bf;
		example = *this;
		return *this;
	}
	FP32(double d) noexcept {
		example = static_cast<float>(d);
		data = reinterpret_cast<uint32_t*>(&example)[0];
	}
	operator float() const noexcept {
		uint32_t tmp = data;
		return reinterpret_cast<float*>(&tmp)[0];
	}
	FP32(uint32_t t) noexcept {
		data = t;
		example = *this;
	}
	void print() const noexcept {
		cout << bitset<32>(data) << endl;
	}

	FP32 add2(const FP32 l, const FP32 r) noexcept { //Add of normals is correct
		FP32 res;
		res.example = l.example + r.example;
		uint32_t sign;
		FP32 fpinf;
		FP32 fpnan;
		fpinf.data = 0x7f800000 + (l.getsign() << 31);
		fpnan.data = 0x7f800001;
		fpinf.example = res.example;
		fpnan.example = res.example;
		if (l.isfpnan() || r.isfpnan()) {
			return fpnan;
		}
		if (l.isfpinf() || r.isfpinf()) {
			if (l.getsign() ^ r.getsign()) return fpnan;
			else return fpinf;
		}
		// 0 +- 0 interactions?

		int32_t el = l.getexp(), er = r.getexp(), eres;
		int32_t ml = l.getmantissa() + (int32_t(l.getexp() > 0) << 23), mr = r.getmantissa() + (int32_t(r.getexp() > 0) << 23), mres;
		//          (l.getmantissa() + (uint64_t(l.getexp() > 0) << 23))    (r.getmantissa() + (uint64_t(r.getexp() > 0) << 23)); //bad type (m1*m2)/2^23 = m1/2^23 * m2
		if (l.getsign()) ml = -ml;
		if (r.getsign()) mr = -mr;

		if (el > er) {
			eres = el;
			mres = ml + roundDiv32(mr, el - er);
		}

		else if (er > el) {
			eres = er;
			//cout << dec << "2 " << mr << " " << ml << " " << roundDiv32(ml, er - el) << endl;
			mres = mr + roundDiv32(ml, er - el);
			//cout << mres << endl;
		}
		else {
			eres = el;
			mres = ml + mr;
		}
		//передавать в roundDiv число - разность, а потом отбрасывать деление на числа с длиной битов меньше 1? Что это? Кто это написал? Аа, идея ясна но она фигня

		if (mres < 0) {
			res.data = int32_t(1) << 31;
			mres = -mres;
		}
		else /*if (mres > 0)*/ {
			res.data = 0;
		}
		//else {
		//	mres = 1; // sad..
		//}
		if (eres == 0 && mres >= (int32_t(1) << 23)) { // from subnormals to normals //
			cout << "1" << endl;
			//cout << hex << mres << endl;
			++eres;
			//mres -= (int32_t(1) << 23);
			//mres >>= 1;
			//mres += (int32_t(1) << 23);
		}

		while (mres >= (int32_t(1) << 24)) { //instruction to count 00001***mant zeros can be used //one operation! // just an if
			cout << "2" << endl;
			++eres;
			//if (eres > 0) // non correct
			mres >>= 1;
			//mres = roundDiv(mres, 1); //works correct with simple div
		}

		while (mres < (int32_t(1) << 23) && eres > 0) { // just an if
			cout << "3" << endl;
			--eres;
			if (eres > 0) //subnormals
				mres <<= 1;
		}
		
		if (eres >= 0xFF) {
			res.data += 0xFF << 23;
			return res;
		}
		res.data += eres << 23;
		if (eres > 0)
			res.data += mres - (int32_t(1) << 23);
		else
			res.data += mres;

		return res;
	}

	FP32 mul2(const FP32 l, const FP32 r) const noexcept { //last bit error in subnormals
		FP32 res; 
		//res.data = (l.getsign() ^ r.getsign()) << 31;
		res.data = 0;
		res.example = l.example * r.example;
		FP32 fpinf;
		FP32 fpnan;
		fpinf.data = 0x7f800000 + ((l.getsign() ^ r.getsign()) << 31);
		fpnan.data = 0x7f800001 + ((l.getsign() ^ r.getsign()) << 31);
		fpinf.example = res.example;
		//if (l.getsign() ^ r.getsign()) fpinf.example = -INFINITY;
		fpnan.example = res.example;

		if (l.isfpinf() || r.isfpinf()) {
			return fpinf;
		}
		if (l.isfpnan() || r.isfpnan()) {
			return fpnan;
		}
		if ((l.data >> 1) == 0 || (r.data >> 1) == 0) {
			res.data = (l.getsign() ^ r.getsign()) << 31; //bad
			return data;
		}

		uint64_t mres = (l.getmantissa() + (uint64_t(l.getexp() > 0) << 23)) * (r.getmantissa() + (uint64_t(r.getexp() > 0) << 23)); //bad type (m1*m2)/2^23 = m1/2^23 * m2
		int32_t eres = l.getexp() + r.getexp() - int32_t(127) + (l.getexp() == 0 || r.getexp() == 0); // subnormal
		while ((mres < (uint64_t(1) << 46)) && (eres >= 0)) { // from subnormal to normal ZEEEEEEEE00000*****... if one arg is subnormal and m3 < 2^23
			--eres;
			mres <<= 1;
		}
		mres = roundDiv64(mres, 23);
		
		while (mres >= (uint64_t(1) << 24)) { //instruction to count 00001***mant zeros can be used. mb only 1 shift?
		// just an if
			++eres;
			mres >>= 1;
			//mres = roundDiv(mres, 1); //works correct with simple div
		}

		if (eres > 0) {
			mres -= uint64_t(1) << 23; //as normals
			res.data += mres;
			res.data += eres << 23;
		}
		else if (eres >= -23) {
			res.data += mres/* >> -eres*/;
			res.data >>= -(eres - 1); 
			//res.data = roundDiv(res.data, -eres); //? correct? NO
		}
		else {
			res.data += (l.getsign() ^ r.getsign()) << 31; //bad
			return res;
		}
		if (res.data >> 31 || res.isfpnan()) {
			return fpinf;
		}
		res.data += (l.getsign() ^ r.getsign()) << 31;
		return res;
	}

	FP32& operator+= (const FP32& bf) noexcept {
		uint16_t sign;
		if (getsign() && bf.getsign()) sign = 1;
		else if (getsign()) {
			sign = uint16_t(1) << 15;
			*this -= bf;
			data |= sign;
			return *this;
		}
		else if (bf.getsign()) {
			*this -= bf;
			return *this;
		}
	}
	FP32& operator-= (const FP32& bf) noexcept {
		uint16_t sign;
		if (getsign() && bf.getsign()) sign = 1;
		else if (getsign()) {
			sign = uint16_t(1) << 15;
			*this -= bf;
			data |= sign;
			return *this;
		}
		else if (bf.getsign()) {
			*this -= bf;
			return *this;
		}
	}
	FP32& operator*= (const FP32& bf) noexcept {

	}
	FP32& operator/= (const FP32& bf) noexcept {

	}
	FP32 operator+ (const FP32& bf) const noexcept {
		FP32 res = *this;
		res += bf;
		return res;
	}
	FP32 operator- (const FP32& bf) const noexcept {
		FP32 res = *this;
		res -= bf;
		return res;
	}
	FP32 operator* (const FP32& bf) const noexcept {
		FP32 res = *this;
		res *= bf;
		return res;
	}
	FP32 operator/ (const FP32& bf) const noexcept {
		FP32 res = *this;
		res /= bf;
		return res;
	}
};

class Alltests {
	FP32 l;
	FP32 r;
public:
	Alltests() = default;
	bool equal_prec(uint32_t a, uint32_t b, uint32_t prec) {
		//if (a == 0 && b == 0x80000000 || a == 0x80000000 && b == 0) return true;
		for (uint32_t p = 0; p <= prec; ++p) {
			if (a - p == b || a + p == b) return true;
		}
		return false;
	}
	bool run_specific() {
		vector<uint32_t> vl = {};
		vector<uint32_t> vr = {};
		size_t from = 0;
		for (size_t i = from; i < vl.size(); ++i) {
			cout << hex << vl[i] << ", " << vr[i];
			l = vl[i];
			r = vr[i];
			//l = l.add(l, r);
			//l = l.mul2(l, r);
			l = l.add2(l, r);
			//cout << endl << float(BF16(l)) << " " << float(BF16(r)) << endl;
			//cout << l.example << " " << float(l) << endl;
			if (/*!isnan(l.example)*/ (l.example == l.example) && !equal_prec(FP32(l.example).data, l.data, 1)) {
				//cout << " ADD ERROR\n";
				cout << " ERROR\n";
				cout << l.example << " expected, " << float(l) << " instead\n";
				FP32(l.example).print();
				l.print();
				return false;
			}
			cout << " PASSED\n";
		}
		return true;
	}
	bool run() {
		uint64_t lc, rc;
		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 76543) {
		//#pragma omp parallel for
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 76543) {
				//cout << hex << lc << ", " << rc;
				//if (lc < 0x00800000 || rc < 0x00800000) continue;
				l = uint32_t(lc);
				r = uint32_t(rc);
				l = l.add2(l, r);
				//l = l.mul2(l, r);
				//cout << endl << float(BF16(l)) << " " << float(BF16(r)) << endl;
				//cout << l.example << " " << float(l) << endl;
				if (/*!isnan(l.example)*/ (l.example == l.example) && !equal_prec(FP32(l.example).data, l.data, 1)) {
					//cout << " ADD ERROR\n";
					cout << hex << endl << lc << ", " << rc;
					cout << " ERROR\n";
					cout << l.example << " expected, " << float(l) << " instead\n";
					FP32(l.example).print();
					l.print();
					return false;
				}
				//cout << " PASSED\n";
			}
		}
		return true;
	}
};

template <typename T>
vector <vector <T> > mmul(const vector <vector <T> >& A, const vector <vector <T> >& B) {
	size_t m = A.size(), n = B.size(), p = B[0].size(); // A: m x n, B: n x p
	vector <vector <T> > C (m, vector<T>(p)); 

	for (size_t i = 0; i < m; ++i) {

	for (size_t j = 0; j < p; ++j) {

	for (size_t k = 0; k < n; ++k) {

	C[i][j] += A[i][k] * B[k][j];

	}
	}
	}

	return C;
}

vector <vector <double> > mmulv2(const vector <vector <double> >& A, const vector <vector <double> >& B) {
	size_t m = A.size(), n = B.size(), p = B[0].size(); // A: m x n, B: n x p
	vector <vector <double> > C(m, vector<double>(p));

	for (size_t i = 0; i < m; ++i) {
		for (size_t k = 0; k < n; ++k) {
//#pragma omp parallel for shedul static 6
			for (size_t j = 0; j < p; ++j) {


				C[i][j] += A[i][k] * B[k][j];

			}
		}
	}

	return C;
}

int main() {
	//using mat = vector <vector <double> >;
	//size_t n = 10000;
	//auto l = [](const mat& m) {
	//	for (size_t i = 0; i < m.size(); ++i) {
	//		for (size_t j = 0; j < m[i].size(); ++j) {
	//			cout << setw(10) << m[i][j];
	//		}
	//		cout << endl;
	//	}
	//};
	//auto r = [&n]() {
	//	mat m = mat(n, vector<double>(n));
	//	//srand(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
	//	//srand(5);
	//	//default_random_engine generator;
	//	random_device rd;
	//	mt19937 mt(rd());
	//	uniform_real_distribution<double> dist(-100.0, 100.0);
	//	for (size_t i = 0; i < n; ++i) {
	//		for (size_t j = 0; j < n; ++j) {
	//			m[i][j] = dist(mt);
	//		}
	//	}
	//	return m;
	//};
	//
	//mat A = {
	//	{1, 5, 6},
	//	{6, 3, 8},
	//	{2, 2, 2},
	//	{1, 0, -4}
	//};
	//mat B = {
	//	{3, -8, 0, 0},
	//	{-2, -3, -3, 1},
	//	{12, 11, 6, -1}
	//};
	//
	//A = r();
	//B = r();
	//
	//cout << "ready\n";
	//auto start = chrono::high_resolution_clock::now();
	//
	//mat C = mmulv2(A, B);
	//
	//auto stop = chrono::high_resolution_clock::now();
	//auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	//cout << duration.count() << endl;
	//mat C = mmul(B, A);
	//l(C);


	bool flag = true;
	Alltests tests;
	//if (flag) flag = tests.run_specific();
	if (flag) flag = tests.run();
	cout << endl << "ENDED" << endl;


	//int a = -100;
	//cout << a << endl;
	//cout << (a >> 1) << endl;
	//cout << (a >> 2) << endl;
	//cout << (a >> 3) << endl;
	//cout << (a >> 4) << endl;
	//cout << (-99 % 16) << endl;
	//cout << (-99 >> 4) << endl;


	return 0;
}