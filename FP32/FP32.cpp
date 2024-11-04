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
	uint32_t getmantissa() const noexcept {
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
		if (qbits >= 32 || a == 0) return 0; // a == 0 bad
		if (qbits == 0) return a;
		//return roundDiv32v2(a, qbits); //
		int32_t unit;
		if (a >= 0) unit = 1;
		else unit = -1;
		//                              '>' NOT '>=' mathematically for '-'
		int32_t tmp = (a % (1 << qbits)) >= (unit * (1 << (qbits - 1)));
		if ((a % (1 << qbits)) == 0) tmp = 0; // very bad
		return (a >> qbits) + tmp;
	}
	int32_t roundDiv32v2(int32_t a, int32_t qbits) const noexcept {
		a >>= qbits;
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
	operator double() const noexcept {
		uint32_t tmp = data;
		return double(reinterpret_cast<float*>(&tmp)[0]);
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
		if (l.data == 0x80000000 && r.data == 0x80000000) { // bad
			return l;
		}

		int32_t el = l.getexp(), er = r.getexp(), eres;
		int32_t ml = l.getmantissa() + (int32_t(l.getexp() > 0) << 23), mr = r.getmantissa() + (int32_t(r.getexp() > 0) << 23), mres;
		if (l.getsign()) ml = -ml;
		if (r.getsign()) mr = -mr;

		// sum of a normal and of a subnormal

		if (el > er) { // here
			if (er == 0) // subnormals
				++er;
			eres = el;
			mres = ml + roundDiv32(mr, el - er);
			//cout << mres;
			//if (el - er <= 32) mres = ml + mr / (1 << (el - er)); //why?
			//	if ((ml % (1 << (el - er)) >= )
			//else mres = ml;
			//mres = roundDiv32((ml << (el - er)) + mr, el - er);
		}

		else if (er > el) {
			if (el == 0) // subnormals
				++el;
			eres = er;
			//mres = mr + roundDiv32(ml, er - el);
			cout << mres << endl;
			//mres = roundDiv32((mr << (er - el)) + ml, er - el);
		}
		else {
			eres = el;
			mres = ml + mr;
		}
		

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
		if (mres >= (int32_t(1) << 24)) { //instruction to count 00001***mant zeros can be used //one operation! // just an if
			if (eres > 0) { // subnormals
				//mres >>= 1; // it seems like I should keep one more bit to increase precision. I should rewrite all the code
				mres = roundDiv32(mres, 1);
			}
			++eres;
			//mres = roundDiv(mres, 1); //works correct with simple div
		}
			
		while (mres < (int32_t(1) << 23) && eres > 0) { // mb error is here or it is not
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

	// masks:
	// sign: 0x80000000
	// exp:  0x7F800000
	// mant: 0x007FFFFF

	static uint32_t add3(uint32_t l, uint32_t r, float& example) noexcept {
		example = float(FP32(l)) + float(FP32(r)); //
		uint32_t res;
		uint32_t el = (l & 0x7F800000) >> 23;
		uint32_t er = (r & 0x7F800000) >> 23;
		uint32_t eres;
		int32_t ml = (l & 0x007FFFFF + ((el > 0) << 23)) << 1; // one extra bit
		int32_t mr = (r & 0x007FFFFF + ((er > 0) << 23)) << 1;
		int32_t mres;
		if (l >> 31) ml = -ml; // how to make it faster?
		if (r >> 31) mr = -mr;

		if (el == 0x7F800000 || er == 0x7F800000) { // nan and inf checking
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return l;
			if (mr != 0 && er == 0x7F800000) // if right is nan
				return r;
			if (l == r)
				return l;
			else
				return r + 1; // return nan if inf - inf
		}

		if (el > er) { // calulate exponent and mantissa making exponents equal each other
			er += (er == 0);
			eres = el;
			if (el - er >= 32) mres = ml;
			else mres = ml + (mr >> (el - er));
		}
		else {
			el += (el == 0);
			eres = er;
			if (er - el >= 32) mres = mr;
			else mres = mr + (ml >> (er - el));
		}

		if (mres < 0) { // calculate mantissa and a sign
			res = 0x80000000;
			mres = -mres;
		}
		else {
			res = 0;
		}

		if (mres >= 0x02000000) { // if mres is greater than a 2^23 (2^24)
			mres >>= (eres > 0); // subnormals
			++eres;
		}
		while (mres < 0x01000000 && eres > 0) { // if mres is less than a 2^23 (2^24)
			--eres;
			mres <<= (eres > 0); // subnormals
		}

		mres = (mres + (mres & 0x00000001)) >> 1; // last bit stored. How to include it in mres >= 2^25 section??? 
		if (mres >= 0x01000000) { // mres is greater than 2^23
			mres >>= (eres > 0); // subnormals
			++eres;
		}

		if (eres >= 0xFF) {
			return res + 0x7F800000;
		}
		res += eres << 23;
		if (eres > 0)
			return res + mres - 0x00800000;
		else
			return res + mres;

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
			return res;
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
			//mres = roundDiv64(mres, 1); //works correct with simple div // it will break my while
		}

		if (eres > 0) {
			mres -= uint64_t(1) << 23; //as normals
			res.data += mres;
			res.data += eres << 23;
		}
		else if (eres >= -23) {
			res.data += mres/* >> -eres*/;
			res.data >>= -eres + 1; 
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

	// masks:
	// sign: 0x80000000
	// exp:  0x7F800000
	// mant: 0x007FFFFF
	
	static uint32_t mul3(uint32_t l, uint32_t r, float& example) noexcept {
		example = float(FP32(l)) * float(FP32(r)); //
		uint32_t res = (l ^ r) & 0x80000000;
		uint32_t el = l & 0x7F800000;
		uint32_t er = r & 0x7F800000;
		int32_t eres;
		uint32_t ml = l & 0x007FFFFF;
		uint32_t mr = r & 0x007FFFFF;
		uint32_t mres;

		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res | l;
			else if (er == 0x7F800000) // if right is inf or nan
				return res | r;
			return res | l;
		}

		eres = ((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x00800000 * (el > 0)) >> 11) * ((mr + 0x00800000 * (er > 0)) >> 11); // calculating 1 bit extended mantissa

		while ((mres < 0x01000000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1; 
			mres <<= 1;
		}
		while (mres >= 0x02000000) { // mres is greater than 2^23 (2^25). Can it be speeded up?
			eres += 1;
			mres >>= 1;
		}

		if (eres > 0) { // if normal
			mres = (mres + (mres & 0x00000001)) >> 1; // last bit stored
			if (mres >= 0x01000000) { // mres is greater than 2^23 (2^24)
				eres += 1;
				mres >>= 1;
			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F800000;
			}
			else {
				return res + mres - 0x00800000 + (eres << 23); // calculating result
			}
		}
		else if (eres >= -23) { // if subnormal
			mres >>= (-eres + 1);
			return res + ((mres + (mres & 0x00000001)) >> 1); // last bit stored
		}

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
		return *this;
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
		return *this;
	}
	FP32& operator*= (const FP32& bf) noexcept {
		return *this;
	}
	FP32& operator/= (const FP32& bf) noexcept {
		return *this;
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
	public:
	FP32 l;
	FP32 r;
	Alltests() = default;
	bool equal_prec(uint32_t a, uint32_t b, uint32_t prec) {
		//if (a == 0 && b == 0x80000000 || a == 0x80000000 && b == 0) return true;
		for (uint32_t p = 0; p <= prec; ++p) {
			if (a - p == b || a + p == b) return true;
		}
		return false;
	}
	bool run_specific() {
		vector<uint32_t> vl = {0x380fa};
		vector<uint32_t> vr = {0x17f1570};
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
			if (/*!isnan(l.example)*/ (l.example == l.example) && !equal_prec(FP32(l.example).data, l.data, 0)) {
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
		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 76542) {
		//#pragma omp parallel for
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 76542) {
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
					cout << hex << endl << lc << ", " << rc << " , that is " << FP32(uint32_t(lc)).example << ", " << FP32(uint32_t(rc)).example;
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
	void run2() {
		int input;
		uint64_t lc, rc;
		uint32_t res;
		float f;
		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 69632) { // 6528
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 69632) {
				res = FP32::add3(uint32_t(lc), uint32_t(rc), f);
				if (f == f && res != FP32(f).data && res - 1 != FP32(f).data && res + 1 != FP32(f).data) {
				//if (f == f && res != FP32(f).data) { // 2940 3641 4341 -> 701 700
					cout << hex << endl << lc << ", " << rc << " , that is " << FP32(uint32_t(lc)).example << ", " << FP32(uint32_t(rc)).example << " ERROR\n";
					cout << f << " expected, " << float(FP32(res)) << " instead\n";
					FP32(f).print();
					cout << bitset<32>(res) << endl;
					cout << endl << "Continue? 1 - yes, 0 - no\n";
					cin >> input;
					if (input) continue;
					return;
				}
			}
		}		
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
	//if (flag) flag = tests.run();
	tests.run2();
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