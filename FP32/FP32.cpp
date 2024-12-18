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
#include <unordered_map>
#include <string>
#include <fstream>

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
		int64_t ml = ((uint64_t(l) & 0x007FFFFF) + ((uint64_t(el > 0)) << 23)) << 24; // 24 extra bit
		int64_t mr = ((uint64_t(r) & 0x007FFFFF) + ((uint64_t(er > 0)) << 23)) << 24;
		int64_t mres;
		ml *= -(2 * int32_t(l >> 31) - 1);
		mr *= -(2 * int32_t(r >> 31) - 1);

		if (el == 0xFF || er == 0xFF) { // nan and inf checking
			if (el == 0xFF) {
				if (ml != 0) return l;
				if (er == 0xFF) {
					if (mr != 0) return r;
					if ((l ^ r) == 0x80000000) return r + 1;
					else return r;
				}
				return l;
			}
			else
				return r;
		}

		if (el > er) { // calulate exponent and mantissa making exponents equal each other
			er += (er == 0);
			eres = el;
			if (el - er >= 64) mres = ml; // 25
			else mres = ml + (mr >> (el - er));
		}
		else {
			er += (er == 0);
			el += (el == 0);
			eres = er;
			if (er - el >= 64) mres = mr;
			else mres = mr + (ml >> (er - el));
		}

		res = 0x80000000 * (mres < 0);
		mres *= -(2 * (mres < 0) - 1);

		//		mres = (mres >> 24) + ((mres & 0xFF'FFFF) > 0x80'0000) + ((mres & 0x1FF'FFFF) == 0x180'0000); // last 24 bit stored 
		if (mres >= 0x1'0000'0000'0000 /*0x100'0000*/) { // if mres is greater than a 2^23 (2^48) // make a non-if code
			mres >>= (eres > 0); // subnormals
			++eres;
		}
		else {
			while (mres < 0x8000'0000'0000 /*0x80'0000*/ && eres > 0) { // if mres is less than a 2^23 (2^24) // can add mres == 0
				--eres;
				mres <<= (eres > 0); // subnormals
			}
		}

		mres = (mres >> 24) + ((mres & 0xFF'FFFF) > 0x80'0000) + ((mres & 0x1FF'FFFF) == 0x180'0000); // last 24 bit stored 
		if (mres >= 0x01000000) { // mres is greater than 2^23 // add this into >= 2^48
			mres >>= (eres > 0); // subnormals
			++eres;
		}

		if (eres >= 0xFF) { // in inf
			return res + 0x7F800000;
		}
		res += eres << 23; // calcuating result
		if (eres > 0)
			return res + mres - 0x00800000;
		else
			return res + mres;

		return res;
	}

	static uint32_t sub(uint32_t l, uint32_t r, float& example) noexcept { // check me!
		r = (r & 0x7FFF'FFFF) + (~r & 0x8000'0000);
		return FP32::add3(l, r, example);
	}

	FP32 mul_not_use(const FP32 l, const FP32 r) const noexcept { //last bit error in subnormals
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

	static uint32_t mul(uint32_t l, uint32_t r, float& example) noexcept { // CSR register, rounding to closest odd number // CRlibn // Increase bit stored!
		example = float(FP32(l)) * float(FP32(r)); //
		uint32_t res = (l ^ r) & 0x80000000;
		uint32_t el = l & 0x7F800000;
		uint32_t er = r & 0x7F800000;
		int32_t eres;
		uint64_t ml = l & 0x007FFFFF;
		uint64_t mr = r & 0x007FFFFF;
		uint64_t mres;

		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res + (l & 0x7FFFFFFF);
			else if (er == 0x7F800000) // if right is inf or nan
				return res + (r & 0x7FFFFFFF);
			return res + (l & 0x7FFFFFFF);
		}

		eres = ((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x00800000 * (el > 0))) * ((mr + 0x00800000 * (er > 0))); // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

		while ((mres < 0x4000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x8000'0000'0000) { // mres is greater than 2^23 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1;
		}

		if (eres > 0) { // if normal
			mres = (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored 
			if (mres >= 0x01000000) { // mres is greater than 2^23 (2^24) // should I?
				eres += 1;
				mres >>= 1;
			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F800000;
			}
			else {
				return res + uint32_t(mres) - 0x00800000 + (eres << 23); // calculating result
			}
		}
		else if (eres >= -23) { // if subnormal
//			mres >>= (-eres + 1);
			uint64_t shift = 1ull << (23 - eres + 1); 
//			return res + (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored

//			cout << dec << endl << eres << hex << " " << shift - 1 << " " << (shift >> 1) << " " << (shift << 1) - 1 << " " << shift + (shift >> 1) << endl;
			return res + (mres >> (23 - eres + 1)) + 
			( (mres & (shift - 1) ) > (shift >> 1) ) + 
			( ( (mres & ( (shift << 1) - 1) ) ) == (shift + (shift >> 1) ) );
		}

		return res;
	}

	static uint32_t newton_iter(uint32_t x, uint32_t d) { // x = x * (2 - r*x) or x = x + x(1 - dx) // 2 sum algorithm
		uint32_t res;
		float dummy = 0.0;
		int outCounter = 0;
//		cout << endl << ++outCounter << ": " << hex << x << " " << d << endl; // input data

		uint64_t mx = x & 0x007FFFFF;
		uint64_t md = d & 0x007FFFFF;
//		res = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(d, x, dummy), dummy), dummy);
		uint32_t etmp = ((((x & 0x7F800000) + (d & 0x7F800000)) >> 23) - 127);
		uint64_t mtmp = (mx * md + (mx << 23) + (md << 23)) << 2; // 1
//		cout << ++outCounter << ": " <<  mtmp << " " << etmp << endl; // correct

		mtmp >>= (etmp == 125);
//		mtmp -= 0x4000'0000'0000 * (etmp == 125); // mtmp + 2^46 << 1 >> 1 = 2^46; mtmp - 2^47: leading 2^46 bit << 1
		mtmp -= 0x8000'0000'0000 * (etmp == 125);
		etmp += (etmp == 125);
		mtmp = 0x2'0000'0000'0000 * (etmp == 126) + 0x1'0000'0000'0000 - mtmp; // 2 - dx; md = ...
//		mtmp = 0x1'0000'0000'0000 * (etmp == 126) + 0x8000'0000'0000 - mtmp;
//		md = mtmp; //
//		cout << ++outCounter << ": " <<  mtmp << " " << etmp << endl; // correct

//		cout << endl << hex << etmp << " " << mtmp << endl; //
		
		// x * tmp left

		// x is 24 bit len, tmp is 48 bit len (including leading bit) (49)
		// the x*tmp result is 72 bit len, should be reduced to 24 bit - 8 bit longes than 64 (9)
		// the idea is to multiply x as 40 bit len number, x*tmp result is 48+40 = 88 bit, that is 24 bit longer than 64 (39)
		// the result is close enough to 1.0, that's why it cannot be 2.0 after rounding - 1 if can be used
		// that's a lie! you should make an if to the exponent!
		// x <<= 40 - 24 = 16 (39 - 24 = 15)
		
		if (mtmp < 0x1'0000'0000'0000) { // mres if less than 2^23; witout if etmp -= (mtmp < 0x8000'0000'0000) (0x8000'0000'0000)
			etmp -= 1;
			mtmp <<= 1; //
		}
		if (mtmp >= 0x2'0000'0000'0000) { // mres is greater than 2^24; without if (0x1'0000'0000'0000)
			etmp += 1;
			mtmp >>= 1; //
		}
		
		etmp = (etmp << 23) + (x & 0x7F800000) - 0x3F800000; // eres
//		cout << ++outCounter << ": " << etmp << endl; // correct

		// mx * mtmp = (x1*2^32 + y1)*(x2*2^32 + y2) = x1x2*2^64 + y1y2 + x1y2*2^32 + y1x2*2^32 = mx*2^64 + md
		// increasing x up to 40 bit len
		mx += 0x0080'0000;
		mx <<= 16; // 17
//		cout << ++outCounter << ": " << mx << endl; //

		// counting using 2^52 mod
		uint64_t x1 = mx >> 26, y1 = mx & 0x3FF'FFFF, x2 = mtmp >> 26, y2 = mtmp & 0x3FF'FFFF;
//		cout << ++outCounter << ": " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
		md = y1 * y2 + ((x1 * y2) & 0x3FF'FFFF) + ((y1 * x2) & 0x3FF'FFFF); // lover half
		mx = x1 * x2 + ((x1 * y2) >> 26) + ((y1 * x2) >> 26) + (md >> 52); // upper half, + md / 2^52
		md &= 0x000F'FFFF'FFFF'FFFF; // % 2^52, lover half
//		cout << ++outCounter << ": " << mx << " " << md << endl; // error is here i think

		// translation into 2^64 mod
		res = mx >> 12; // res = mx / 2^12
		md += ((mx & 0x0FFF) << 52); // md = mx + mx & 2^12 (mx = mx'*2^52) I DONT LIKE THIS LINE
//		cout << ++outCounter << ": " << res << " " << md << endl;

		// rounding
		res += (md > 0x8000'0000'0000'0000) + (md == 0x8000'0000'0000'0000) * ((res & 0x1) == 0x1);
		etmp += (res >= 0x100'0000) << 23;
		res >>= (res >= 0x100'0000);
//		cout << ++outCounter << ": " << res << endl;

		// making result
		
		res = etmp + res - 0x0080'0000; // -?
//		res = etmp + res;
 //		cout << ++outCounter << ": " << res << endl;
		

//		cout << etmp << " " << mtmp << endl; //
		/*
		mtmp = (mtmp >> 24) + ((mtmp & 0xFF'FFFF) > 0x80'0000) + ((mtmp & 0x1FF'FFFF) == 0x180'0000); // last 24 bit stored
		if (mtmp >= 0x0100'0000) { // mres is greater than 2^23 (2^24) // should I?
				etmp += 1;
				mtmp >>= 1;
			}
//		cout << mtmp << endl; //

		res = uint32_t(mtmp) - 0x0080'0000 + (etmp << 23);
		cout << res << endl; //
		res = FP32::mul3(x, res, dummy);
		*/
		return res;
	}

	static uint32_t newton_iter2(uint32_t x, uint32_t d) {
		uint32_t res;
		int outCounter = 0;
		bool printflag = false;
		uint64_t shifted;
		if (printflag) cout << endl << ++outCounter << ": " << hex << x << " " << d << endl; // input data

		uint64_t mx = x & 0x007FFFFF;
		uint64_t md = d & 0x007FFFFF;
		uint32_t etmp = ((((x & 0x7F800000) + (d & 0x7F800000)) >> 23) - 127);
		uint64_t mtmp = (mx * md + (mx << 23) + (md << 23)) << 1;
		if (printflag) cout << ++outCounter << ": " <<  mtmp << " " << etmp << endl; //

		mtmp >>= (etmp == 125);
		mtmp -= 0x4000'0000'0000 * (etmp == 125); // mtmp + 2^46 << 1 >> 1 = 2^46; mtmp - 2^47: leading 2^46 bit << 1
		etmp += (etmp == 125);
		mtmp = 0x1'0000'0000'0000 * (etmp == 126) + 0x8000'0000'0000 - mtmp; //  24 extra bit, 24 usual bit

		if (printflag) cout << ++outCounter << ": " <<  mtmp << " " << etmp << endl; // 

		if (mtmp < 0x8000'0000'0000) { // mres if less than 2^23; witout if etmp -= (mtmp < 0x8000'0000'0000)
			etmp -= 1;
//		etmp -= (mtmp < 0x8000'0000'0000);
			mtmp <<= 1; //
 //		mtmp <<= (mtmp < 0x8000'0000'0000);
		}
		if (mtmp >= 0x1'0000'0000'0000) { // mres is greater than 2^24; without if
			etmp += 1; 
//		etmp += (mtmp >= 0x1'0000'0000'0000);
			mtmp >>= 1; // bad
//		mtmp >>= (mtmp >= 0x1'0000'0000'0000);
		}
		etmp = (etmp << 23) + (x & 0x7F800000) - 0x3F800000; // eres

		mx += 0x0080'0000; // mx is 24 bit long, etmp should be 40 bit long

		mtmp = (mtmp >> 8) + ((mtmp & 0x00FF) > 0x0080) + ((mtmp & 0x01FF) == 0x0180); // this one
		if (printflag) cout << ++outCounter << ": " << mtmp << " " << mx << endl; //

		mtmp = mtmp * mx;
		if (printflag) cout << ++outCounter << ": " << mtmp << endl; //

		res = (mtmp >> 39) + ((mtmp & 0x7F'FFFF'FFFF) > 0x40'0000'0000) + ((mtmp & 0x0FF'FFFF'FFFF) == 0x0C0'0000'0000);
		if (printflag) cout << ++outCounter << ": " << res << endl; //
		// overflow
		etmp += ((res >= 0x100'0000) << 23);
		res >>= (res >= 0x100'0000);
		res = etmp + res - 0x0080'0000;

		if (printflag) cout << ++outCounter << ": " << res << endl; //
		return res;
		
	}

	static uint64_t special_newton_iter_rouding(uint32_t x, uint32_t d) {
		uint64_t res;

		uint64_t mx = x & 0x007FFFFF;
		uint64_t md = d & 0x007FFFFF;
		uint32_t etmp = ((((x & 0x7F800000) + (d & 0x7F800000)) >> 23) - 127);
		uint64_t mtmp = (mx * md + (mx << 23) + (md << 23)) << 1;
		// check on mtmp == 0

		mtmp >>= (etmp == 125);
		mtmp -= 0x4000'0000'0000 * (etmp == 125); // mtmp + 2^46 << 1 >> 1 = 2^46; mtmp - 2^47: leading 2^46 bit << 1
		etmp += (etmp == 125);
		mtmp = 0x1'0000'0000'0000 * (etmp == 126) + 0x8000'0000'0000 - mtmp; // 2 - dx; 24 extra bit, 24 usual bit

		if (mtmp < 0x8000'0000'0000) { // mres if less than 2^23; witout if etmp -= (mtmp < 0x8000'0000'0000)
			etmp -= 1;
			mtmp <<= 1; //
		}
		if (mtmp >= 0x1'0000'0000'0000) { // mres is greater than 2^24; without if
			etmp += 1;
			mtmp >>= 1; // bad
		}
		etmp = (etmp << 23) + (x & 0x7F800000) - 0x3F800000; // eres

		mx += 0x0080'0000; // mx is 24 bit long, etmp should be 40 bit long
		mtmp = (mtmp >> 8) + ((mtmp & 0x00FF) > 0x0080) + ((mtmp & 0x01FF) == 0x0180); // this one
		mtmp = mtmp * mx; // x * (2-dx)

		res = (mtmp >> 39) + ((mtmp & 0x7F'FFFF'FFFF) > 0x40'0000'0000) + ((mtmp & 0x0FF'FFFF'FFFF) == 0x0C0'0000'0000); // 64 bit to 24 bit (max 2^64 - 1 to 2^24 - 1)

		// overflow
		etmp += ((res >= 0x100'0000) << 23);
		res >>= (res >= 0x100'0000);
		
		// binary calculations
		mtmp = res * (md + 0x0080'0000);
		if (mtmp < 0x8000'0000'0000) {
			uint64_t diff1 = 0x8000'0000'0000 - mtmp;
			uint64_t diff2 = (res + 1) * (md + 0x0080'0000) - 0x8000'0000'0000;
			if (diff2 < diff1) ++res;
		}

		res = etmp + res - 0x0080'0000;
		// binary +1 to mantissa
		return res;
	}

	static uint64_t special_newton_iter_no_rounding(uint32_t x, uint32_t d) {
		uint64_t res;

		uint64_t mx = x & 0x007FFFFF;
		uint64_t md = d & 0x007FFFFF;
		uint32_t etmp = ((((x & 0x7F800000) + (d & 0x7F800000)) >> 23) - 127);
		uint64_t mtmp = (mx * md + (mx << 23) + (md << 23)) << 1;

		mtmp >>= (etmp == 125);
		mtmp -= 0x4000'0000'0000 * (etmp == 125); // mtmp + 2^46 << 1 >> 1 = 2^46; mtmp - 2^47: leading 2^46 bit << 1
		etmp += (etmp == 125);
		mtmp = 0x1'0000'0000'0000 * (etmp == 126) + 0x8000'0000'0000 - mtmp; // 2 - dx; 24 extra bit, 24 usual bit

		if (mtmp < 0x8000'0000'0000) { // mres if less than 2^23; witout if etmp -= (mtmp < 0x8000'0000'0000)
			etmp -= 1;
			mtmp <<= 1; //
		}
		if (mtmp >= 0x1'0000'0000'0000) { // mres is greater than 2^24; without if
			etmp += 1;
			mtmp >>= 1; // bad
		}
		etmp = (etmp << 23) + (x & 0x7F800000) - 0x3F800000; // eres

		mx += 0x0080'0000; // mx is 24 bit long, etmp should be 40 bit long
		mtmp = (mtmp >> 8) + ((mtmp & 0x00FF) > 0x0080) + ((mtmp & 0x01FF) == 0x0180); // this one
		mtmp = mtmp * mx; // x * (2-dx)

		//		res = (mtmp >> 39) + ((mtmp & 0x7F'FFFF'FFFF) > 0x40'0000'0000) + ((mtmp & 0x0FF'FFFF'FFFF) == 0x0C0'0000'0000); // 64 bit to 24 bit (max 2^64 - 1 to 2^24 - 1)
		res = (mtmp >> 23) + ((mtmp & 0x7F'FFFF) > 0x40'0000) + ((mtmp & 0x0FF'FFFF) == 0x0C0'0000); // 64 bit to 40 bit (max 2^64 - 1 to 2^40 - 1)

		// overflow
		etmp += ((res >= 0x100'0000'0000) << 23);
		res >>= (res >= 0x100'0000'0000);

		// binary calculations
		mtmp = res * (md + 0x0080'0000);
		if (mtmp < 0x8000'0000'0000'0000) {
			uint64_t diff1 = 0x8000'0000'0000'0000 - mtmp;
			uint64_t diff2 = (res + 1) * (md + 0x0080'0000) - 0x8000'0000'0000'0000;
			if (diff2 < diff1) ++res;
		}

		mtmp = etmp; // not to overflow for etmp
		res = (mtmp << 16) + res - 0x0080'0000'0000;
		// binary +1 to mantissa
		return res;
	}

	static uint32_t fma(uint32_t a, uint32_t b, uint32_t c, float& example) { // a*b + c;
		float dummy; //
		bool coutflag = true; //
		if (coutflag) cout << endl << hex << a << " " << b << " " << c << endl; //
		if (coutflag) cout << float(FP32(a)) << " " << float(FP32(b)) << " " << float(FP32(c)) << endl;//
		example = std::fmaf(float(FP32(a)), float(FP32(b)), float(FP32(c))); //
		uint32_t res = (a ^ b) & 0x80000000;
		int32_t ea = a & 0x7F800000;
		int32_t eb = b & 0x7F800000;
		int32_t ec = c & 0x7F800000;
		int32_t eres;
		uint64_t ma = a & 0x007FFFFF;
		uint64_t mb = b & 0x007FFFFF;
		int64_t mc = c & 0x007FFFFF;
		int64_t mres;

		// NAN AND INF
		if ((ea == 0x7F800000) || (eb == 0x7F800000)) { // nan and inf checking for mul
			if (ma != 0 && eb == 0x7F800000) // if left is nan
				return res + (a & 0x7FFFFFFF);
			else if (eb == 0x7F800000) // if right is inf or nan
				return res + (b & 0x7FFFFFFF);
			return res + (a & 0x7FFFFFFF);
		}
		if (ec == 0x7F800000) { // inf - inf error(((
			if (mc != 0) return c;
			else return c; // inf - inf
		}

		// MULTIPLYING
		eres = ((ea + eb) >> 23) + (ea == 0) + (eb == 0); // calculating exponent NEW!!!
		eres -= 127;
		mres = (((ma + 0x00800000 * (ea > 0))) * ((mb + 0x00800000 * (eb > 0)))) << 1; // calculating 24 bit extended mantissa
		if (coutflag) cout << "mres: " << mres << ", eres: " << eres << endl; //

		// MUL RESULT TO NORMAL
		while ((mres < 0x8000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^24. Can it be speeded up? // >= 0
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x1'0000'0000'0000) { // mres is greater than 2^24 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1; // just 1? (eres > 0)
		}
		if (mres >= 0x8000'0000'0000 && eres <= 0) { // new
			mres >>= 1;
		}
		if (coutflag) cout << "mres: " << mres << ", eres: " << eres << ", true: " << FP32::mul(a, b, dummy) << endl; //

		// ADDING PREPARATIONS
		mc = (mc + (uint64_t(ec > 0) << 23)) << 24;
		mres *= -(2 * int32_t(res >> 31) - 1);
		mc   *= -(2 * int32_t(c   >> 31) - 1);
		eb = eres;
		ec >>= 23; // ec > 0??
		if (coutflag) cout << "mres: " << mres << ", eb: " << eb << endl; //
		if (coutflag) cout << "mc: " << mc << ", ec: " << ec << endl; //

		// ADDING
		if (eb > ec) { // calulate exponent and mantissa making exponents equal each other
//			ec += (ec == 0);
			eres = eb;
			if (eb - ec >= 64) mres = mres; // 25
			else mres = mres + (mc >> (eb - ec));
		}
		else {
//			ec += (ec == 0);
//			eb += (eb <= 0); // new
			eres = ec;
			if (ec - eb >= 64) mres = mc;
			else mres = mc + (mres >> (ec - eb));
		}
		res = 0x80000000 * (mres < 0); // abs of mres, creating res sign
		mres *= -(2 * (mres < 0) - 1);
		if (coutflag) cout << "mres: " << mres << ", eb: " << eb << endl; //

		// ADD RES NORMALIZATION
		if (mres >= 0x1'0000'0000'0000) { // if mres is greater than a 2^23 (2^48) // make a non-if code
			mres >>= (eres > 0); // subnormals
			++eres;
		}
		else {
			while (mres < 0x8000'0000'0000 && eres > 0) { // if mres is less than a 2^23 (2^24) // can add mres == 0
				--eres;
				mres <<= (eres > 0); // subnormals
			}
		}
		if (coutflag) cout << "mres: " << mres << ", eb: " << eb << endl; //

		// ROUNDING
		mres = (mres >> 24) + ((mres & 0xFF'FFFF) > 0x80'0000) + ((mres & 0x1FF'FFFF) == 0x180'0000); // last 24 bit stored 
		if (mres >= 0x01000000) { // mres is greater than 2^23 // add this into >= 2^48
			mres >>= (eres > 0); // subnormals
			++eres;
		}

		if (eres >= 0xFF) { // in inf
			return res + 0x7F800000;
		}
		res += eres << 23; // calcuating result
		if (eres > 0)
			return res + mres - 0x00800000;
		else
			return res + mres;

//		return res;


		//if (eres > 0) { // if normal
		//	mres = (mres >> 24) + ((mres & 0xFF'FFFF) > 0x80'0000) + ((mres & 0x1FF'FFFF) == 0x180'0000); // last 23 bit stored 
		//	if (mres >= 0x01000000) { // mres is greater than 2^23 (2^24) // should I?
		//		eres += 1;
		//		mres >>= 1;
		//	}
		//	if (eres >= 0xFF) { // if infinity
		//		return res | 0x7F800000;
		//	}
		//	else {
		//		return res + uint32_t(mres) - 0x00800000 + (eres << 23); // calculating result
		//	}
		//}
		//else if (eres >= -23) { // if subnormal
		//	uint64_t shift = 1ull << (23 - eres + 1); 
		//
		//	return res + (mres >> (23 - eres + 1)) + 
		//	( (mres & (shift - 1) ) > (shift >> 1) ) + 
		//	( ( (mres & ( (shift << 1) - 1) ) ) == (shift + (shift >> 1) ) );
		//}
		
		/*
		mres = (mres >> 24) + ((mres & 0xFF'FFFF) > 0x80'0000) + ((mres & 0x1FF'FFFF) == 0x180'0000); // last 24 bit stored 
		if (mres >= 0x01000000) { // mres is greater than 2^23
			mres >>= (eres > 0); // subnormals
			++eres;
		}

		if (eres >= 0xFF) { // in inf
			return res + 0x7F800000;
		}
		res += eres << 23; // calcuating result
		if (eres > 0)
			return res + mres - 0x00800000;
		else
			return res + mres;

		return res;
		*/
		/* 
		if (el == 0xFF || er == 0xFF) { // nan and inf checking for add
			if (el == 0xFF) {
				if (ml != 0) return l;
				if (er == 0xFF) {
					if (mr != 0) return r;
					if ((l ^ r) == 0x80000000) return r + 1;
					else return r;
				}
				return l;
			}
			else
				return r;
		}
		*/

		return res;
	}

	static uint32_t usub(uint32_t a) {
		return (a & 0x7FFF'FFFF) + (~a & 0x8000'0000);
	}

	static uint64_t special_mul_32_64_64rouding(uint32_t l, uint64_t r) { // double value as 1 8 39: 16 bit extended mantissa
		uint64_t res = 0;
		uint32_t el = l & 0x7F800000;
		uint32_t er = (r & 0x7F80'0000'0000) >> 16;
		int64_t eres;
		uint64_t ml = l & 0x007F'FFFF;
		uint64_t mr = r & 0x007F'FFFF'FFFF;
		uint64_t mres;


		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res + (l & 0x7FFFFFFF);
			else if (er == 0x7F800000) // if right is inf or nan
				return res + ((r >> 16) & 0x7FFFFFFF);
			return res + (l & 0x7FFFFFFF);
		}

		eres = int((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x0080'0000 * (el > 0))) * ((mr + 0x0080'0000'0000 * (er > 0))); // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

		while ((mres < 0x4000'0000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x8000'0000'0000'0000) { // mres is greater than 2^23 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1;
		}


		if (eres > 0) { // if normal
			mres = (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored, rounding
			if (mres >= 0x0100'0000'0000) { // mres is greater than 2^23 (2^24) // should I?
				eres += 1;
				mres >>= 1;
			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F80'0000'0000;
			}
			else {
				return res + mres - 0x0080'0000'0000 + (eres << (23 + 16)); // calculating result
			}
		}
		else if (eres >= - 23 - 16) { // if subnormal
//			mres >>= (-eres + 1);
			uint64_t shift = 1ull << (23 - eres + 1);
//			return res + (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored

//			cout << dec << endl << eres << hex << " " << shift - 1 << " " << (shift >> 1) << " " << (shift << 1) - 1 << " " << shift + (shift >> 1) << endl;
			return res + (mres >> (23 - eres + 1)) +
				((mres & (shift - 1)) > (shift >> 1)) +
				(((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
		}

		return res;
	}

	static uint32_t special_mul_32_64_32rouding(uint32_t l, uint64_t r) { // double value as 1 8 39: 16 bit extended mantissa
		uint32_t res = 0;
		uint32_t el = l & 0x7F800000;
		uint32_t er = (r & 0x7F80'0000'0000) >> 16;
		int64_t eres;
		uint64_t ml = l & 0x007F'FFFF;
		uint64_t mr = r & 0x007F'FFFF'FFFF;
		uint64_t mres;


		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res + (l & 0x7FFFFFFF);
			else if (er == 0x7F800000) // if right is inf or nan
				return res + ((r >> 16) & 0x7FFFFFFF);
			return res + (l & 0x7FFFFFFF);
		}

		eres = int((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x0080'0000 * (el > 0))) * ((mr + 0x0080'0000'0000 * (er > 0))); // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

		while ((mres < 0x4000'0000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x8000'0000'0000'0000) { // mres is greater than 2^23 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1;
		}


		if (eres > 0) { // if normal
			mres = (mres >> 39) + ((mres & 0x7F'FFFF'FFFF) > 0x40'0000'0000) + ((mres & 0xFF'FFFF'FFFF) == 0xC0'0000'0000); // last 23 bit stored, rounding
			if (mres >= 0x0100'0000) { // mres is greater than 2^23 (2^24) // should I?
				eres += 1;
				mres >>= 1;
			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F800000;
			}
			else {
				return res + uint32_t(mres) - 0x00800000 + (eres << 23); // calculating result
			}
		}
		else if (eres >= - 23) { // if subnormal
			//			mres >>= (-eres + 1);
			uint64_t shift = 1ull << (23 - eres + 1 + 16); //+
			//			return res + (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored

			//			cout << dec << endl << eres << hex << " " << shift - 1 << " " << (shift >> 1) << " " << (shift << 1) - 1 << " " << shift + (shift >> 1) << endl;
			return res + (mres >> (23 - eres + 1 + 16)) +
				((mres & (shift - 1)) > (shift >> 1)) +
				(((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
		}

		return res;
	}

	static uint64_t special_mul_32_32_64norouding(uint32_t l, uint32_t r) { // double value as 1 8 46: 16 bit extended mantissa
		uint64_t res = (l ^ r) & 0x7FFF'FFFF;
		uint32_t el = l & 0x7F80'0000;
		uint32_t er = r & 0x7F80'0000;
		int64_t eres;
		uint64_t ml = l & 0x007F'FFFF;
		uint64_t mr = r & 0x007F'FFFF;
		uint64_t mres;


		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res + (l & 0x7FFFFFFF);
			else if (er == 0x7F800000) // if right is inf or nan
				return res + (r & 0x7FFFFFFF);
			return res + (l & 0x7FFFFFFF);
		}

		eres = int((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x0080'0000 * (el > 0))) * ((mr + 0x0080'0000 * (er > 0))); // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

		while ((mres < 0x4000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x8000'0000'0000) { // mres is greater than 2^23 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1; // is it rounding?
		}


		if (eres > 0) { // if normal
//			mres = (mres >> 39) + ((mres & 0x7F'FFFF'FFFF) > 0x40'0000'0000) + ((mres & 0xFF'FFFF'FFFF) == 0xC0'0000'0000); // last 23 bit stored, rounding
//			if (mres >= 0x0100'0000) { // mres is greater than 2^23 (2^24) // should I?
//				eres += 1;
//				mres >>= 1;
//			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F800000;
			}
			else {
				return res + uint64_t(mres) - 0x4000'0000'0000 + (eres << 46); // calculating result
			}
		}
		else if (eres >= -23) { // if subnormal
			//			mres >>= (-eres + 1);
			uint64_t shift = 1ull << (- eres + 1); //
			//			return res + (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored

			//			cout << dec << endl << eres << hex << " " << shift - 1 << " " << (shift >> 1) << " " << (shift << 1) - 1 << " " << shift + (shift >> 1) << endl;
			return res + (mres >> (- eres + 1)) +
				((mres & (shift - 1)) > (shift >> 1)) +
				(((mres & ((shift << 1) - 1))) == (shift + (shift >> 1))); // this is actually rouding
		}

		return res;
	}

	static uint64_t special_mul_32_64_64norouding(uint32_t l, uint64_t r) { // double value as 1 8 39: 16 bit extended mantissa
		uint64_t res = 0;
		uint32_t el = l & 0x7F800000;
		uint32_t er = (r & 0x7F80'0000'0000) >> 16;
		int64_t eres;
		uint64_t ml = l & 0x007F'FFFF;
		uint64_t mr = r & 0x007F'FFFF'FFFF;
		uint64_t mres;


		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res + (l & 0x7FFFFFFF);
			else if (er == 0x7F800000) // if right is inf or nan
				return res + ((r >> 16) & 0x7FFFFFFF);
			return res + (l & 0x7FFFFFFF);
		}

		eres = int((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x0080'0000 * (el > 0))) * ((mr + 0x0080'0000'0000 * (er > 0))); // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

		while ((mres < 0x4000'0000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x8000'0000'0000'0000) { // mres is greater than 2^23 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1;
		}


		if (eres > 0) { // if normaд
			mres = (mres >> 7) + ((mres & 0x7F) > 0x40) + ((mres & 0xFF) == 0xC0); // last 23 bit stored, rounding
			if (mres >= 0x0100'0000'0000'0000) { // mres is greater than 2^23 (2^24) // should I?
				eres += 1;
				mres >>= 1;
			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F80'0000'0000'0000;
			}
			else {
				return res + mres - 0x0080'0000'0000'0000 + (eres << (23 + 32)); // calculating result
			}
		}
		else if (eres >= - 23 - 32) { // if subnormal
			//			mres >>= (-eres + 1);
			uint64_t shift = 1ull << (7 - eres + 1);
			//			return res + (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored

			//			cout << dec << endl << eres << hex << " " << shift - 1 << " " << (shift >> 1) << " " << (shift << 1) - 1 << " " << shift + (shift >> 1) << endl;
			return res + (mres >> (7 - eres + 1)) +
				((mres & (shift - 1)) > (shift >> 1)) +
				(((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
		}

		return res;
	}

	static uint64_t special_mul_32_64_64noExp(uint32_t l, uint64_t r) { // double value as 1 8 39: 16 bit extended mantissa
		uint64_t res = 0;
		uint32_t el = l & 0x7F800000;
		uint32_t er = (r & 0x7F80'0000'0000) >> 16;
		int64_t eres;
		uint64_t ml = l & 0x007F'FFFF;
		uint64_t mr = r & 0x007F'FFFF'FFFF;
		uint64_t mres;


		if ((el == 0x7F800000) || (er == 0x7F800000)) { // nan and inf
			if (ml != 0 && el == 0x7F800000) // if left is nan
				return res + (l & 0x7FFFFFFF);
			else if (er == 0x7F800000) // if right is inf or nan
				return res + ((r >> 16) & 0x7FFFFFFF);
			return res + (l & 0x7FFFFFFF);
		}

		eres = int((el + er) >> 23) - 127 + (el == 0) + (er == 0); // calculating exponent
		mres = ((ml + 0x0080'0000 * (el > 0))) * ((mr + 0x0080'0000'0000 * (er > 0))); // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

		while ((mres < 0x4000'0000'0000'0000) && (eres >= 0)) { // mres has no leading bit 2^23. Can it be speeded up?
			eres -= 1;
			mres <<= 1;
		}
		while (mres >= 0x8000'0000'0000'0000) { // mres is greater than 2^23 (2^?). Can it be speeded up?
			eres += 1;
			mres >>= 1;
		}


		if (eres > 0) { // if normaд
//			mres = (mres >> 7) + ((mres & 0x7F) > 0x40) + ((mres & 0xFF) == 0xC0); // last 23 bit stored, rounding
			if (mres >= 0x8000'0000'0000'0000) { // mres is greater than 2^23 (2^24) // should I?
				eres += 1;
				mres >>= 1;
			}
			if (eres >= 0xFF) { // if infinity
				return res | 0x7F80'0000'0000'0000;
			}
			else {
//				return res + mres - 0x0080'0000'0000'0000 + (eres << (23 + 32)); // calculating result
				return mres - 0x4000'0000'0000'0000;
			}
		}
		else if (eres >= -62) { // if subnormal
			//			mres >>= (-eres + 1);
			uint64_t shift = 1ull << (- eres + 1);
			//			return res + (mres >> 23) + ((mres & 0x7F'FFFF) > 0x40'0000) + ((mres & 0xFF'FFFF) == 0xC0'0000); // last 23 bit stored

			//			cout << dec << endl << eres << hex << " " << shift - 1 << " " << (shift >> 1) << " " << (shift << 1) - 1 << " " << shift + (shift >> 1) << endl;
			return res + (mres >> (- eres + 1)) +
				((mres & (shift - 1)) > (shift >> 1)) +
				(((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
		}

		return res;
	}

	static uint32_t div(uint32_t l, uint32_t r, float& example) noexcept {
		uint32_t savedl = l, savedr = r;
		cout << endl << l << " " << r << endl;
		float dummy; //
		example = float(FP32(l)) / float(FP32(r)); //
		uint32_t rCopied = r; //
		uint32_t lCopied = l; //
		uint32_t res = (l ^ r) & 0x80000000;
		uint32_t el = l & 0x7F800000;
		uint32_t er = r & 0x7F800000;
		int32_t eres = el + 127 - er;
		uint64_t ml = l & 0x007FFFFF;
		uint64_t mr = r & 0x007FFFFF;
		uint64_t mres;

		if ((el == 0x7F800000) || (er == 0x7F800000) || (r == 0) || (r == 0x80000000) || (l == 0) || (l == 0x80000000)) { // nan and inf
			if (el == 0x7F800000) { // left is inf or nan
				if (ml != 0) // left is nan
					return res + (l & 0x7FFFFFFF);
				else if (er == 0x7F800000) // right is inf or nan
					return res + 0x7F800001; // inf / inf = nan
				else
					return res + (l & 0x7FFFFFFF); // inf / finite = inf
			}
			if (er == 0x7F800000) { // right is inf or nan
				if (mr != 0) // is nan
					return res + (r & 0x7FFFFFFF);
				else // is inf
					return res; // zero
			}

			if ((r == 0) || (r == 0x80000000)) { // right is zero
				if ((l == 0) || (l == 0x80000000)) // left is zero
					return res | 0x7F800001;
				else
					return res | 0x7F800000;
			}

			if ((l == 0) || (l == 0x80000000)) {
				return res;
			}
		}

//		eres = (el >> 23) + 126 + (mr == 0); // diapason (0.5, 1]
		eres = (el >> 23) + 126; // diapason [0.5, 1)
		if (er == 0) {
			--eres; // er += 1 if right is subnormal
			while (mr < 0x80'0000) { // right to normal
				++eres; // --er
				mr <<= 1;
			}
			mr -= 0x80'0000;
		}
		eres -= (er >> 23);
		if (el == 0) {
			++eres; // el += 1 if left is subnormal
			while (ml < 0x80'0000) { // left to normal
				--eres; // --el
				ml <<= 1;
			}
			ml -= 0x80'0000;
		}

		// calculations 
//		r = mr + ((126 + (mr == 0)) << 23); // diapason (0.5, 1]
		r = mr + (126 << 23); // diapason [0.5, 1)

		//eres = (el >> 23);
		//eres -= (er >> 23) - 126; // if el is subnormal? 
//		cout << hex << endl << "r " << r << endl;
		uint32_t x = FP32::add3(0x4034b4b5, FP32::mul(0xbff0f0f1, r, dummy), dummy); // 48/17 - 32/17 * d

		/*
		x = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(r, x, dummy), dummy), dummy);
		x = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(r, x, dummy), dummy), dummy);
		x = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(r, x, dummy), dummy), dummy);
//		x = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(r, x, dummy), dummy), dummy); // change
		uint32_t a, b;
		a = FP32::sub(r, 0x3f000000, dummy); // a = r - 0.5 // 0x3f000000 = 0.5
		a = FP32::mul3(a, x, dummy); // a = a * x = (r - 0.5) * x
//		cout << hex << endl << a << endl;
		b = FP32::mul3(0x3f000000, x, dummy); // b = 0.5 * x
//		cout << b << endl;
		uint32_t s = FP32::add3(a, b, dummy); // s = a + b
//		cout << endl << s << endl;
		uint32_t t = FP32::sub(FP32::sub(s, b, dummy), a, dummy); // t = (s - b) - a
//		cout << t << endl;
//		x = FP32::mul3(x, FP32::sub(0x40000000, s, dummy), dummy); // x = x * (2 - s)
		uint32_t tmp = FP32::mul3(x, 0x40000000, dummy);
		tmp = FP32::sub(tmp, FP32::mul3(x, s, dummy), dummy);
		tmp = FP32::sub(tmp, FP32::mul3(x, t, dummy), dummy);
		x = tmp;
//		x = FP32::add3(x, FP32::mul3(x, t, dummy), dummy); // x += x * t // -t? */

//		cout << hex << endl << x << " " << r << endl; //
		// but it is possible than there will be an answer on 1st or 2nd iteration. Make a check on each one!
//		cout << hex << x << endl; //
		x = newton_iter2(x, r); // x = x * (2 - r*x) or x = x + x(1 - dx)
//		cout << hex << x << endl; //
		x = newton_iter2(x, r); // x = x * (2 - r*x)
//		cout << hex << x << endl; //
//		cout << FP32::mul3(x, r, dummy) << endl;
		x = special_newton_iter_rouding(x, r); // x = x * (2 - r*x) // make 2 binary search operations? embed this check into newton_iter 
//		cout << hex << x << endl; //
//		cout << hex << FP32::mul3(x, r, dummy) << endl; //
//		if (FP32::mul3(x, r, dummy) > 0x3f800000ul) --x;
//		else if (FP32::mul3(x, r, dummy) < 0x3f800000ul) ++x;
//		x = newton_iter2(x, r); // ? нужна 4ая итерация? fma needed
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);
//		x = newton_iter(x, r);

//		x = newton_iter2(x, r);
//		x = newton_iter2(x, r);
//		x = newton_iter2(x, r);
//		x = newton_iter2(x, r);

// 
//		x = FP32::sub( FP32::mul3( 0x40000000, x, dummy ), FP32::mul3 ( r, FP32::mul3( x, x, dummy ), dummy), dummy );
//		x = FP32::sub( FP32::mul3( 0x40000000, x, dummy ), FP32::mul3 ( r, FP32::mul3( x, x, dummy ), dummy), dummy );
//		x = FP32::sub( FP32::mul3( 0x40000000, x, dummy ), FP32::mul3 ( r, FP32::mul3( x, x, dummy ), dummy), dummy );
//		x = FP32::sub( FP32::mul3( 0x40000000, x, dummy ), FP32::mul3 ( r, FP32::mul3( x, x, dummy ), dummy), dummy );
//		x = FP32::sub( FP32::mul3( 0x40000000, x, dummy ), FP32::mul3 ( r, FP32::mul3( x, x, dummy ), dummy), dummy );

//		x = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(r, x, dummy), dummy), dummy); //
//		x = FP32::mul3(x, FP32::sub(0x40000000, FP32::mul3(r, x, dummy), dummy), dummy); //
//		cout << hex << "r^-1 " << x << endl;

//		++x; //
//		cout << dec << eres << hex << endl;

		if (eres >= 255) 
			return res + 0x7F800000;
		else if (eres > 0) {
			l = (eres << 23) + ml;

			cout << l << endl;
			uint32_t y = FP32::mul(l, x, dummy); // y = l*x = l * 1/r
			cout << res + y << endl;
//			uint32_t d = FP32::sub(l, FP32::mul3(r, y, dummy), dummy); // d = l - r*y = l - r*l*1/r
			r = usub(r);
			float df = std::fmaf(float(FP32(savedr)), float(FP32(y)), float(FP32(savedl)));
			cout << setprecision(10) << df << endl;
			cout << hex << FP32(df).data << " " << y << " " << l << " " << r << endl;
			y = FP32(std::fmaf(df, float(FP32(x)), float(FP32(y)))).data;
			/*
			df = std::fmaf(float(FP32(r)), float(FP32(y)), float(FP32(l)));
			y = FP32(std::fmaf(df, float(FP32(x)), float(FP32(y)))).data;
			*/
			r = usub(r);
//			cout << hex << d << " " << y << endl;
			cout << hex << FP32(df).data << " " << y << endl;
//			y = FP32::add3(y, FP32::mul3(d, x, dummy), dummy); // y = y + d*x = y + d*1/a
			
			cout << y << endl; 
			

//			res += FP32::mul3(l, x, dummy); // + 0x80'0000; // l * r mantissa
			res += y;
		}
		else if (eres >= -23) {
			l = ml + 0x80'0000;
//			l = FP32(float(FP32(l)) / 2.0f).data;
//			x -= (-eres + 1) << 23; // x_exp - shift of exponent during the calculations to [0.5, 1) NO L SHIFT

			cout << x << endl;
			uint32_t y = FP32::mul(l, x, dummy); // y = l*x = l * 1/r
			cout << y << endl;
//			uint32_t d = FP32::sub(l, FP32::mul3(r, y, dummy), dummy); // d = l - r*y = l - r*l*1/r
			r = usub(r);
			float df = std::fmaf(float(FP32(savedr)), float(FP32(y)), float(FP32(savedl)));
			cout << setprecision(20) << df << endl;
			cout << hex << FP32(df).data << " " << float(FP32(y)) << " " << float(FP32(l)) << " " << float(FP32(r)) << endl;
			y = FP32(std::fmaf(df, float(FP32(x)), float(FP32(y)))).data;
			/*
			df = std::fmaf(float(FP32(r)), float(FP32(y)), float(FP32(l)));
			y = FP32(std::fmaf(df, float(FP32(x)), float(FP32(y)))).data;
			*/
			r = usub(r);
			cout << hex << FP32(df).data << " " << y << endl;
//			y = FP32::add3(y, FP32::mul3(d, x, dummy), dummy); // y = y + d*x = y + d*1/a

			cout << y << endl; 
			

//			res += FP32::mul3(l, x, dummy);
			res += y;
		}
//		else if (eres >= -24)
//			res += 1; // ?????
//			return res;
		else 
			return res;

		return res;

//		cout << FP32::mul3(rCopied, res, dummy) << endl; //
//		cout << endl; //
		// check here what num is more suitable
		/*
		uint64_t mres1, mres2, mres3, diff1, diff2, diff3;
		mres = lCopied; // ????????????????
		mres <<= 23; // exanding l to 2^46   !!!23!!!
//		cout << lCopied << " " << mres << endl; //
		mres1 = special_mul(rCopied, res); // counting r*res without rouding
//		if ((res & 0x7fffffff) == 0x0) mres1 += 0x400000000000; // ??
//		cout << res << endl;
		if ((res & 0x7fffffff) < 0x7f7fffff) {
			mres2 = special_mul(rCopied, res + 1);
		}
		else mres2 = mres1;
		if ((res & 0x7fffffff) > 0x0) {
			mres3 = special_mul(rCopied, res - 1);
//			if ((res & 0x7fffffff) == 0x1) mres3 += 0x400000000000; // ??
		}
		else mres3 = mres1;
//		cout << mres1 << " " << mres2 << " " << mres3 << endl;
//		cout << res << " " << res + 1 << " " << res - 1 << endl;
//		if ((mres1 & 0x3F'FFFF'FFFF'FFFF) == 0 || (mres1 & 0x3F'FFFF'FFFF'FFFF) == 0x3F'C000'0000'0000) diff1 = 0xFFFF'FFFF'FFFF'FFFF;
		if (mres > mres1) diff1 = mres - mres1; // else
		else diff1 = mres1 - mres;
//		if ((mres2 & 0x3F'FFFF'FFFF'FFFF) == 0 || (mres2 & 0x3F'FFFF'FFFF'FFFF) == 0x3F'C000'0000'0000) diff2 = 0xFFFF'FFFF'FFFF'FFFF;
		if (mres > mres2) diff2 = mres - mres2; // else
		else diff2 = mres2 - mres;
//		if ((mres3 & 0x3F'FFFF'FFFF'FFFF) == 0 || (mres3 & 0x3F'FFFF'FFFF'FFFF) == 0x3F'C000'0000'0000) diff3 = 0xFFFF'FFFF'FFFF'FFFF;
		if (mres > mres3) diff3 = mres - mres3; // else
		else diff3 = mres3 - mres;

//		cout << diff1 << " " << diff2 << " " << diff3 << endl;

		if (diff1 == diff2) {
			if (res % 2 == 0) return res;
			else return res + 1;
		}
		else if (diff1 < diff2) {
			if (diff3 < diff1) return res - 1;
			else if (diff3 == diff1) {
				if (res % 2 == 0) return res;
				else return res - 1;
			}
			else return res;
		}
		else if (diff2 < diff1) {
			if (diff3 < diff2) return res - 1;
			else return res + 1;
		}
		*/


//		if (FP32::mul3(rCopied, res, dummy) < l) ++res; // bad

		return res;
	}

	static uint32_t zagryadskov_iter(uint64_t l, uint32_t r, uint64_t y, int iterCount = 5) { // l\r = y; r*y = l
		uint64_t tmp;
//		l &= 0x7F80'0000;
		l <<= 32 + 7; // 16
		cout << endl << "input: " << hex << " l: " << l << " r: " << r << " y: " << y << endl;

		uint64_t currentBit = 0x8000'0000'0000'0000; // idea is to substract a one from bit. E.g, res < 0x123456 => r += curbit, res > 0x123456 => r -= curbit
		currentBit >>= 56; // 56 bits are correct. 16 bits are zero, 1 bit is sign (0), 8 bits are exp, 23 bits are mantissa (48 by this moment), 8 bit idk
		for (int _ = 0; _ < iterCount; ++_) {
//			tmp = r * y; // r, y - 64 bit
			tmp = special_mul_32_64_64noExp(r, y); // extend tmp and l up to the 64 bit long 
			cout << hex << "r*y: " << tmp << " y: " << y << endl;
			if (tmp < l) y += currentBit;
			else if (tmp > l) y -= currentBit;
			currentBit >>= 1;
		}

		y = (y >> 16) + ((y & 0xFFFF) > 0x8000) + ((y & 0x1'FFFF) == 0x1'8000); // last 16 bit stored 
//		y = (y >> 16) + ((y & 0xFFFF) >= 0x8000);
		if ((y & 0x7F80'0000) == 0x7F80'0000) return 0x7F80'0000;
		return y;
	}

	static uint32_t closedNum(uint64_t l, uint32_t r, uint32_t y) {
		uint32_t lefty = y, righty = y; // trying to find closest value
		if ((y & 0x7FFF'FFFF) != 0x7FFF'FFFF) righty += 1; // if y not an inf
		if ((y & 0x7FFF'FFFF) != 0x0000'0000) lefty  -= 1; // if y not a zero

		uint64_t leftl, midl, rightl; // calculating approximations
		leftl  = FP32::special_mul_32_32_64norouding(r, lefty );
		midl   = FP32::special_mul_32_32_64norouding(r, y     );
		rightl = FP32::special_mul_32_32_64norouding(r, righty);
		l <<= 23;

		cout << hex << endl << leftl << " " << midl << " " << rightl << " " << l << endl;

		uint64_t diffl, diffm, diffr; // calculating differences
		if (leftl < l) diffl = l - leftl;
		else diffl = leftl - l;
		if (midl < l) diffm = l - midl;
		else diffm = midl - l;
		if (rightl > l) diffr = rightl - l;
		else diffr = l = rightl;

		cout << hex << endl << diffl << " " << diffm << " " << diffr << endl;

		if (diffl < diffm) {
			if (diffl < diffr) return lefty;
			if (diffl == diffr) return y;
			if (diffl > diffr) return righty;
		}
		if (diffl == diffm) {
			if (diffr < diffl) return righty;
			if (diffr == diffl) return y;
			if (diffr > diffl) {
				if ((y & 0x1) == 0x0) return y;
				if ((y & 0x1) == 0x1) return lefty;
			}
		}
		if (diffl > diffm) {
			if (diffr < diffm) return righty;
			if (diffr == diffm) {
				if ((y & 0x1) == 0x0) return y;
				if ((y & 0x1) == 0x1) return righty;
			}
			if (diffr > diffm) return y;
		}

		return y;
	}

	static uint32_t div2(uint32_t l, uint32_t r, float& example) noexcept {
//		cout << endl << l << " " << r << endl; //
		float dummy; //
		example = float(FP32(l)) / float(FP32(r)); //
		uint32_t rCopied = r & 0x7FFF'FFFF; //
		uint32_t lCopied = l & 0x7FFF'FFFF; //
		uint32_t res = (l ^ r) & 0x8000'0000;
		uint32_t el = l & 0x7F80'0000;
		uint32_t er = r & 0x7F80'0000;
		int32_t eres = el + 127 - er;
		uint64_t ml = l & 0x007F'FFFF;
		uint64_t mr = r & 0x007F'FFFF;
		uint64_t mres;


		if ((el == 0x7F800000) || (er == 0x7F800000) || (r == 0) || (r == 0x80000000) || (l == 0) || (l == 0x80000000)) { // nan and inf
			if (el == 0x7F800000) { // left is inf or nan
				if (ml != 0) // left is nan
					return res + (l & 0x7FFFFFFF);
				else if (er == 0x7F800000) // right is inf or nan
					return res + 0x7F800001; // inf / inf = nan
				else
					return res + (l & 0x7FFFFFFF); // inf / finite = inf
			}
			if (er == 0x7F800000) { // right is inf or nan
				if (mr != 0) // is nan
					return res + (r & 0x7FFFFFFF);
				else // is inf
					return res; // zero
			}

			if ((r == 0) || (r == 0x80000000)) { // right is zero
				if ((l == 0) || (l == 0x80000000)) // left is zero
					return res | 0x7F800001;
				else
					return res | 0x7F800000;
			}

			if ((l == 0) || (l == 0x80000000)) {
				return res;
			}
		}


		eres = (el >> 23) + 126; // diapason [0.5, 1)
		if (er == 0) {
			--eres; // er += 1 if right is subnormal
			while (mr < 0x80'0000) { // right to normal
				++eres; // --er
				mr <<= 1;
			}
			mr -= 0x80'0000;
		}
		eres -= (er >> 23);
		if (el == 0) {
			++eres; // el += 1 if left is subnormal
			while (ml < 0x80'0000) { // left to normal
				--eres; // --el
				ml <<= 1;
			}
			ml -= 0x80'0000;
		}
		r = mr + (126 << 23); // diapason [0.5, 1)

		uint32_t x = FP32::add3(0x4034b4b5, FP32::mul(0xbff0f0f1, r, dummy), dummy); // 48/17 - 32/17 * d
		x = newton_iter2(x, r); // x = x * (2 - r*x) or x = x + x(1 - dx)
		x = newton_iter2(x, r); // x = x * (2 - r*x)
//		uint64_t bigx = special_newton_iter_no_rounding(x, r); // x = x * (2 - r*x) // make 2 binary search operations? embed this check into newton_iter 
		x = special_newton_iter_rouding(x, r);

		if (eres >= 255) // l / r = y => l = r*y
			return res + 0x7F800000;
		else if (eres > 0) {
//			cout << endl << "normal branch" << endl;
//			l = (eres << 23) + ml; // shift l here to diapason [0.5, 1). USUAL
			l = (126 << 23) + ml; // FMA
			uint32_t y = FP32::mul(l, x, dummy); // FMA
//			uint32_t y = FP32::special_mul_32_64_32rouding(l, bigx);

//			uint64_t y = FP32::special_mul_32_64_64rouding(l, bigx);
//			y = zagryadskov_iter(l, r, y);
// 
			cout << endl;
			if ((y & 0x7FFF'FFFF) != 0x7F80'0000) {
				rCopied = usub(rCopied);
//				cout << FP32(rCopied) << " " << FP32(lCopied) << endl;
				float df = std::fmaf(FP32(rCopied), FP32(y), FP32(lCopied));
//				cout << FP32(y) << " " << df << endl;
				y = FP32(std::fmaf(df, 1.0f/float(FP32(rCopied)), FP32(y))).data;
//				cout << FP32(y) << " " << FP32(x) << endl;
			}
//			y -= (126 - eres) << 23; //!!!
			/*
			eres = -126 + ((y & 0x7F80'0000) >> 23) + eres;
			int eresx = -126 + ((x & 0x7F80'0000) >> 23) + eres;

			if (eres >= 0xFF)
				y = 0xFF << 23;
			else if (eres > 0) {
				y &= 0x807F'FFFF;
				y += eres << 23;
			}
			else if (eres >= -23) {
				y = (y & 0x007F'FFFF) + 0x0080'0000;
				uint64_t shift = 1ull << (-eres + 1);
				y = (y >> (-eres + 1)) +
					((y & (shift - 1)) > (shift >> 1)) +
					(((y & ((shift << 1) - 1))) == (shift + (shift >> 1)));
			}
			else
				y = 0;

			if (eresx >= 0xFF) {
				x = 0xFF << 23;
			}
			else if (eresx > 0) {
				x &= 0x807F'FFFF;
				x += eresx << 23;
			}
			else if (eresx >= -23) {
				x = (x & 0x007F'FFFF) + 0x0080'0000;
//				cout << y << endl;  // 0000 0000 1111 0100 1000 1001 1000 1101
				uint64_t shift = 1ull << (-eres + 1);
//				y >>= -eres + 1;
				x = (x >> (-eres + 1)) + ((x & (shift - 1)) > (shift >> 1)) + (((x & ((shift << 1) - 1))) == (shift + (shift >> 1)));
//				cout << y << endl; // 0000 0000 0111 1010 0100 0100 1100 0110
			}
			else
				x = 0;
			
			y = FP32(std::fmaf(df, FP32(x), FP32(y))).data;
			cout << FP32(y) << " " << FP32(x) << endl;
			} */
			res += uint32_t(y);
		}
		else if (eres >= -30) {
			//			cout << endl << "subnormal branch" << endl;
			int teres = eres;
			l = ml + 0x80'0000;
			l += 126 << 23; // FMA
			x -= (-eres + 1) << 23; // x_exp - shift of exponent during the calculations to [0.5, 1) NO L SHIFT
//			bigx -= (uint64_t(-eres + 1) << (23 + 16)); // 16
			r += (-eres + 1) << 23; // FMA	
			uint32_t y = FP32::mul(l, x, dummy); // y = l*x = l * 1/r // FMA

//			uint64_t y = FP32::special_mul_32_64_64rouding(l, bigx);
//			y = zagryadskov_iter(l, r, y);
			
			
			if ((y & 0x7FFF'FFFF) != 0x7F80'0000) {
				r = usub(r);
//				cout << hex << endl << FP32(r) << " " << FP32(l) << endl; //
				float df = std::fmaf(FP32(r), FP32(y), FP32(l));
//				cout << FP32(y) << " " << df << endl; //
				eres = -126 + ((y & 0x7F80'0000) >> 23);
				int eresx = -126 + ((x & 0x7F80'0000) >> 23);
			//}
//			cout << hex << endl << y << endl; // 0011 1111 0111 0100 1000 1001 1000 1101
			if (eres >= 0xFF) {
				y = 0xFF << 23;
			}
			else if (eres > 0) {
				y &= 0x807F'FFFF;
				y += eres << 23;
			}
			else if (eres >= -23) {
				y = (y & 0x007F'FFFF) + 0x0080'0000;
//				cout << y << endl;  // 0000 0000 1111 0100 1000 1001 1000 1101
				uint64_t shift = 1ull << (-eres + 1);
//				y >>= -eres + 1;
				y = (y >> (-eres + 1)) + ((y & (shift - 1)) > (shift >> 1)) + (((y & ((shift << 1) - 1))) == (shift + (shift >> 1)));
//				cout << y << endl; // 0000 0000 0111 1010 0100 0100 1100 0110
			}
			else
				y = 0;

			if (eresx >= 0xFF) {
				x = 0xFF << 23;
			}
			else if (eresx > 0) {
				x &= 0x807F'FFFF;
				x += eresx << 23;
			}
			else if (eresx >= -23) {
				x = (x & 0x007F'FFFF) + 0x0080'0000;
//				cout << y << endl;  // 0000 0000 1111 0100 1000 1001 1000 1101
				uint64_t shift = 1ull << (-eres + 1);
//				y >>= -eres + 1;
				x = (x >> (-eres + 1)) + ((x & (shift - 1)) > (shift >> 1)) + (((x & ((shift << 1) - 1))) == (shift + (shift >> 1)));
//				cout << y << endl; // 0000 0000 0111 1010 0100 0100 1100 0110
			}
			else
				x = 0;

			y = FP32(std::fmaf(df, FP32(x), FP32(y))).data; 
			}
/*			x += 7 << 23; // FMA	
			cout << "f " << endl;
			rCopied = usub(rCopied);
			cout << hex << endl << FP32(rCopied) << " " << FP32(lCopied) << endl; //
			float df = std::fmaf(FP32(rCopied), FP32(y), FP32(lCopied));
			cout << FP32(y) << " " << df << endl; //
			y = FP32(std::fmaf(df, FP32(x), FP32(y))).data;
			cout << FP32(y) << " " << FP32(x) << endl; //
			eres = -126 + ((y & 0x7F80'0000) >> 23); */

			res += uint32_t(y);
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
	friend std::ostream& operator<<(std::ostream& s, const FP32& f) {
		s << f.example;
		return s;
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
				//l = l.add2(l, r);
	//			l = l.mul_32_32_32(l, r);
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
		int input = 0;
		uint64_t lc, rc;
		uint32_t res;
		float f;
		size_t from = 2;

		vector<uint32_t> vl = { 0x11000, 0x11000, 0x11000, 0x11000, 0x11000, 0x11000, 0x3c900000, 0x1980, 0x1980, 0x1980, 0x11000, 0x11000, 0x11000, 0x3c911000, 0x11000 }; // {0x11000, 0x40011000, 0x811000, 0x811000, 0xaec000, 0xb85000, 0x14ffd180, 0x17ffe800, 0x2e7fd180, 0x317fe800, 0x47ffd180, 0x4affe800, 0x11000, 0x11000};
		vector<uint32_t> vr = { 0x3f00c000, 0x42719000, 0x3c0e6000, 0x42f11000, 0x3c801000, 0x48004000, 0x80009000, 0x27533793, 0x3eac3ed3, 0x30d69c11, 0x231000, 0x383f0000, 0x3c05e000, 0x80009000, 0x3c091000 }; // {0x231000, 0x80231000, 0x511d000, 0x520b000, 0x6fdc000, 0x7ecd000, 0xb47fdc3a, 0x98ffeffe, 0xb47fdc3a, 0x98ffeffe, 0xb47fdc3a, 0x98ffeffe, 0x1166000, 0x48004000 };
		for (size_t i = from; i < 3; ++i) {
			cout << hex << vl[i] << ", " << vr[i] << endl;
			res = FP32::div2(uint32_t(vl[i]), uint32_t(vr[i]), f);
//			res = FP32::fma3(vl[i], vr[i], 0x0, f);
			if (f == f && res != FP32(f).data) {
					cout << hex << endl << vl[i] << ", " << vr[i] << " , that is " << FP32(uint32_t(vl[i])).example << ", " << FP32(uint32_t(vr[i])).example << " ERROR\n";
					cout << f << " expected, " << float(FP32(res)) << " instead\n";
					FP32(f).print();
					cout << bitset<32>(res) << endl;
					input = 1;
					continue;
			}
		}
		if (input) {
			cout << "\nNOT PASSED\n";
			return;
		}

//		cout << "Add\n";

//	for (uint64_t abcd = 0; abcd <= 0xFFFFFFFF; abcd += 69632)
		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 69632) { // 6528 69632 0x11000
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 69632) {
//				res = FP32::add3(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::sub(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::mul3(uint32_t(lc), uint32_t(rc), f);
				res = FP32::div2(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::fma3(uint32_t(lc), uint32_t(rc), uint32_t(abcd), f);
//				if (f == f && res != FP32(f).data && res - 1 != FP32(f).data && res + 1 != FP32(f).data) {
				if (f == f && res != FP32(f).data) {
//					cout << hex << endl << abcd << ", " << lc << ", " << rc << " , that is " << FP32(uint32_t(abcd)).example << ", " << FP32(uint32_t(lc)).example << ", " << FP32(uint32_t(rc)).example << " ERROR\n";
					cout << hex << endl << lc << ", " << rc << " , that is " << FP32(uint32_t(lc)).example << ", " << FP32(uint32_t(rc)).example << " ERROR\n";
					cout << f << " expected, " << float(FP32(res)) << " instead\n";
					FP32(f).print();
					cout << bitset<32>(res) << endl;
					cout << endl << "Continue? 1 - yes, 0 - no\n";
					cin >> input;
//					input = 1;
					if (input) continue;
					return;
				}
			}
		}		
		/*
		cout << "Sub\n";

		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 5417) { // 6528 69632 0x11000
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 5852) {
//				res = FP32::add3(uint32_t(lc), uint32_t(rc), f);
				res = FP32::sub(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::mul3(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::div(uint32_t(lc), uint32_t(rc), f);
//				if (f == f && res != FP32(f).data && res - 1 != FP32(f).data && res + 1 != FP32(f).data) {
				if (f == f && res != FP32(f).data) {
					cout << hex << endl << lc << ", " << rc << " , that is " << FP32(uint32_t(lc)).example << ", " << FP32(uint32_t(rc)).example << " ERROR\n";
					cout << f << " expected, " << float(FP32(res)) << " instead\n";
					FP32(f).print();
					cout << bitset<32>(res) << endl;
					cout << endl << "Continue? 1 - yes, 0 - no\n";
//					cin >> input;
					input = 1;
					if (input) continue;
					return;
				}
			}
		}

		cout << "Mul\n";

		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 5417) { // 6528 69632 0x11000
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 5852) {
//				res = FP32::add3(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::sub(uint32_t(lc), uint32_t(rc), f);
				res = FP32::mul3(uint32_t(lc), uint32_t(rc), f);
//				res = FP32::div(uint32_t(lc), uint32_t(rc), f);
//				if (f == f && res != FP32(f).data && res - 1 != FP32(f).data && res + 1 != FP32(f).data) {
				if (f == f && res != FP32(f).data) {
					cout << hex << endl << lc << ", " << rc << " , that is " << FP32(uint32_t(lc)).example << ", " << FP32(uint32_t(rc)).example << " ERROR\n";
					cout << f << " expected, " << float(FP32(res)) << " instead\n";
					FP32(f).print();
					cout << bitset<32>(res) << endl;
					cout << endl << "Continue? 1 - yes, 0 - no\n";
//					cin >> input;
					input = 1;
					if (input) continue;
					return;
				}
			}
		} */
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
//	cout << hex << 0xa5effull * 0x6cd000ull + (0xa5effull << 23) + (0x6cd000ull << 23) << endl;
//	cout << hex << 0x8a5effull * 0xecd000ull << endl;

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

void foo() {
	int records = 5;
	vector<int> ids;

	std::unordered_map<int, std::string> mp;
	for (int i = 0; i < records; ++i) {
		mp.insert({ids[i], "hash_" + std::to_string(i) + ".txt"});
	}

    // do smth

	int id = 5;
	string str = mp.find(id)->second;

	std::ofstream file;
	string line;
	file.open(str);
	file << line;
}