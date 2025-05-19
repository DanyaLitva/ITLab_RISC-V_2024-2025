#include "FP32.h"

#if !__GNUC__ && !__clang__
int32_t __builtin_clz(uint32_t num) {
	\
		if (num == 0) return 32; \
			int32_t res = 0; \
			while (!(num & 0x8000'0000)) {
				\
					num <<= 1; \
					++res; \
			} \
				return res; \
}

int32_t __builtin_clzll(uint64_t num) {
	\
		if (num == 0) return 64; \
			int32_t res = 0; \
			while (!(num & 0x8000'0000'0000'0000)) {
				\
					num <<= 1; \
					++res; \
			} \
				return res; \
}
#endif

uint32_t FP32::add(uint32_t l, uint32_t r) const noexcept {
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

	if (mres >= 0x1'0000'0000'0000) { // if mres is greater than a 2^23 (2^48) // make a non-if code
		mres >>= (eres > 0); // subnormals
		++eres;
	}
	else if (mres == 0)
		eres = 0;
	else if (__builtin_clzll(mres) > 16 && eres > 0 && (__builtin_clzll(mres) - 16) < eres) { // if mres is less than a 2^23 (2^24) // can add mres == 0
		eres -= (__builtin_clzll(mres) - 16);
		mres <<= (__builtin_clzll(mres) - 16);
	}
	else if (__builtin_clzll(mres) > 16 && eres > 0 && (__builtin_clzll(mres) - 16) >= eres) {
		mres <<= eres - 1;
		eres = 0;
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
	return res + mres - (eres > 0) * 0x00800000;
}

uint32_t FP32::sub(uint32_t l, uint32_t r) const noexcept  {
	r = (r & 0x7FFF'FFFF) + (~r & 0x8000'0000);
	return FP32::add(l, r);
}	

uint32_t FP32::usub(uint32_t a) const noexcept {
	return (a & 0x7FFF'FFFF) + (~a & 0x8000'0000);
}

uint32_t FP32::mul(uint32_t l, uint32_t r) const noexcept  {
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

	if (mres == 0) return res;
	int32_t etmp = eres;
	eres -= std::max(std::min((__builtin_clzll(mres) - 17), etmp), 0);
	mres <<= std::max(std::min((__builtin_clzll(mres) - 17), etmp), 0);
	eres += (17 - __builtin_clzll(mres)) * (mres >= 0x8000'0000'0000);
	mres >>= (17 - __builtin_clzll(mres)) * (mres >= 0x8000'0000'0000);

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
		uint64_t shift = 1ull << (23 - eres + 1);
		return res + (mres >> (23 - eres + 1)) +
			((mres & (shift - 1)) > (shift >> 1)) +
			(((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
	}

	return res;
}

uint32_t FP32::fma(uint32_t a, uint32_t b, uint32_t c) const noexcept {
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
	if ((ea == 0x7F800000) || (eb == 0x7F800000) || (ec == 0x7F800000)) {
		return FP32::add(FP32::mul(a, b), c);
	}
	// MULTIPLYING
	eres = (ea >> 23) + (eb >> 23) + (ea == 0) + (eb == 0); // calculating exponent
	eres -= 127;
	mres = (((ma + 0x00800000 * (ea > 0))) * ((mb + 0x00800000 * (eb > 0)))) << 1; // calculating 24 bit extended mantissa
	
	// MUL RESULT TO NORMAL
	if (mres == 0) { 
		eres = 0;
	}

	// ADDING PREPARATIONS
	mc = (mc + (int64_t(ec > 0) << 23)) << 24;
	mres *= -(2 * int32_t(res >> 31) - 1); // make them signed
	mc *= -(2 * int32_t(c >> 31) - 1);
	eb = eres;
	ec >>= 23;
	bool ecDenormal = (ec == 0);
	ec += ecDenormal;

	// ADDING
	uint64_t mresLShift = 0;
	uint8_t LShift = 0;
	bool mresLShiftSign = 0;
	if (eb >= ec) { // calulate exponent and mantissa making exponents equal each other
		eres = eb;
		if (eb - ec >= 64) mres = mres; // 25
		else { 
			LShift = eb - ec;
			mresLShift = mc & ((1ull << (eb - ec)) - 1);
			mresLShiftSign = (mc < 0); // mc
			mres = mres + (mc >> (eb - ec));
		}
	}
	else {
		eres = ec;
		if (ec - eb >= 64) mres = mc;
		else {
			LShift = ec - eb;
			mresLShift = mres & ((1ull << (ec - eb)) - 1); // mresLShift is signed, you should calculate this
			mresLShiftSign = (mres < 0); // mres
			mres = mc + (mres >> (ec - eb)); // ERROR! NO ROUNDING PLEASE
		}
	}
	res = 0x80000000 * (mres < 0); // creating res sign
	bool resSign = (mres < 0);
	mres *= -(2 * (mres < 0) - 1); // abs of mres

	// ADD RES NORMALIZATION
	if (mres >= 0x1'0000'0000'0000) { // if mres is greater than a 2^23 (2^48) // make a non-if code 0x1'0000'0000'0000
		if ((mres >> 48) == 2) {
			mresLShift += ((mres & 0x3ull) << LShift);
			mres >>= 2;
			eres += 2;
			LShift += 2;
		}
		else {
			mresLShift += ((mres & 0x1ull) << LShift);
			mres >>= 1;
			eres += 1;
			LShift += 1;
		}
	}
	else if (mres == 0)
		eres = 0;
	else if (__builtin_clzll(mres) > 16 && eres > 0 && (__builtin_clzll(mres) - 16) < eres) { // if mres is less than a 2^23 (2^24) // can add mres == 0
		eres -= (__builtin_clzll(mres) - 16);
		mres <<= (__builtin_clzll(mres) - 16);
	}
	else if (__builtin_clzll(mres) > 16 && eres > 0 && (__builtin_clzll(mres) - 16) >= eres) {
		mres <<= eres - 1;
		eres = 0;
	}

	if (eres >= 0xFF) { // in inf
		return res + 0x7F80'0000;
	}
	else if (eres > 0) { // use here all your saved stuff (mresLShift)
		uint32_t r1 = ((mres & 0xFF'FFFF) > 0x80'0000) + ((mres & 0xFF'FFFF) == 0x80'0000) * (mresLShift != 0x0) * (resSign == mresLShiftSign);
		uint32_t r2 = ((mres & 0x1FF'FFFF) == 0x180'0000) * (mresLShift == 0x0);
		mres = (mres >> 24) + r1 + r2;
		return res + mres - 0x0080'0000 + (eres << 23);
	}
	else if (eres >= -23) {
		uint64_t shift = 1ull << (24 - eres);
		uint32_t r1 = ((mres & (shift - 1)) > (shift >> 1)) + ((mres & (shift - 1)) == (shift >> 1)) * (mresLShift != 0x0) * (resSign == mresLShiftSign);
		uint32_t r2 = (((mres & ((shift << 1) - 1))) == (shift + (shift >> 1))) * (mresLShift == 0x0);
		return res + (mres >> (24 - eres)) + r1 + r2;
	}
	else {
		return res;
	}
	return res;
}

uint32_t FP32::div(uint32_t l, uint32_t r) const noexcept {
	uint32_t res = (l ^ r) & 0x8000'0000;
	uint32_t el = l & 0x7F80'0000;
	uint32_t er = r & 0x7F80'0000;
	int32_t eres = 0;
	uint64_t ml = l & 0x007F'FFFF;
	uint64_t mr = r & 0x007F'FFFF;
	bool coutflag = false;


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


	eres = (el >> 23) + 126; // diapason [1, 2)
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
	r = mr + (127ul << 23); // diapason [1, 2)
	l = ml + (127ul << 23);

	uint32_t x = fma(0xbef0f0f1ul, r, 0x3fb4b4b5ul); // 24/17 - 8/17 * d.
	r = usub(r);
	uint32_t e;

	e = fma(r, x, 0x3F800000ul); // 1 iteration
	x = fma(x, e, x);
	e = fma(r, x, 0x3F800000ul); // 2 iteration
	x = fma(x, e, x);
	e = fma(r, x, 0x3F800000ul); // 3 iteration AND RN(1/a) (EXCEPT FOR UNTEGRAL SIGNIFICAND 11...11)
	x = fma(x, e, x);

	uint32_t y = mul(l, x); // Probably, 1 ulp error. Testing required. y = l / r
	e = fma(r, y, l);
	y = fma(e, x, y);
	e = fma(r, y, l);
	uint32_t tmpy = fma(e, x, y);
	bool dy = y != tmpy; // for possible double-roundings
	bool de = (e >> 31) == 0x1;
	y = tmpy;

	r = usub(r);
	eres += (y >> 23) - 126; // overflow if y >= 1

	if (eres >= 255)
		res += 0x7F800000;
	else if (eres > 0) {
		res += (eres << 23) + (y & 0x007F'FFFF);
	}
	else if (eres >= -23) { // PAGE 137, r shows the rounding
		y = (y & 0x007F'FFFF) + 0x0080'0000;
		uint32_t shift = 1ull << (-eres + 1);
		if (e == 0x0 || ((y & (shift - 1)) != (shift >> 1))) { // round-to-nearest ties-to-even if theorem case or midpoint
			res += (y >> (-eres + 1)) +
				((y & (shift - 1)) > (shift >> 1)) +
				(((y & ((shift << 1) - 1))) == (shift + (shift >> 1)));
		}
		
		else if (dy && de || !dy && !de) { // underestimates
			res += (y >> (-eres + 1)) + 1;
		}
		else if (dy && !de || !dy && de) { // overestimates
			res += (y >> (-eres + 1)) + 0;
		}
		
	}

	return res;
}
FP32::FP32(uint32_t t): data(t) {}
FP32::FP32(float f) noexcept {
	data = reinterpret_cast<uint32_t*>(&f)[0];
}
FP32::FP32(double f) noexcept {
	float ff = f;
	data = reinterpret_cast<uint32_t*>(&ff)[0];
}
FP32::FP32(const FP32& fp): data(fp.data) {}
FP32& FP32::operator=(const FP32& fp) noexcept {
	data = fp.data;
	return *this;
}
FP32& FP32::operator=(uint32_t t) noexcept {
	data = t;
	return *this;
}
FP32& FP32::operator=(float f) noexcept  {
	data = FP32(f).data;
	return *this;
}
FP32& FP32::operator=(double f) noexcept {
	data = FP32(f).data;
	return *this;
}
FP32::operator float() const noexcept {
	uint32_t tmp = data;
	return reinterpret_cast<float*>(&tmp)[0];
}
FP32::operator double() const noexcept {
	uint32_t tmp = data;
	return double(reinterpret_cast<float*>(&tmp)[0]);
}

FP32 FP32::operator- () noexcept {
	uint32_t tmp = (data & 0x7FFF'FFFF) + (~data & 0x8000'0000);
	return FP32(tmp);
}
FP32& FP32::operator+= (const FP32& fp) noexcept {
	data = add(data, fp.data);
	return *this;
}
FP32& FP32::operator-= (const FP32& fp) noexcept {
	data = sub(data, fp.data);
	return *this;
}
FP32& FP32::operator*= (const FP32& fp) noexcept {
	data = mul(data, fp.data);
	return *this;
}
FP32& FP32::operator/= (const FP32& fp) noexcept {
	data = div(data, fp.data);
	return *this;
}
FP32 FP32::operator+ (const FP32& fp) const noexcept {
	FP32 tmp = *this;
	tmp += fp;
	return tmp;
}
FP32 FP32::operator- (const FP32& fp) const noexcept {
	FP32 tmp = *this;
	tmp -= fp;
	return tmp;
}
FP32 FP32::operator* (const FP32& fp) const noexcept {
	FP32 tmp = *this;
	tmp *= fp;
	return tmp;
}
FP32 FP32::operator/ (const FP32& fp) const noexcept {
	FP32 tmp = *this;
	tmp /= fp;
	return tmp;
}
FP32 FP32::fma(const FP32& a, const FP32& b) const noexcept {
	FP32 tmp = *this;
	tmp = fma(data, a.data, b.data);
	return tmp;
}

std::ostream& operator << (std::ostream& stream, const FP32& fp) {
	float tmp = fp;
	stream << tmp;
	return stream;
}
std::istream& operator >> (std::istream& stream, FP32& fp) {
	float tmp;
	stream >> tmp;
	fp = tmp;
	return stream;
}

bool FP32::operator == (const FP32& fp) const noexcept {
	if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
	return data == fp.data;
};
bool FP32::operator < (const FP32& fp) const noexcept {
	if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
	FP32 tmp = *this - fp;
	return tmp.data > 0x80000000ul;
}
bool FP32::operator > (const FP32& fp) const noexcept {
	if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
	FP32 tmp = *this - fp;
	return (tmp.data < 0x80000000ul) && (tmp.data != 0ul);
}
bool FP32::operator >= (const FP32& fp) const noexcept {
	if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
	FP32 tmp = *this - fp;
	return (tmp.data <= 0x80000000ul);
}
bool FP32::operator <= (const FP32& fp) const noexcept {
	if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
	FP32 tmp = *this - fp;
	return (tmp.data >= 0x80000000ul) && (tmp.data == 0ul);
}
bool FP32::operator != (const FP32& fp) const noexcept {
	if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return true;
	return data != fp.data;
}
