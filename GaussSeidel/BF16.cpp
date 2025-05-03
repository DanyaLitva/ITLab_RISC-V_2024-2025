#include "BF16.h"

uint16_t BF16::add(uint16_t l, uint16_t r) const noexcept {
	uint16_t res;
	uint16_t el = (l & 0x7F80) >> 7;
	uint16_t er = (r & 0x7F80) >> 7;
	uint16_t eres;
	int32_t ml = ((uint32_t(l) & 0x007F) + ((uint32_t(el > 0)) << 7)) << 16; // 16 extra bit
	int32_t mr = ((uint32_t(r) & 0x007F) + ((uint32_t(er > 0)) << 7)) << 16;
	int32_t mres;
	ml *= -(2 * int32_t(l >> 15) - 1);
	mr *= -(2 * int32_t(r >> 15) - 1);

	if (el == 0xFF || er == 0xFF) { // nan and inf checking
		if (el == 0xFF) {
			if (ml != 0) return l;
			if (er == 0xFF) {
				if (mr != 0) return r;
				if ((l ^ r) == 0x8000) return r + 1;
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
		if (el - er >= 32) mres = ml; // 25
		else mres = ml + (mr >> (el - er));
	}
	else {
		er += (er == 0);
		el += (el == 0);
		eres = er;
		if (er - el >= 32) mres = mr;
		else mres = mr + (ml >> (er - el));
	}

	res = 0x8000 * (mres < 0);
	mres *= -(2 * (mres < 0) - 1);

	if (mres >= 0x100'0000) { // if mres is greater than a 2^23 (2^48) // make a non-if code
		mres >>= (eres > 0); // subnormals
		++eres;
	}
	else if (mres == 0)
		eres = 0;
	else if (__builtin_clz(uint32_t(mres)) > 8 && eres > 0 && (__builtin_clz(uint32_t(mres)) - 8) < eres) { // if mres is less than a 2^23 (2^24) // can add mres == 0
		eres -= (__builtin_clz(uint32_t(mres)) - 8);
		mres <<= (__builtin_clz(uint32_t(mres)) - 8);
	}
	else if (__builtin_clz(uint32_t(mres)) > 8 && eres > 0 && (__builtin_clz(uint32_t(mres)) - 8) >= eres) {
		mres <<= eres - 1;
		eres = 0;
	}

	mres = (mres >> 16) + ((mres & 0xFFFF) > 0x8000) + ((mres & 0x1'FFFF) == 0x1'8000); // last 24 bit stored 
	if (mres >= 0x0100) { // mres is greater than 2^23 // add this into >= 2^48
		mres >>= (eres > 0); // subnormals
		++eres;
	}

	if (eres >= 0xFF) { // in inf
		return res + 0x7F80;
	}
	res += eres << 7; // calcuating result
	return res + mres - (eres > 0) * 0x0080;
	
}

uint16_t BF16::usub(uint16_t r) const noexcept {
	return (r & 0x7FFFu) + (~r & 0x8000u);
}

uint16_t BF16::sub(uint16_t l, uint16_t r) const noexcept {
	r = usub(r);
	return add(l, r);
}

uint16_t BF16::mul(uint16_t l, uint16_t r) const noexcept {

	uint16_t res = (l ^ r) & 0x8000;
	uint16_t el = l & 0x7F80;
	uint16_t er = r & 0x7F80;
	int32_t eres;
	uint32_t ml = l & 0x007F;
	uint32_t mr = r & 0x007F;
	uint32_t mres;

	if ((el == 0x7F80) || (er == 0x7F80)) { // nan and inf
		if (ml != 0 && el == 0x7F80) // if left is nan
			return res + (l & 0x7FFF);
		else if (er == 0x7F80) // if right is inf or nan
			return res + (r & 0x7FFF);
		return res + (l & 0x7FFF);
	}

	eres = ((el + er) >> 7) - 127 + (el == 0) + (er == 0); // calculating exponent
	mres = ((ml + 0x0080 * (el > 0))) * ((mr + 0x0080 * (er > 0))) << 8; // calculating 23 bit extended mantissa // make mantissa longer up to the 64, << 18

	if (mres == 0) return res;
	int32_t etmp = eres;
	eres -= std::max(std::min((__builtin_clz(mres) - 9), etmp), 0);
	mres <<= std::max(std::min((__builtin_clz(mres) - 9), etmp), 0);
	eres += (9 - __builtin_clz(mres)) * (mres >= 0x80'0000);
	mres >>= (9 - __builtin_clz(mres)) * (mres >= 0x80'0000);

	if (eres > 0) { // if normal
		mres = (mres >> 15) + ((mres & 0x7FFF) > 0x4000) + ((mres & 0xFFFF) == 0xC000); // last 23 bit stored 
		if (mres >= 0x0100) { // mres is greater than 2^23 (2^24) // should I?
			eres += 1;
			mres >>= 1;
		}
		if (eres >= 0xFF) { // if infinity
			return res | 0x7F80;
		}
		else {
			return res + uint32_t(mres) - 0x0080 + (eres << 7); // calculating result
		}
	}
	else if (eres >= -7) { // if subnormal
		uint32_t shift = 1u << (15 - eres + 1);
		return res + (mres >> (15 - eres + 1))
			+ ((mres & (shift - 1)) > (shift >> 1))
			+ (((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
	}

	return res;
}

uint16_t BF16::fma(uint16_t a, uint16_t b, uint16_t c) const noexcept {

	uint16_t res = (a ^ b) & 0x8000;
	int16_t ea = a & 0x7F80;
	int16_t eb = b & 0x7F80;
	int16_t ec = c & 0x7F80;
	int16_t eres;
	uint32_t ma = a & 0x007F;
	uint32_t mb = b & 0x007F;
	int32_t mc = c & 0x007F;
	int32_t mres;

	// NAN AND INF
	if ((ea == 0x7F80) || (eb == 0x7F80) || (ec == 0x7F80)) {
		return add(mul(a, b), c);
	}
	// MULTIPLYING
	eres = (ea >> 7) + (eb >> 7) + (ea == 0) + (eb == 0); // calculating exponent
	eres -= 127;
	mres = (((ma + 0x0080 * (ea > 0))) * ((mb + 0x0080 * (eb > 0)))) << 13; // calculating +12 bit extended mantissa

	// MUL RESULT TO NORMAL
	if (mres == 0) { 
		eres = 0;
	}

	// ADDING PREPARATIONS
	mc = (mc + (int32_t(ec > 0) << 7)) << 20;
	mres *= -(2 * int32_t(res >> 15) - 1); // make them signed
	mc *= -(2 * int32_t(c >> 15) - 1);
	eb = eres;
	ec >>= 7;
	bool ecDenormal = (ec == 0);
	ec += ecDenormal;

	// ADDING
	if (eb >= ec) { // calulate exponent and mantissa making exponents equal each other
		eres = eb;
		if (eb - ec >= 32) mres = mres; // 25
		else { 
			mres = mres + (mc >> (eb - ec));
		}
	}
	else {
		eres = ec;
		if (ec - eb >= 32) mres = mc;
		else {
			mres = mc + (mres >> (ec - eb));
		}
	}
	res = 0x8000 * (mres < 0); // creating res sign
	bool resSign = (mres < 0);
	mres *= -(2 * (mres < 0) - 1); // abs of mres

	// ADD RES NORMALIZATION
	if (mres >= 0x1000'0000) { // if mres is greater than a 2^23 (2^48) // make a non-if code 0x1'0000'0000'0000
		if ((mres >> 28) == 2) {
			mres >>= 2;
			eres += 2;
		}
		else {
			mres >>= 1;
			eres += 1;
		}
	}
	else if (mres == 0)
		eres = 0;
	else if (__builtin_clz(mres) > 4 && eres > 0 && (__builtin_clz(mres) - 4) < eres) { // if mres is less than a 2^23 (2^24) // can add mres == 0
		eres -= (__builtin_clz(mres) - 4);
		mres <<= (__builtin_clz(mres) - 4);
	}
	else if (__builtin_clz(mres) > 4 && eres > 0 && (__builtin_clz(mres) - 4) >= eres) {
		mres <<= eres - 1;
		eres = 0;
	}

	if (eres >= 0xFF) { // in inf
		return res + 0x7F80;
	}
	else if (eres > 0) { // use here all your saved stuff (mresLShift)
		mres = (mres >> 20) + ((mres & 0xF'FFFF) > 0x8'0000) + ((mres & 0x1F'FFFF) == 0x18'0000);
		return res + mres - 0x0080 + (eres << 7);
	}
	else if (eres >= -7) {
		uint32_t shift = 1u << (20 - eres);
		return res + (mres >> (20 - eres)) + ((mres & (shift - 1)) > (shift >> 1)) + (((mres & ((shift << 1) - 1))) == (shift + (shift >> 1)));
	}
	else {
		return res;
	}
	return res;
}

uint16_t BF16::div(uint16_t l, uint16_t r) const noexcept {
	uint16_t res = (l ^ r) & 0x8000;
	uint16_t el = l & 0x7F80;
	uint16_t er = r & 0x7F80;
	int16_t eres = 0;
	uint16_t ml = l & 0x007F;
	uint16_t mr = r & 0x007F;

	if ((el == 0x7F80) || (er == 0x7F80) || (r == 0) || (r == 0x8000) || (l == 0) || (l == 0x8000)) { // nan and inf
		if (el == 0x7F80) { // left is inf or nan
			if (ml != 0) // left is nan
				return res + (l & 0x7FFF);
			else if (er == 0x7F80) // right is inf or nan
				return res + 0x7F81; // inf / inf = nan
			else
				return res + (l & 0x7FFF); // inf / finite = inf
		}
		if (er == 0x7F80) { // right is inf or nan
			if (mr != 0) // is nan
				return res + (r & 0x7FFF);
			else // is inf
				return res; // zero
		}

		if ((r == 0) || (r == 0x8000)) { // right is zero
			if ((l == 0) || (l == 0x8000)) // left is zero
				return res | 0x7F81;
			else
				return res | 0x7F80;
		}

		if ((l == 0) || (l == 0x8000)) {
			return res;
		}
	}


	eres = (el >> 7) + 126; // diapason [1, 2)
	if (er == 0) {
		--eres; // er += 1 if right is subnormal
		eres += __builtin_clz(uint32_t(mr)) - 24;
		mr <<= __builtin_clz(uint32_t(mr)) - 24;
		/*
		while (mr < 0x80) { // right to normal
			++eres; // --er
			mr <<= 1;
		}
		*/
		mr -= 0x80;
	}
	eres -= (er >> 7);
	if (el == 0) {
		++eres; // el += 1 if left is subnormal
		eres -= __builtin_clz(uint32_t(ml)) - 24;
		ml <<= __builtin_clz(uint32_t(ml)) - 24;
		ml -= 0x80;
	}
	r = mr + (127u << 7); // diapason [1, 2)
	l = ml + (127u << 7);

	uint16_t x = fma((uint16_t)0xbef1u, r, (uint16_t)0x3fb5u); // 24/17 - 8/17 * d.
	r = usub(r);
	uint16_t e;

	e = fma(r, x, (uint16_t)0x3F80u); // 1 iteration
	x = fma(x, e, x);
	e = fma(r, x, (uint16_t)0x3F80u); // 2 iteration
	x = fma(x, e, x);
	e = fma(r, x, (uint16_t)0x3F80u); // 3 iteration AND RN(1/a) (EXCEPT FOR UNTEGRAL SIGNIFICAND 11...11)
	x = fma(x, e, x);

	uint16_t y = mul(l, x); // Probably, 1 ulp error. Testing required. y = l / r
	e = fma(r, y, l);
	y = fma(e, x, y);
	e = fma(r, y, l);
	uint16_t tmpy = fma(e, x, y);
	bool dy = y != tmpy; // for possible double-roundings
	bool de = (e >> 15) == 0x1;
	y = tmpy;

	r = usub(r);
	eres += (y >> 7) - 126; // overflow if y >= 1

	if (eres >= 255)
		res += 0x7F80;
	else if (eres > 0) {
		res += (eres << 7) + (y & 0x007F);
	}
	else if (eres >= -7) { // PAGE 137, r shows the rounding
		y = (y & 0x007F) + 0x0080;
		uint32_t shift = 1u << (-eres + 1);
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

BF16::BF16(uint16_t t) noexcept : data(t) {}
BF16::BF16(float f) noexcept {
	uint32_t tmp = *reinterpret_cast<uint32_t*>(&f);
	data = (tmp >> 16) + ((tmp & 0xFFFF) > 0x8000) + ((tmp & 0x1'FFFF) == 0x1'8000);
}
BF16::BF16(double d) noexcept {
	float f = d;
	uint32_t tmp = *reinterpret_cast<uint32_t*>(&f);
	data = (tmp >> 16) + ((tmp & 0xFFFF) > 0x8000) + ((tmp & 0x1'FFFF) == 0x1'8000);
}
BF16::BF16(const BF16& fp) noexcept : data(fp.data) {}
BF16& BF16::operator=(const BF16& fp) noexcept {
	data = fp.data;
	return *this;
}
BF16& BF16::operator=(uint16_t t) noexcept {
	data = t;
	return *this;
}
BF16& BF16::operator=(float f) noexcept {
	data = BF16(f).data;
	return *this;
}
BF16& BF16::operator=(double f) noexcept {
	data = BF16(f).data;
	return *this;
}
BF16::operator float() const noexcept {
	uint32_t tmp = data;
	tmp <<= 16;
	return reinterpret_cast<float*>(&tmp)[0];
}
BF16::operator double() const noexcept {
	uint32_t tmp = data;
	tmp <<= 16;
	return double(reinterpret_cast<float*>(&tmp)[0]);
}

BF16 BF16::operator- () noexcept {
	data = usub(data);
	return *this;
}
BF16& BF16::operator+= (const BF16& fp) noexcept {
	data = add(data, fp.data);
	return *this;
}
BF16& BF16::operator-= (const BF16& fp) noexcept {
	data = sub(data, fp.data);
	return *this;
}
BF16& BF16::operator*= (const BF16& fp) noexcept {
	data = mul(data, fp.data);
	return *this;
}
BF16& BF16::operator/= (const BF16& fp) noexcept {
	data = div(data, fp.data);
	return *this;
}
BF16 BF16::operator+ (const BF16& fp) const noexcept {
	BF16 tmp = *this;
	tmp += fp;
	return tmp;
}
BF16 BF16::operator- (const BF16& fp) const noexcept {
	BF16 tmp = *this;
	tmp -= fp;
	return tmp;
}
BF16 BF16::operator* (const BF16& fp) const noexcept {
	BF16 tmp = *this;
	tmp *= fp;
	return tmp;
}
BF16 BF16::operator/ (const BF16& fp) const noexcept {
	BF16 tmp = *this;
	tmp /= fp;
	return tmp;
}
BF16 BF16::fma (const BF16& a, const BF16& b) const noexcept {
	BF16 tmp = *this;
	tmp = fma(data, a.data, b.data);
	return tmp;
}

std::ostream& operator << (std::ostream& stream, const BF16& fp) {
	float tmp = fp;
	stream << tmp;
	return stream;
}
std::istream& operator >> (std::istream& stream, BF16& fp) {
	float tmp;
	stream >> tmp;
	fp = tmp;
	return stream;
}

bool BF16::operator == (const BF16& fp) const noexcept {
	if (((data & 0x7FFF) > 0x7F80) || ((fp.data & 0x7FFF) > 0x7F80)) return false;
	return data == fp.data;
}
bool BF16::operator < (const BF16& fp) const noexcept {
	if (((data & 0x7FFF) > 0x7F80) || ((fp.data & 0x7FFF) > 0x7F80)) return false;
	BF16 tmp = *this - fp;
	return tmp.data > 0x8000u;
}
bool BF16::operator > (const BF16& fp) const noexcept {
	if (((data & 0x7FFF) > 0x7F80) || ((fp.data & 0x7FFF) > 0x7F80)) return false;
	BF16 tmp = *this - fp;
	return (tmp.data < 0x8000u) && (tmp.data != 0u);
}
bool BF16::operator >= (const BF16& fp) const noexcept {
	if (((data & 0x7FFF) > 0x7F80) || ((fp.data & 0x7FFF) > 0x7F80)) return false;
	BF16 tmp = *this - fp;
	return (tmp.data <= 0x8000u);
}
bool BF16::operator <= (const BF16& fp) const noexcept {
	if (((data & 0x7FFF) > 0x7F80) || ((fp.data & 0x7FFF) > 0x7F80)) return false;
	BF16 tmp = *this - fp;
	return (tmp.data >= 0x8000u) && (tmp.data == 0u);
}
bool BF16::operator != (const BF16& fp) const noexcept {
	if (((data & 0x7FFF) > 0x7F80) || ((fp.data & 0x7FFF) > 0x7F80)) return true;
	return data != fp.data;
}