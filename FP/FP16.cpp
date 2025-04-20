#include "FP16.h"

bool FP16::IsSubnormal() const noexcept { return (exp == 0) && (man != 0); }
bool FP16::IsInf() const noexcept { return exp == ((1 << expLength) - 1) && (man == 0); }
bool FP16::IsNull() const noexcept { return (exp == 0) && (man == 0); }
bool FP16::IsNan() const noexcept { return ((exp == ((1 << expLength) - 1)) && (man != 0)); }
void FP16::printFP16_bites() const noexcept{
    std::cout << sign << "_";
    for (char i = 0; i < expLength; ++i) std::cout << ((exp >> (expLength - i - 1)) & 1);
    std::cout << "_";
    for (char i = 0; i < manLength; ++i) std::cout << ((man >> (manLength - i - 1)) & 1);
}

/*FP16::operator float16_t() const noexcept {
    unsigned short bits = (sign << (manLength + expLength)) | (exp << manLength) | man;
    return *reinterpret_cast<float16_t*>(&bits);
}*/

FP16& FP16::operator=(FP16& N) noexcept {
    if (this != &N) {
        sign = N.sign;
        exp = N.exp;
        man = N.man;
    }
    return *this;
}

void FP16::set_from_uint(uint16_t temp) noexcept{
    sign = (temp >> 15) & 0x1;
    exp = (temp >> 10) & 0x1F;
    man = temp & 0x3FF;
}

uint16_t FP16::get_int() const noexcept {
    uint16_t temp = (sign << (manLength + expLength)) | (exp << manLength) | man;
    return temp;
}

FP16::FP16(uint16_t temp) noexcept {
    sign = temp >> (manLength + expLength);
    exp = (temp >> manLength) & ((1 << expLength) - 1);
    man = temp & ((1 << manLength) - 1);
}

FP16 FP16::operator+(FP16 right) const noexcept {
    if (IsNull()) return right;
    if (right.IsNull()) return *this;

    FP16 Res;
    uint16_t temp = 0;
    int8_t diff = exp - right.exp;
    int16_t shift;
    uint8_t max_exp;
    uint32_t lost_bit = 0;

    if (diff < 0) max_exp = right.exp;
    else max_exp = exp;

    if ((exp == right.exp)) {
        shift = 1 << (manLength);
        if (man >= right.man) {
            temp = (1 << manLength) + man + (1 - int(2 * ((this->sign + right.sign) % 2))) * (shift + right.man);
            Res.sign = sign;
        }
        else {
            temp = (1 << manLength) + right.man + (1 - int(2 * ((this->sign + right.sign) % 2))) * (man + shift);
            Res.sign = right.sign;
        }
    }

    if (this->exp > right.exp) {
        shift = 1 << (manLength - diff);
        if ((manLength - diff) < 0) shift = 0;
        temp = (1 << manLength) + this->man + (1 - int(2 * ((this->sign + right.sign) % 2))) * (shift + (right.man >> diff));
        lost_bit = right.man & (((1 << (diff)) - 1));
        lost_bit += (1 << manLength) & (((1 << (diff)) - 1));
        Res.sign = sign;
    }

    if (this->exp < right.exp) {
        diff *= -1;
        shift = 1 << (manLength - diff);
        if (manLength - diff < 0) shift = 0;
        temp = (1 << manLength) + right.man + (1 - int(2 * ((this->sign + right.sign) % 2))) * (shift + (this->man >> diff));
        lost_bit = man & (((1 << (diff)) - 1));
        lost_bit += (1 << manLength) & (((1 << (diff)) - 1));
        Res.sign = right.sign;
    }

    if (IsSubnormal() || right.IsSubnormal()) {
        if (IsSubnormal() && right.IsSubnormal()) {
            if (man >= right.man) {
                temp = man + (1 - (2 * ((this->sign + right.sign) % 2))) * right.man;
                Res.sign = sign;
            }
            else {
                temp = right.man + (1 - (2 * ((this->sign + right.sign) % 2))) * man;
                Res.sign = right.sign;
            }
        }
        else {
            //left subnormal
            if (IsSubnormal()) {
                diff = right.exp - 1;
                temp = (1 << manLength) + right.man + (1 - (2 * ((this->sign + right.sign) % 2))) * (this->man >> diff);
                lost_bit = man & (((1 << (diff)) - 1));
                Res.sign = right.sign;
                max_exp = right.exp;
            }
            //right subnormal
            else {
                diff = exp - 1;
                temp = (1 << manLength) + this->man + (1 - (2 * ((this->sign + right.sign) % 2))) * (right.man >> diff);
                lost_bit = right.man & (((1 << (diff)) - 1));
                Res.sign = sign;
                max_exp = exp;
            }
        }
    }

    if (sign != right.sign && lost_bit != 0) {
        temp--;
        lost_bit = (1 << (diff)) - lost_bit;
    }


    while (temp < (1 << manLength) && max_exp>1) {
        temp <<= 1;
        max_exp--;
        temp += (lost_bit >> (diff - 1)) & 1;
        lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
        diff--;
    }

    if (temp == 0) {
        Res.sign = 0;
        Res.exp = 0;
        Res.man = 0;
        return Res;
    }

    if (temp >= (1 << (manLength + 1)) && (max_exp > 0)) {
        lost_bit += ((temp & 1) << (diff));
        diff++;
        temp = temp >> 1;
        max_exp++;
    }


    if ((lost_bit >> (diff - 1)) & 1) {
        lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
        if (lost_bit > 0) {
            temp++;
        }
        else {
            temp += (temp & 1);
        }
        if ((temp >= (1 << (manLength + 1))) && (max_exp >= 1)) {
            temp = temp >> 1;
            max_exp++;
        }
    }


    if (temp < (1 << manLength) && max_exp == 1) {
        max_exp = 0;
    }

    if (temp >= (1 << manLength) && (max_exp == 0)) {
        max_exp = 1;
        temp -= (1 << manLength);
    }

    if ((max_exp >= ((1 << expLength) - 1)) || IsInf() || right.IsInf() || IsNan() || right.IsNan()) {
        //inf +- inf
        if (IsInf() && right.IsInf()) {
            if (sign != right.sign) {
                Res.exp = (1 << expLength);
                Res.man = 1;
                return Res;
            }
            else {
                return *this;
            }
        }

        if (IsInf()) return *this;

        if (right.IsInf()) return right;

        if (IsNan()) return *this;

        if (right.IsNan()) return right;

        //if result is inf
        if (max_exp >= ((1 << expLength) - 1)) {
            Res.exp = ((1 << expLength) - 1);
            Res.man = 0;
            return Res;
        }
    }


    Res.exp = max_exp;
    Res.man = (temp - (1 << manLength));
    if (max_exp == 0) Res.man = temp;
    return Res;
}

FP16 FP16::operator-(FP16 right) const noexcept {
    right.sign += 1;
    return *this + right;
}

FP16 FP16::operator*(FP16 right) const noexcept {
    FP16 Res;
    Res.sign = sign ^ right.sign;
    uint16_t temp_man;
    int16_t temp_exp;
    uint32_t lost_bit;
    uint8_t diff = manLength;
    if (IsNull() || right.IsNull()) {
        Res.exp = 0;
        Res.man = 0;
        return Res;
    }
    if (!IsSubnormal() && !right.IsSubnormal()) {
        temp_exp = exp + right.exp - shiftExp - shiftExp;
        temp_man = (1 << manLength) + man + right.man + uint16_t((uint32_t(man) * right.man) >> manLength);
    }
    else {
        if (IsSubnormal() && right.IsSubnormal()) {
            Res.exp = 0;
            Res.man = 0;
            return Res;
        }
        else {
            temp_exp = exp - shiftExp + right.exp - shiftExp + 1;
            if (IsSubnormal()) {
                temp_man = man + uint16_t((uint32_t(man) * right.man) >> manLength);
            }
            else {
                temp_man = right.man + uint16_t((uint32_t(man) * right.man) >> manLength);
            }
        }
    }

    lost_bit = (uint32_t(man) * right.man) & ((1 << manLength) - 1);

    if (temp_exp < -(shiftExp - 1)) {
        while (temp_exp < -(shiftExp - 1)) {
            lost_bit += (temp_man & 1) << diff;
            diff++;
            temp_man = temp_man >> 1;
            temp_exp++;
        }
        if (temp_man < (1 << manLength)) temp_exp = -shiftExp;
    }

    while (temp_man < (1 << manLength) && temp_exp>-(shiftExp - 1)) {
        temp_man <<= 1;
        temp_exp--;

        temp_man += (lost_bit >> (diff - 1)) & 1;
        lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
        diff--;
    }

    if ((temp_man >= (1 << (manLength + 1))) && (temp_exp >= -(shiftExp - 1))) {
        lost_bit += ((temp_man & 1) << (diff));
        diff++;
        temp_man = temp_man >> 1;
        temp_exp++;
    }

    if ((lost_bit >> (diff - 1)) & 1) {
        lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
        if (lost_bit > 0) {
            temp_man++;
        }
        else {
            temp_man += (temp_man & 1);
        }
        if ((temp_man >= (1 << (manLength + 1))) && (temp_exp >= -(shiftExp - 1))) {
            temp_man = temp_man >> 1;
            temp_exp++;
        }
    }

    if (temp_exp == -shiftExp && temp_man >= (1 << manLength)) {
        temp_exp++;
    }

    if ((temp_exp >= ((1 << expLength) - 1 - shiftExp)) || IsInf() || right.IsInf() || IsNan() || right.IsNan()) {
        //inf +- inf
        if (IsInf() && right.IsInf()) {
            if (sign != right.sign) {
                Res.exp = (1 << expLength);
                Res.man = 1;
                return Res;
            }
            else {
                return *this;
            }
        }

        if (IsInf()) return *this;

        if (right.IsInf()) return right;

        if (IsNan()) return *this;

        if (right.IsNan()) return right;

        //if result is inf
        if (temp_exp >= (((1 << expLength) - 1 - shiftExp))) {
            Res.exp = (1 << expLength) - 1;
            Res.man = 0;
            return Res;
        }
    }

    if ((temp_man < (1 << manLength)) && (temp_exp == -(shiftExp - 1))) {
        Res.exp = 0;
        Res.man = temp_man;
        return Res;
    }

    Res.exp = temp_exp + shiftExp;
    Res.man = temp_man - (1 << manLength);

    return Res;
}

uint16_t fma_half(uint16_t a, uint16_t b, uint16_t c) noexcept {
    FP16 A, B, C;
    A.set_from_uint(a);
    B.set_from_uint(b);
    C.set_from_uint(c);
    return (fma4(A, B, C).get_int());
}

uint16_t mul_half(uint16_t a, uint16_t b) noexcept {
    FP16 A, B;
    A.set_from_uint(a);
    B.set_from_uint(b);
    return ((A * B).get_int());
}

uint16_t usub_half(uint16_t temp) noexcept {
    if ((temp >> 15) == 1) temp &= ((1 << 15) - 1);
    else temp |= (1 << 15);
    return temp;
}


FP16 FP16::operator/(FP16 right) const noexcept {
    FP16 Res;
    uint16_t l = get_int();
    uint16_t r = right.get_int();
    uint16_t res = (l ^ r) & 0x8000;
    uint16_t el = l & 0x7C00;
    uint16_t er = r & 0x7C00;
    int16_t eres = 0;
    uint16_t ml = l & 0x03FF;
    uint16_t mr = r & 0x03FF;
    bool coutflag = false;

    // Handle special cases (NaN, Inf, zero)
    if ((el == 0x7C00) || (er == 0x7C00) || (r == 0) || (r == 0x8000) || (l == 0) || (l == 0x8000)) {
        if (el == 0x7C00) { // left is inf or nan
            if (ml != 0) {// left is nan
                Res.set_from_uint(res + (l & 0x7FFF));
                return Res;
            }
            else if (er == 0x7C00) { // right is inf or nan
                Res.set_from_uint(res + 0x7C01); // inf / inf = nan
                return Res;
            }
            else {
                Res.set_from_uint(res + (l & 0x7FFF));
                return Res;
            }
        }
        if (er == 0x7C00) { // right is inf or nan
            if (mr != 0) { // is nan
                Res.set_from_uint(res + (r & 0x7FFF));
                return Res;
            }
            else { // is inf
                Res.set_from_uint(res);
                return Res; // zero
            }
        }

        if ((r == 0) || (r == 0x8000)) { // right is zero
            if ((l == 0) || (l == 0x8000)) { // left is zero
                Res.set_from_uint(res | 0x7C01);
                return Res;
            }
            else
            {
                Res.set_from_uint(res | 0x7C00);
                return Res;
            }
        }

        if ((l == 0) || (l == 0x8000)) {
            Res.set_from_uint(res);
            return Res;
        }
    }


    eres = (el >> 10) + 14;
    if (er == 0) {
        --eres; // adjust for subnormal
        while (mr < 0x0400) { // normalize right
            ++eres;
            mr <<= 1;
        }
        mr -= 0x0400;
    }
    eres -= (er >> 10);
    if (el == 0) {
        ++eres; // adjust for subnormal
        while (ml < 0x0400) { // normalize left
            --eres;
            ml <<= 1;
        }
        ml -= 0x0400;
    }


    r = mr + (15u << 10); // FP16 bias = 15
    l = ml + (15u << 10);


    r = usub_half(r);
    uint16_t x = fma_half(14216, r, 15782); // 24/17 - 8/17*d
    uint16_t e;

    e = fma_half(r, x, (uint16_t)0x3C00); // 1 iteration
    x = fma_half(x, e, x);
    e = fma_half(r, x, (uint16_t)0x3C00); // 2 iteration
    x = fma_half(x, e, x);


    uint16_t y = mul_half(l, x);
    e = fma_half(r, y, l);
    y = fma_half(e, x, y);
    e = fma_half(r, y, l);
    uint16_t tmpy = fma_half(e, x, y);
    bool dy = y != tmpy;
    bool de = (e >> 15) == 0x1;
    y = tmpy;


    r = usub_half(r);
    eres += (y >> 10) - 14;


    if (eres >= 31)
        res += 0x7C00;
    else if (eres > 0) {
        res += (eres << 10) + (y & 0x03FF);
    }
    else if (eres >= -10) {
        y = (y & 0x03FF) + 0x0400;
        uint16_t shift = 1u << (-eres + 1);
        if (e == 0x0 || ((y & (shift - 1)) != (shift >> 1))) {
            res += (y >> (-eres + 1)) +
                ((y & (shift - 1)) > (shift >> 1)) +
                (((y & ((shift << 1) - 1))) == (shift + (shift >> 1)));
        }
        else if (dy && de || !dy && !de) {
            res += (y >> (-eres + 1)) + 1;
        }
        else if (dy && !de || !dy && de) {
            res += (y >> (-eres + 1)) + 0;
        }
    }

    Res.set_from_uint(res);
    if ((Res.exp != 31) && (((this->exp == 0) && (right.exp != 0)) || (this->man == 0)) && (right.man == 1023) && (Res.man == 0) && (Res.exp != 0)) Res.man++;
    return Res;
}

bool FP16::operator == (const FP16& right) const noexcept {
    return ((this->sign == right.sign) && (this->exp == right.exp) && (this->man == right.man));
}

bool FP16::operator != (const FP16& right) const noexcept {
    return (!((*this) == right));
}

bool FP16::operator < (const FP16& right) const noexcept {
    if (sign != right.sign) return(sign == 1);
    if (exp != right.exp) return (exp < right.exp);
    return(man < right.man);
}
bool FP16::operator > (const FP16& right) const noexcept {
    if (sign != right.sign) return(sign == 0);
    if (exp != right.exp) return (exp > right.exp);
    return(man > right.man);
}
bool FP16::operator >= (const FP16& right) const noexcept {
    return !((*this) < right);
}
bool FP16::operator <= (const FP16& right) const noexcept {
    return !((*this) > right);
}

FP16& FP16::operator+=(const FP16& right) noexcept {
    FP16 temp = (*this + right).sign;
    *this = temp;
    return *this;
}

FP16& FP16::operator-=(const FP16& right) noexcept {
    FP16 temp = (*this - right).sign;
    *this = temp;
    return *this;
}

FP16& FP16::operator*=(const FP16& right) noexcept {
    FP16 temp = (*this * right).sign;
    *this = temp;
    return *this;
}

FP16& FP16::operator/=(const FP16& right) noexcept {
    FP16 temp = (*this / right).sign;
    *this = temp;
    return *this;
}

FP16 FP16::operator-() noexcept {
    FP16 temp(*this);
    temp.sign = (temp.sign) ^ 1;
    return temp;
}

FP16 fma4(FP16 a, FP16 b, FP16 c) noexcept {
    FP16 Res;
    uint32_t mul_man;
    int16_t mul_exp;
    int64_t shift = manLength;
    uint64_t lostbit;
    Res.sign = a.sign ^ b.sign;

    mul_exp = a.exp + b.exp - shiftExp;
    if (a.exp != 0 && b.exp != 0) {
        mul_man = (((1 << (manLength)) + (a.man + b.man)) << manLength) + (uint64_t(a.man) * b.man);
    }
    else {
        if (a.exp == 0 && b.exp == 0) {
            mul_exp = 0;
            mul_man = 0;
            return c;
        }
        else {
            mul_exp = a.exp - shiftExp + b.exp + 1;
            if (a.exp == 0) {
                mul_man = (a.man << manLength) + (uint64_t(a.man) * b.man);
            }
            else {
                mul_man = (b.man << manLength) + (uint64_t(a.man) * b.man);
            }
        }
    }

    if (mul_man == 0) return c;

    if (mul_exp < 1) {
        while (mul_exp < 1) {
            shift++;
            mul_exp++;
        }
        if ((mul_man << shift) < (1 << manLength)) mul_exp = 0;
    }

    if ((mul_man >> shift) >= (1 << (1 + manLength))) {
        shift++;
        mul_exp++;
    }

    while ((mul_man >> shift) < (1 << (manLength)) && (mul_exp > 1)) {
        shift--;
        mul_exp--;
    }

    if (mul_exp == 1 && ((mul_man >> shift) < (1 << (manLength)))) {
        mul_exp = 0;
    }

    if (mul_exp >= ((1 << expLength) - 1)) {
        Res.exp = ((1 << expLength) - 1);
        Res.man = 0;
        return Res;
    }
    bool mul_sign = Res.sign;
    uint32_t g_exp, l_exp, g_man, l_man;
    int64_t temp_man;
    uint16_t temp_exp;
    uint32_t temp;
    temp_man = mul_man >> shift;
    bool flag = 0;
    if (mul_exp > c.exp) {
        g_exp = mul_exp;
        g_man = mul_man;
        l_man = c.man;
        if (c.exp != 0) l_man += (1 << manLength);
        l_man = l_man << shift;
        l_exp = c.exp;
    }
    else {
        if (mul_exp < c.exp) {
            l_exp = mul_exp;
            l_man = mul_man;
            g_man = c.man;
            if (c.exp != 0) g_man += (1 << manLength);
            g_man = g_man << shift;
            g_exp = c.exp;
            Res.sign = c.sign;
        }
        else {
            temp = c.man;
            if (c.exp > 0) temp += (1 << manLength);
            temp = temp << shift;
            if (temp > mul_man) {
                l_exp = mul_exp;
                l_man = mul_man;
                g_man = temp;
                g_exp = c.exp;
                Res.sign = c.sign;
            }
            else {
                g_exp = mul_exp;
                g_man = mul_man;
                l_man = temp;
                l_exp = c.exp;
            }
        }
    }

    if (l_exp == 0) l_exp = 1;
    if (g_exp == 0) g_exp = 1;
    temp_exp = g_exp;
    if (g_exp != l_exp)   temp_man = (int64_t(g_man) << (1 << expLength)) + ((int(1) - int(2 * ((mul_sign + c.sign) % 2))) * (int64_t(l_man) << ((1 << expLength) - (g_exp - l_exp))));
    else temp_man = (int64_t(g_man) << (1 << expLength)) + ((int(1) - int(2 * ((mul_sign + c.sign) % 2))) * (int64_t(l_man) << (1 << expLength)));
    shift += (1 << expLength);
    while (((temp_man >> shift) < (1 << manLength)) && temp_exp > 1) {
        temp_exp--;
        shift--;
    }
    if (((temp_man >> shift) >= (1 << (manLength + 1)) && (temp_exp > 0))) {
        shift++;
        temp_exp++;
    }
    lostbit = temp_man & ((uint64_t(1) << (shift)) - 1);
    temp_man >>= shift;
    if ((lostbit >> (shift - 1)) & 1) {
        lostbit = lostbit & ((int64_t(1) << (shift - 1)) - 1);
        if (lostbit > 0) temp_man++;
        else {
            temp_man += (temp_man & 1);
        }
        if (temp_man == (1 << (manLength + 1))) {
            temp_man >>= 1;
            temp_exp++;
        }
    }

    if ((temp_exp == 1) && (temp_man < (1 << manLength))) temp_exp = 0;

    Res.man = temp_man;
    Res.exp = temp_exp;

    if (temp_exp >= ((1 << expLength) - 1)) {
        Res.man = 0;
        Res.exp = ((1 << expLength) - 1);
    }

    if (temp_man == 0) {
        Res.exp = 0;
        Res.man = 0;
        Res.sign = c.sign;
    }

    return Res;
}

float fp16_to_float(uint16_t h) noexcept {
    uint32_t sign = (h >> 15) & 0x1;
    uint32_t exponent = (h >> 10) & 0x1F;
    uint32_t mantissa = h & 0x3FF;
    uint32_t f;
    if (exponent == 0) {
        if (mantissa == 0) {
            f = sign << 31;
        }
        else {
            exponent = 127 - 14;
            while ((mantissa & 0x400) == 0) {
                mantissa <<= 1;
                exponent--;
            }
            mantissa &= 0x3FF;
            f = (sign << 31) | (exponent << 23) | (mantissa << 13);
        }
    }
    else if (exponent == 0x1F) {
        f = (sign << 31) | 0x7F800000 | (mantissa << 13);
    }
    else {
        exponent += (127 - 15);
        f = (sign << 31) | (exponent << 23) | (mantissa << 13);
    }
    return *reinterpret_cast<float*>(&f);
}

std::ostream& operator<<(std::ostream& os, FP16 num) noexcept {
    float f = fp16_to_float(num.get_int());
    os << f;
    return os;
}


