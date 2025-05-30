#include "FP16.h"

bool FP16::IsSubnormal() const noexcept { return (exp == 0) && (man != 0); }
bool FP16::IsInf() const noexcept { return exp == ((1 << 5) - 1) && (man == 0); }
bool FP16::IsNull() const noexcept { return (exp == 0) && (man == 0); }
bool FP16::IsNan() const noexcept { return ((exp == ((1 << 5) - 1)) && (man != 0)); }
void FP16::printFP16_bites() const noexcept{
    std::cout << sign << "_";
    for (char i = 0; i < 5; ++i) std::cout << ((exp >> (5 - i - 1)) & 1);
    std::cout << "_";
    for (char i = 0; i < 10; ++i) std::cout << ((man >> (10 - i - 1)) & 1);
}

/*FP16::operator float16_t() const noexcept {
    unsigned short bits = (sign << (10 + 5)) | (exp << 10) | man;
    return *reinterpret_cast<float16_t*>(&bits);
}*/

FP16& FP16::operator=(const FP16& N) noexcept {
    if (this != &N) {
        sign = N.sign;
        exp = N.exp;
        man = N.man;
    }
    return *this;
}

void FP16::set_from_uint(const uint16_t& temp) noexcept{
    sign = (temp >> 15) & 0x1;
    exp = (temp >> 10) & 0x1F;
    man = temp & 0x3FF;
}

uint16_t FP16::get_int() const noexcept {
    uint16_t temp = (sign << (10 + 5)) | (exp << 10) | man;
    return temp;
}

//FP16::FP16(uint16_t temp) noexcept {
//    sign = temp >> (10 + 5);
//    exp = (temp >> 10) & ((1 << 5) - 1);
//    man = temp & ((1 << 10) - 1);
//}


FP16 FP16::operator-(FP16 right) const noexcept {
    right.sign += 1;
    return *this + right;
}

FP16 FP16::operator*(const FP16& right) const noexcept {
    FP16 Res;
    Res.sign = sign ^ right.sign;
    uint16_t temp_man;
    int16_t temp_exp;
    uint32_t lost_bit;
    uint8_t diff = 10;
    if (IsNull() || right.IsNull()) {
        Res.exp = 0;
        Res.man = 0;
        return Res;
    }
    if (!IsSubnormal() && !right.IsSubnormal()) {
        temp_exp = exp + right.exp - 15 - 15;
        temp_man = (1 << 10) + man + right.man + uint16_t((uint32_t(man) * right.man) >> 10);
    }
    else {
        if (IsSubnormal() && right.IsSubnormal()) {
            Res.exp = 0;
            Res.man = 0;
            return Res;
        }
        else {
            temp_exp = exp - 15 + right.exp - 15 + 1;
            if (IsSubnormal()) {
                temp_man = man + uint16_t((uint32_t(man) * right.man) >> 10);
            }
            else {
                temp_man = right.man + uint16_t((uint32_t(man) * right.man) >> 10);
            }
        }
    }

    lost_bit = (uint32_t(man) * right.man) & ((1 << 10) - 1);

    if (temp_exp < -(15 - 1)) {
        while (temp_exp < -(15 - 1)) {
            lost_bit += (temp_man & 1) << diff;
            diff++;
            temp_man = temp_man >> 1;
            temp_exp++;
        }
        if (temp_man < (1 << 10)) temp_exp = -15;
    }

    if ((temp_man >= (1 << (10 + 1))) && (temp_exp >= -(15 - 1))) {
        lost_bit += ((temp_man & 1) << (diff));
        diff++;
        temp_man = temp_man >> 1;
        temp_exp++;
    }
    else {
        while (temp_man < (1 << 10) && temp_exp>-(15 - 1)) {
            temp_man <<= 1;
            temp_exp--;

            temp_man += (lost_bit >> (diff - 1)) & 1;
            lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
            diff--;
        }
    }


    if ((lost_bit >> (diff - 1)) & 1) {
        lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
        if (lost_bit > 0) {
            temp_man++;
        }
        else {
            temp_man += (temp_man & 1);
        }
        if ((temp_man >= (1 << (10 + 1))) && (temp_exp >= -(15 - 1))) {
            temp_man = temp_man >> 1;
            temp_exp++;
        }
    }

    if (temp_exp == -15 && temp_man >= (1 << 10)) {
        temp_exp++;
    }

    if ((temp_exp >= (32 - 1 - 15)) || IsInf() || right.IsInf() || IsNan() || right.IsNan()) {
        //inf +- inf
        if (IsInf() && right.IsInf()) {
            if (sign != right.sign) {
                Res.exp = (1 << 5) - 1;
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
        if (temp_exp >= (((1 << 5) - 1 - 15))) {
            Res.exp = (1 << 5) - 1;
            Res.man = 0;
            return Res;
        }
    }

    if ((temp_man < (1 << 10)) && (temp_exp == -(15 - 1))) {
        Res.exp = 0;
        Res.man = temp_man;
        return Res;
    }

    Res.exp = temp_exp + 15;
    Res.man = temp_man - (1 << 10);

    return Res;
}

uint16_t fma_half(const uint16_t& a, const uint16_t& b, const uint16_t& c) noexcept {
    FP16 A, B, C;
    A.set_from_uint(a);
    B.set_from_uint(b);
    C.set_from_uint(c);
    return (fma4(A, B, C).get_int());
}

uint16_t mul_half(const uint16_t& a, const uint16_t& b) noexcept {
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


FP16 FP16::operator/(const FP16& right) const noexcept {
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
    FP16 temp = (*this + right);
    *this = temp;
    return *this;
}

FP16& FP16::operator-=(const FP16& right) noexcept {
    FP16 temp = (*this - right);
    *this = temp;
    return *this;
}

FP16& FP16::operator*=(const FP16& right) noexcept {
    FP16 temp = (*this * right);
    *this = temp;
    return *this;
}

FP16& FP16::operator/=(const FP16& right) noexcept {
    FP16 temp = (*this / right);
    *this = temp;
    return *this;
}

FP16 FP16::operator-() noexcept {
    FP16 temp(*this);
    temp.sign = (temp.sign) ^ 1;
    return temp;
}

FP16 fma4(const FP16& a, const FP16& b, const FP16& c) noexcept {
    FP16 Res;
    uint32_t mul_man;
    int16_t mul_exp;
    int64_t shift = 10;
    uint64_t lostbit;
    Res.sign = a.sign ^ b.sign;

    mul_exp = a.exp + b.exp - 15;
    if (a.exp != 0 && b.exp != 0) {
        mul_man = ((1024 + (uint64_t(a.man) + uint64_t(b.man))) << 10) + (uint64_t(a.man) * b.man);
    }
    else {
        if (a.exp == 0 && b.exp == 0) {
            mul_exp = 0;
            mul_man = 0;
            return c;
        }
        else {
            mul_exp = a.exp - 15 + b.exp + 1;
            if (a.exp == 0) {
                mul_man = (uint64_t(a.man) << 10) + (uint64_t(a.man) * b.man);
            }
            else {
                mul_man = (uint64_t(b.man) << 10) + (uint64_t(a.man) * b.man);
            }
        }
    }

    if (mul_man == 0) return c;

    if (mul_exp < 1) {
        while (mul_exp < 1) {
            shift++;
            mul_exp++;
        }
        if ((mul_man << shift) < (1 << 10)) mul_exp = 0;
    }

    if ((mul_man >> shift) >= (1 << (1 + 10))) {
        shift++;
        mul_exp++;
    }
    else {
        while ((mul_man >> shift) < (1 << (10)) && (mul_exp > 1)) {
            shift--;
            mul_exp--;
        }
    }

    if (mul_exp == 1 && ((mul_man >> shift) < (1 << (10)))) {
        mul_exp = 0;
    }

    if (mul_exp >= ((1 << 5) - 1)) {
        Res.exp = ((1 << 5) - 1);
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
        if (c.exp != 0) l_man += (1 << 10);
        l_man = l_man << shift;
        l_exp = c.exp;
    }
    else {
        if (mul_exp < c.exp) {
            l_exp = mul_exp;
            l_man = mul_man;
            g_man = c.man;
            if (c.exp != 0) g_man += (1 << 10);
            g_man = g_man << shift;
            g_exp = c.exp;
            Res.sign = c.sign;
        }
        else {
            temp = c.man;
            if (c.exp > 0) temp += (1 << 10);
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
    if (g_exp != l_exp)   temp_man = (int64_t(g_man) << (1 << 5)) + ((int(1) - int(2 * ((mul_sign + c.sign) % 2))) * (int64_t(l_man) << ((1 << 5) - (g_exp - l_exp))));
    else temp_man = (int64_t(g_man) << (1 << 5)) + ((int(1) - int(2 * ((mul_sign + c.sign) % 2))) * (int64_t(l_man) << (1 << 5)));
    shift += (1 << 5);
    while (((temp_man >> shift) < (1 << 10)) && temp_exp > 1) {
        temp_exp--;
        shift--;
    }
    if (((temp_man >> shift) >= (1 << (10 + 1)) && (temp_exp > 0))) {
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
        if (temp_man == (1 << (10 + 1))) {
            temp_man >>= 1;
            temp_exp++;
        }
    }

    if ((temp_exp == 1) && (temp_man < (1 << 10))) temp_exp = 0;

    Res.man = temp_man;
    Res.exp = temp_exp;

    if (temp_exp >= ((1 << 5) - 1)) {
        Res.man = 0;
        Res.exp = ((1 << 5) - 1);
    }

    if (temp_man == 0) {
        Res.exp = 0;
        Res.man = 0;
        Res.sign = c.sign;
    }

    return Res;
}

float fp16_to_float(const uint16_t& h) noexcept {
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


FP16& FP16::operator=(const float& N) noexcept {
    if (N == 0.0f) {
        sign = std::signbit(N) ? 1 : 0;
        exp = 0;
        man = 0;
        return *this;
    }

    uint32_t float_bits;
    std::memcpy(&float_bits, &N, sizeof(float));
    
    sign = (float_bits >> 31) & 0x1;
    int32_t float_exp = (float_bits >> 23) & 0xFF;
    uint32_t float_man = float_bits & 0x7FFFFF;


    if (float_exp == 0xFF) {
        exp = (1 << 5) - 1;
        man = float_man ? 1 : 0;
        return *this;
    }


    if (float_exp == 0) {
        float_exp = -126;

    } else {
        float_exp -= 127;
        float_man |= 0x800000;
    }

    int32_t new_exp = float_exp + 15;
    uint32_t new_man = float_man >> 13;


    if (new_exp <= 0) {
        if (new_exp < -10) {
            exp = 0;
            man = 0;
            return *this;
        }


        exp = 0;
        uint32_t shift = 1 - new_exp;
        new_man = float_man >> (shift + 13);
        

        uint32_t round_bit = (float_man >> (shift + 12)) & 1;
        uint32_t sticky = (float_man << (32 - (shift + 12))) != 0;
        new_man += (round_bit & (new_man | sticky));

        if (new_man >= (1 << 10)) {
            new_man = 0;
            exp = 1;
        }

        man = new_man;
        return *this;
    }


    exp = new_exp;
    man = new_man;


    uint32_t round_bit = (float_man >> 12) & 1;
    uint32_t sticky = (float_man & 0xFFF) != 0;
    man += (round_bit & (man | sticky));


    if (man >= (1 << 10)) {
        man >>= 1;
        exp++;
        if (exp >= ((1 << 5) - 1)) {
            exp = (1 << 5) - 1;
            man = 0;
        }
    }

    return *this;
}

FP16& FP16::operator=(const double& N) noexcept{
    *this = float(N);
    return *this;
}

FP16::operator double() const noexcept {
    double result = float(*this);
    return result;
}

FP16::operator float() const noexcept {
    float result = fp16_to_float(get_int());
    return result;
}


FP16 FP16::operator+(const FP16& right) const noexcept {
    FP16 Res;
    uint16_t temp = 0;
    int16_t shift;
    uint16_t max_exp,max_man,min_exp,min_man;
    uint32_t lost_bit = 0;
    uint8_t diff;

    if (exp > right.exp) {
        max_exp = exp;
        max_man = man;
        min_exp = right.exp;
        min_man = right.man;
        Res.sign = sign;
        diff = max_exp - min_exp;
    }
    else {
        if (exp == right.exp) {
            if (man >= right.man) {
                max_exp = exp;
                max_man = man;
                min_exp = right.exp;
                min_man = right.man;
                Res.sign = sign;
                diff = max_exp - min_exp;
            }
            else {
                min_exp = exp;
                min_man = man;
                max_exp = right.exp;
                max_man = right.man;
                Res.sign = right.sign;
                diff = max_exp - min_exp;
            }
        }
        else {
            min_exp = exp;
            min_man = man;
            max_exp = right.exp;
            max_man = right.man;
            Res.sign = right.sign;
            diff = max_exp - min_exp;
        }
    }
    

    shift = 1 << (10 - diff);
    if ((10 - diff) < 0) shift = 0;
    temp = (1024) + max_man + (1 - int(2 * ((this->sign + right.sign) % 2))) * (shift + (min_man >> diff));
    lost_bit = min_man & (((1 << (diff)) - 1));
    lost_bit += (1024) & (((1 << (diff)) - 1));


    if ((exp == 0) || (right.exp == 0)) {
        if ((exp == 0) && (right.exp == 0)) {
            temp = max_man + (1 - (2 * ((this->sign + right.sign) % 2))) * min_man;
        }
        else {
                diff = max_exp - 1;
                temp = (1024) + max_man + (1 - (2 * ((this->sign + right.sign) % 2))) * (min_man >> diff);
                lost_bit = min_man & (((1 << (diff)) - 1));
        }
    }


    if ((sign != right.sign) && (lost_bit != 0)) {
        temp--;
        lost_bit = (1 << (diff)) - lost_bit;
    }

    if (temp >= (2048) && (max_exp > 0)) {
        lost_bit += ((temp & 1) << (diff));
        diff++;
        temp = temp >> 1;
        max_exp++;
    }
    else {
        while (temp < (1024) && max_exp>1) {
            temp <<= 1;
            max_exp--;
            temp += (lost_bit >> (diff - 1)) & 1;
            lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
            diff--;
        }
    }
    if (temp == 0) {
        Res.sign = 0;
        Res.exp = 0;
        Res.man = 0;
        return Res;
    }
    if ((lost_bit >> (diff - 1))) {
        lost_bit = lost_bit & ((1 << (diff - 1)) - 1);
        if (lost_bit > 0) {
            temp++;
        }
        else {
            temp += (temp & 1);
        }
        if ((temp >= (1 << (11))) && (max_exp >= 1)) {
            temp = temp >> 1;
            max_exp++;
        }
    }

    if (max_exp <= 1) {
        if (temp < (1024) && max_exp == 1) {
            max_exp = 0;
        }

        if (temp >= (1024) && (max_exp == 0)) {
            max_exp = 1;
            temp -= (1024);
        }
    }

    if ((max_exp >= 31) || (exp == 31) || (right.exp == 31)) {
        //inf +- inf
        if (IsInf() && right.IsInf()) {
            if (sign != right.sign) {
                Res.exp = (1 << 5) - 1;
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
        if (max_exp >= 31) {
            Res.exp = 31;
            Res.man = 0;
            return Res;
        }
    }


    Res.exp = max_exp;
    if (max_exp == 0) Res.man = temp;
    else Res.man = (temp - (1024));
    return Res;
}


