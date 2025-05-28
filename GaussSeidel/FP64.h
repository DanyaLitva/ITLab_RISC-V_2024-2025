//#ifdef _MSC_VER
#include <iostream>
#include <bitset>
//#include <intrin.h>
#include <vector>
#define ull unsigned long long
#define ll long long
#define uint32 (unsigned int)

using namespace std;

#include <vector>
#include <string>
#include <cstring>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <math.h>


//#define USE_MPFR
#ifdef USE_MPFR
#include <gmp.h>
#include <mpfr.h>
#endif

double set_mantissa(double x, uint64_t new_mantissa) {
    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(double));
    constexpr uint64_t mantissa_mask = (1ULL << 52) - 1;
    bits &= ~mantissa_mask;
    bits |= (new_mantissa & mantissa_mask);
    std::memcpy(&x, &bits, sizeof(double));
    return x;
}

long long get_mantissa(double x) {
    if (x == 0)
        return 0;
    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(double));
    uint64_t mantissa_mask = 0x000FFFFFFFFFFFFFULL;
    return (bits & mantissa_mask);
}

double set_exponent(double x, uint16_t new_exponent) {
    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(double));
    constexpr uint64_t exponent_mask = (0x7FFULL << 52);
    bits &= ~exponent_mask;
    bits |= (static_cast<uint64_t>(new_exponent & 0x7FF)) << 52;
    std::memcpy(&x, &bits, sizeof(double));
    return x;
}

int get_exponent(double x) {
    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(double));
    uint64_t exponent_bits = (bits >> 52) & 0x7FF;
    return static_cast<int>(exponent_bits);
}

void printDoubleBits(double value) {
    uint64_t bits = *reinterpret_cast<uint64_t*>(&value);
    uint64_t sign = (bits >> 63) & 0x1;
    int64_t exponent = (bits >> 52) & 0x7FF;
    uint64_t mantissa = bits & 0xFFFFFFFFFFFFF;
    cout << sign << " " << bitset<11>(exponent) << " " << bitset<52>(mantissa) << endl;
}

class FP64 {
    ull frac : 52;
    ull exp : 11;
    ull sign : 1;

    void ConvertToDouble1(double* D) noexcept {
        ull t = 0;
        t |= (((ull)this->sign) << 63);
        t |= (((ull)this->exp) << 52);
        t |= (this->frac);
        *(D) = *((double*)(&t));
    }
    double ConvertToDouble2() noexcept {
        double d = 0;
        ull t = 0;
        t |= (((ull)this->sign) << 63);
        t |= (((ull)this->exp) << 52);
        t |= (this->frac);
        d = *((double*)(&t));
        return d;
    }
    FP64 ConvertToFP64(double A) noexcept {
        FP64 res;
        ull bits = *reinterpret_cast<ull*>(&A);
        res.sign = (bits >> 63) & 1;
        res.exp = (bits >> 52) & 0x7FF;
        res.frac = bits & 0xFFFFFFFFFFFFF;
        return res;
    }
    void to_bin(ull n) {
        if (n > 1)
            to_bin(n >> 1);
        cout << (char)((n & 1) + '0');
    }

    FP64 Mul(FP64 A, FP64 B) noexcept {
        FP64 Res;
        unsigned long long shiftExp = (1 << (11 - 1)) - 1;
        //unsigned long long temp = (N1.in.frac / (1ULL << 26)) * (N2.in.frac / (1ULL << 26)) + N1.in.frac + N2.in.frac + (1ULL << 52);
        Res.exp = ((A.exp - shiftExp) + (B.exp - shiftExp) + shiftExp);
        ull m1strih1 = ((A.frac & 0xffffffff00000000) >> 32ULL);
        ull m1strih2 = (A.frac & 0x00000000ffffffff);
        ull m2strih1 = ((B.frac & 0xffffffff00000000) >> 32ULL);
        ull m2strih2 = (B.frac & 0x00000000ffffffff);
        ull mult = m2strih2 * m1strih2;
        mult /= (1ULL << 32ULL);
        ull temp2 = m1strih1 * m2strih2 + m2strih1 * m1strih2 + mult;
        //tobin(temp2);
        //printf("\n");
        ull temp = (1ULL << 12) * m1strih1 * m2strih1 + ((temp2) / (1ULL << 20ULL)) + A.frac + B.frac + (1ULL << 52); // +2^-20
        //приведение экспоненты
        //unsigned long long flag = temp2 & (1ULL << 20);
        if (temp >= (1ULL << 53)) { // меняем значение сдвига, изначально было 53
          //flag += (temp & 1ULL);
            temp = (temp >> 1ULL);
            Res.exp++;
        }
        // temp2 = 1 + 2^-21
        //округление по отрезанным битам
        ull flag1 = temp2 & (1ULL << 19);
        ull flag2 = temp & 1ULL;
        if (flag2 && flag1)
            temp++;
        if (flag2)
            if (!flag1)
                temp += 0;
        if (!flag2) {
            // 20 битов начиная с 45 бита по 64, если это (2 ^ 20 - 1)
            ull t = ((1ULL << 20) - 1);
            ull tmp = temp2 & t;
            if (tmp > (1ULL << 19))
                temp++;
            // unsigned long long tmp = temp2 - (temp2 >> 20);
             //if (tmp == (1ULL << ))
        }
        /*if (flag) {
          //temp += temp & ((1ULL << 27) - 1);
          temp++;
          printf("hi\n");
        }*/
        Res.frac = temp;
        Res.sign = A.sign + B.sign;
        return Res;
    }
public:
    FP64() noexcept {
        frac = exp = sign = 0;
    }
    FP64(double A) noexcept {
        ull bits = *reinterpret_cast<ull*>(&A);
        this->sign = (bits >> 63) & 1;
        this->exp = (bits >> 52) & 0x7FF;
        this->frac = bits & 0xFFFFFFFFFFFFF;
    }
    FP64(const FP64& A) noexcept {
        this->frac = A.frac;
        this->exp = A.exp;
        this->sign = A.sign;
    }
    FP64& operator=(const FP64& A) {
        this->frac = A.frac;
        this->exp = A.exp;
        this->sign = A.sign;
        return (*this);
    }
    FP64& operator=(double d) {
        (*this) = ConvertToFP64(d);
        return (*this);
    }
    FP64& operator+=(const FP64& A) {
        if (!sign) {
            if (!A.sign)
                (*this) = AddFP64((*this), A);
            else (*this) = SubFP64((*this), A);
        }
        else {
            if (A.sign) {
                (*this) = AddFP64((*this), A);
                (*this).sign = 1;
            }
            else (*this) = SubFP64(A, (*this));
        }
        return (*this);
    }
    FP64 operator+(const FP64& A) noexcept {
        FP64 temp = (*this);
        temp += A;
        return temp;
    }
    FP64& operator-=(const FP64& A) {
        if (sign) {
            if (!A.sign) {
                (*this) = AddFP64((*this), A);
                (*this).sign = 1;
            }
            else (*this) = SubFP64(A, (*this));
        }
        else {
            if (A.sign)
                (*this) = AddFP64((*this), A);
            else (*this) = SubFP64((*this), A);
        }
        return (*this);
    }
    FP64 operator-(const FP64& A) const noexcept {
        FP64 temp = (*this);
        temp -= A;
        return temp;
    }
    FP64& operator*=(const FP64& A) {
        //(*this) = Mul(*(this), A);
        (*this) = MulFP64((*this), A);
        return (*this);
    }
    FP64 operator*(const FP64& A) const noexcept {
        FP64 temp = (*this);
        temp *= A;
        return temp;
    }
    FP64& operator/=(const FP64& A) {
        (*this) = DivFP64((*this), A);
        return (*this);
    }
    FP64 operator/(const FP64& A) /*noexcept?*/ {
        FP64 temp = (*this);
        temp /= A;
        return temp;
    }
    operator double() noexcept {
        double d = this->ConvertToDouble2();
        return d;
    }
    friend std::ostream& operator << (std::ostream& stream, FP64& A) {
        double out = A;
        stream << out;
        return stream;
    }
    friend std::istream& operator >> (std::istream& stream, FP64& A) {
        double in = A;
        stream >> A;
        A = in; //
        return stream;
    }
    void bits() {
        std::cout << sign << " ";
        to_bin(exp);
        std::cout << " ";
        to_bin(frac);
        std::cout << " ";
        std::cout << '\n';
    }
    bool is_subnorm() const noexcept {
        return (exp == 0) && (frac != 0);
    }
    bool is_inf() const noexcept {
        return (exp == 0x7FF) && (frac == 0);
    }
    bool is_NaN() const noexcept {
        return (exp == 0x7FF) && (frac != 0);
    }
    bool is_null() const noexcept {
        return (exp == 0) && (frac == 0);
    }
    bool operator == (const FP64& A) const noexcept {
        if (!(is_NaN() || A.is_NaN())) {
            if (frac == A.frac && exp == A.exp && sign == A.sign) //
                return true;
            else return false;
        }
        else return false;
    }
    // переделать
    bool operator < (FP64& A) noexcept {
        if (!(is_NaN() || A.is_NaN())) {
            FP64 tmp = *this - A;
            return (tmp.GetSign());
        }
        else return false;
    }
    bool operator > (FP64& A) noexcept {
        if (!(is_NaN() || A.is_NaN())) {
            FP64 tmp = A - *this;
            return (tmp.GetSign());
        }
        else return false;
    }
    //
    bool operator != (const FP64& A) const noexcept {
        if (!(is_NaN() || A.is_NaN())) {
            return (!(*this == A));
        }
        else return true; // наверно так
    }
    void SetSign(int s) {
        if (s == 0 || s == 1)
            sign = s;
        else throw std::logic_error("bad sign");
    }
    void SetExp(int e) {
        exp = e; // нужна ли проверка e?
    }
    void SetFrac(ull f) {
        frac = f;
    }
    ull GetSign() {
        return sign;
    }
    ull GetExp() {
        return exp;
    }
    ull GetFrac() {
        return frac;
    }
    FP64 MulFP64(FP64 N1, FP64 N2) {
        FP64 Res;
        unsigned long long shiftExp = 1023;
        //long long tempexp = ((ll(N1.exp) - shiftExp) + (ll(N2.exp) - shiftExp) + shiftExp);
        ll tempexp = ((ll)N1.exp - shiftExp) + (ll)N2.exp;
        //unsigned long long temp = (N1.in.frac / (1ULL << 26)) * (N2.in.frac / (1ULL << 26)) + N1.in.frac + N2.in.frac + (1ULL << 52);
        //Res.exp = ((N1.exp - shiftExp) + (N2.exp - shiftExp) + shiftExp);
        unsigned long long m1strih1 = ((N1.frac & 0xffffffff00000000) >> 32ULL);
        unsigned long long m1strih2 = (N1.frac & 0x00000000ffffffff);
        unsigned long long m2strih1 = ((N2.frac & 0xffffffff00000000) >> 32ULL);
        unsigned long long m2strih2 = (N2.frac & 0x00000000ffffffff);
        unsigned long long mult = m2strih2 * m1strih2;
        //  to_bin(mult);
        //   printf("\n");
        // попробуем сделать такое же округление для mult
        //    ull fflag1 = mult & (1ULL << 31);
        ull charge = 0;
        charge += (mult & 0x00000000ffffffff);
        mult /= (1ULL << 32ULL);
        //     ull fflag2 = mult & 1ULL;
        //
        unsigned long long temp2 = m1strih1 * m2strih2 + m2strih1 * m1strih2 + mult;
        charge += ((temp2 & 0x00000000000fffff) << 32); // необходимо добавить последние 20 битов числа temp2
        //
        /*     if (!fflag1)
         ;
         if (fflag1) {
         ull t = ((1ULL << 32) - 1);
         ull tmp = mult & t;
         if (tmp > (1ULL << 31))
         temp2++;
         else {
         if (fflag2)
         temp2++;
         }
         }*/
         //
         //to_bin(temp2);
        //   printf("\n");
        unsigned long long temp = (1ULL << 12) * m1strih1 * m2strih1 + ((temp2) / (1ULL << 20ULL)) + N1.frac + N2.frac + (1ULL << 52); // +2^-20
        //приведение экспоненты
        //unsigned long long flag = temp2 & (1ULL << 20);
        //ull tmp = 0;//
        ull flag1 = 0; //
        bool f = false;
        if (temp >= (1ULL << 53)) {
            ///  flag1 += (temp & 1ULL);
            // tmp+=(temp & 1ULL);//
            //  int a = temp & 1ULL;
            charge += ((temp & 1ULL) << 52);
            temp = (temp >> 1ULL);
            //      temp += a;
            //Res.exp++;
            tempexp++;
            f = true;
        }
        if (temp == ((1ULL << 53) - 1)) { // считаю, что f == false
            flag1 += charge & (1ULL << 51);
            unsigned long long flag2 = temp & 1ULL;
            if (flag1) {
                ull t;
                t = ((1ULL << 51) - 1);
                ull tmp = charge & t;
                if (tmp)
                    temp++;
                else {
                    if (flag2)
                        temp++;
                }
            }
            if (temp >= (1ULL << 53)) {
                charge += ((temp & 1ULL) << 52);
                temp = (temp >> 1ULL);
                tempexp++;
            }
            ull fordebug = temp - (1ULL << 52);
            Res.exp = (ull)tempexp;
            Res.frac = fordebug;
            Res.sign = N1.sign + N2.sign;
            return Res;
        }
        if (tempexp < 0) {
            Res.exp = 0;
            Res.frac = 0;
            return Res;
        }
        if (tempexp >= 2047 || N1.is_inf() || N2.is_inf() || N1.is_NaN() || N2.is_NaN()) {
            if (N1.is_inf() && N2.is_inf()) {
                if (N1.sign != N2.sign) {
                    Res.exp = 2047;
                    Res.frac = 1;
                    return Res;
                }
                else {
                    return N1;
                }
            }
            if (N1.is_inf()) 
                return N1;
            if (N2.is_inf()) 
                return N2;
            if (N1.is_NaN()) 
                return N1;
            if (N2.is_NaN()) 
                return N2;
            if (tempexp >= 2047) {
                Res.exp = 2047;
                Res.frac = 0;
                return Res;
            }
        }
        // temp2 = 1 + 2^-21
        //округление по отрезанным битам
        // flag1 += temp2 & (1ULL << 19); // 45-ый бит
        if (f) // если f, то число 53 битное, если ~f, то число 52 битное
            flag1 += charge & (1ULL << 52);
        else flag1 += charge & (1ULL << 51);
        //  unsigned long long flag2 = (temp2 << 20ULL) & 1ULL; // менял, 44-ый бит
        unsigned long long flag2 = temp & 1ULL; //- изначально было
        /*      if (flag2 && flag1)
         temp++;
         if (flag2)
         if (!flag1)
         temp += 0;
         if (!flag2) {
         // 20 битов начиная с 45 бита по 64, если это (2 ^ 20 - 1)
         unsigned long long t = ((1ULL << 20) - 1);
         unsigned long long tmp = temp2 & t;
         if (tmp > (1ULL << 19))
         temp++;
         // unsigned long long tmp = temp2 - (temp2 >> 20);
         //if (tmp == (1ULL << ))
         }*/

        if (!flag1)
            ;
        if (flag1) {
            ull t;
            if (f)
                t = ((1ULL << 52) - 1);
            else t = ((1ULL << 51) - 1);
            ull tmp = charge & t;
            if (f) {
                //if (tmp > (1ULL << 52)){
                if (tmp) {
                    temp++;
                }
                else {
                    if (flag2)
                        temp++;
                }
            }
            else {
                // if (tmp > (1ULL << 51))
                if (tmp)
                    temp++;
                else {
                    if (flag2)
                        temp++;
                }
            }
        }
        // temp++;
        ull fordebug = temp - (1ULL << 52);
        Res.exp = (ull)tempexp;
        Res.frac = fordebug;
        Res.sign = N1.sign + N2.sign;
        return Res;
    }
    //   FP64 AddFP64(FP64 A, FP64 B) {
    //       FP64 Res;
    //   //Res->expect = A->expect + B->expect;
    //   long long round;
    //   long long diff = abs(ll(A.exp) - ll(B.exp));
    //   long long MaxExp, MaxSign, MinSign;
    //   /*if (A.exp >= B.exp) {
    //       MaxExp = A.exp;
    //       MaxSign = A.sign;
    //       MinSign = B.sign;
    //       round = (1ULL << 52) + B.frac;
    //   }
    //   else {
    //       MaxExp = B.exp;
    //       MaxSign = B.sign;
    //       MinSign = A.sign;
    //       round = (1ULL << 52) + A.frac;
    //   }*/ // для начала проверим, если экспонента первого не меньше экспоненты второго
    //   MaxExp = A.exp;
    //   MaxSign = A.sign;
    //   MinSign = B.sign;
    //   round = (1ULL << 52) + B.frac;
    //   long long t;
    //   long long sh;
    //   if (53 - diff < 0)
    //       sh = 0;
    //   else { // 53 вместо 52, при 52 могут скзаться на округлении
    //       if (53 - diff == 0)
    //           sh = 1;
    //       else sh = (1LL << (52 - diff)); // LL после 1
    //   }
    //   ull check1, check2;
    //   if (MaxExp - A.exp > 63)
    //       check1 = 0;
    //   else check1 = (A.frac >> (MaxExp - A.exp));
    //   if (diff > 63)
    //       check2 = 0;
    //   else check2 = (round >> diff);
    //   //t = -1LL*(2 * A.sign - 1) * (A.frac >> (MaxExp - A.exp)) - (2 * B.sign - 1) * (B.frac >> (MaxExp - B.exp)) - (2 * MinSign - 1) * sh - (2 * MaxSign - 1) * (1LL << 52LL);
    ////if (A->in.exp >= B->in.exp) {
    //   //   // t = -(2 * MaxSign - 1) * (1LL << 52LL) - 1LL*(2 * A->in.sign - 1) * (A->in.frac >> (MaxExp - //A->in.exp)) - (2 * B->in.sign - 1) * (round >> diff);
    //    // t = (1LL << 52LL) + (A->in.frac >> (MaxExp - A->in.exp)) + (round >> diff);
    // //  }
    //// else t = (1LL << 52LL) + (B->in.frac >> (MaxExp - B->in.exp)) + (round >> diff);
    //  // t = (1LL << (A.exp - MaxExp + 52)) + это строку нужно оставить, вопрос - что добавлять?
    //  // t = -(2 * MaxSign - 1) * (1LL << 52LL) - 1LL * (2 * A.sign - 1) * (A.frac >> (MaxExp - A.exp)) - (2 * B.sign - 1) * (round >> diff);
    //   t = (1LL << 52LL) + check1 + check2;
    //   //t = A.frac + (1LL << 52) + B.frac * (1LL >> diff) + (1LL << (52 - diff));// тоже самое записано, что строчкой выше
    //   if (t < 0) {
    //    t = -t;
    //       Res.sign = 1;
    //   }
    //   else Res.sign = 0;
    //   //if (diff <= 63) { // 63 или 64, пока не понятно
    //   //    long long flag1;
    //   //    //if (diff - 1 > 63)
    //   //      //  flag1 = 0;
    //   //    //else
    //   //    flag1 = round & (1ULL << (diff - 1));
    //   //    long long flag2 = t & 1ULL;
    //   //    if (flag1 && flag2)
    //   //        t++;
    //   //    if (flag2)
    //   //        if (!flag1)
    //   //            t += 0;
    //   //    if (!flag2) {
    //   //        long long tmp = (1ULL << diff) - 1;
    //   //        long long temp = round & tmp;
    //   //        //if (temp > tmp) - здесь наверно ошибка
    //   //        //if (temp >)
    //   //            t++;
    //   //    }
    //   //}
    //   if (diff <= 104) { // 63 или 64, пока не понятно
    //       long long flag1;
    //       //if (diff - 1 > 63)
    //         //  flag1 = 0;
    //       //else
    //       flag1 = round & (1ULL << (diff - 1)); //поменял
    //       //flag1 = round & (1ULL << diff);
    //       long long flag2 = t & 1ULL;
    //       if (!flag1)
    //           ;
    //       if (flag1) {
    //           ull tt;
    //           tt = (1ULL << diff) - 1; // значение свдигов должно совпадать diff и diff - 1
    //           ull tmp = round & tt;
    //           if (tmp)
    //               t++;
    //           else {
    //               if (flag2)
    //                   t++;
    //           }
    //       }
    //
    //       //if (flag1 && flag2)
    //       //    t++;
    //       //if (flag2)
    //       //    if (!flag1)
    //       //        t += 0;
    //       //if (!flag2) {
    //       //    long long tmp = (1ULL << diff) - 1;
    //       //    long long temp = round & tmp;
    //       //    //if (temp > tmp) - здесь наверно ошибка
    //       //    //if (temp >)
    //       //        t++;
    //       //}
    //   }
    //   while (t < (1LL << 52LL) && MaxExp > 0) { // получить коэф умножения, через деление
    //       t *= 2;
    //       MaxExp--;
    //   }
    //   if (t >= (1LL << 53) && (MaxExp < (1LL << 11LL))) {
    //       t = t >> 1LL;
    //       MaxExp++;
    //   }
    //   Res.exp = MaxExp;
    //   Res.frac = (t - (1LL << 52LL));
    //   return Res;
    //   }
    FP64 AddFP64(FP64 A, FP64 B) {
        if (A.is_null())
            return B;
        if (B.is_null())
            return A;
        FP64 Res;
        ull MaxExp = A.exp, MinExp = B.exp, MaxFrac = A.frac, MinFrac = B.frac, MaxSign = A.sign, MinSign = B.sign;
        if (B.exp > A.exp) {
            MaxExp = B.exp;
            MinExp = A.exp;
            MaxFrac = B.frac;
            MinFrac = A.frac;
            MaxSign = B.sign;
            MinSign = A.sign;
        }
        Res.exp = MaxExp;
        ull diff = MaxExp - MinExp;
        if (diff >= 54) {
            Res.frac = MaxFrac;
            Res.sign = MaxSign;
            return Res;
        }
        ull round = ((1ULL << 52) + MinFrac) & ((1ULL << diff) - 1);
        //ull t = (1ULL << 52) + A.frac + (((1ULL << 52) + B.frac) >> (diff));
        long long sh;
        if (53 - diff <= 0) // добавил =
            sh = 0;
        else { // 53 вместо 52, при 52 могут скзаться на округлении
            if (53 - diff == 0)
                sh = 1;
            else sh = (1LL << (52 - diff)); // LL после 1
        }
        //ull deb = B.frac >> diff;
        ull t = MaxFrac + ((MinFrac >> diff) + sh) + (1ULL << 52);
        //
        if (A.is_subnorm() || B.is_subnorm()) {
            if (A.is_subnorm() && B.is_subnorm())
                t = MaxFrac + MinFrac;
            else {
                if (A.is_subnorm()) {
                    diff = B.exp - 1;
                    t = (1LL << 52) + B.frac + (A.frac >> diff);
                    round = A.frac & ((1LL << diff) - 1);
                    Res.sign = B.sign;
                    Res.exp = B.exp;
                }
                else {
                    diff = A.exp - 1;
                    t = (1LL << 52) + A.frac + (B.frac >> diff);
                    round = B.frac & ((1LL << diff) - 1);
                    Res.sign = A.sign;
                    Res.exp = A.exp;
                }
            }
        }
        //
        bool f = false;
        if ((t >= (1ULL << 53)) && (Res.exp > 0)) { //
            round += ((t & 1ULL) << diff);
            t = t >> 1LL;
            Res.exp++;
            f = true;
        }
        ull flag1, flag2;
        flag2 = t & 1ULL;
        if (f)
            flag1 = round & (1ULL << diff);
        else flag1 = round & (1ULL << (diff - 1));
        if (!flag1)
            ;
        if (flag1) {
            ull temp;
            if (f)
                temp = (1ULL << diff) - 1;
            else temp = (1ULL << (diff - 1)) - 1;
            ull tmp = round & temp;
            if (tmp) {
                t++;
            }
            else {
                if (flag2)
                    t++;
            }
        }
        while (t < (1LL << 52LL) && Res.exp > 1) { 
            t = t << 1LL;
            Res.exp--;
        }
        //
        if (t >= (1LL << 52) && (Res.exp == 0)) {
            Res.exp = 1;
            t -= (1LL << 52);
        }
        //
        if (Res.exp >= 2047 || A.exp == 2047 || B.exp == 2047) {
            if (A.is_inf() && B.is_inf()) {
                if (A.sign != B.sign) {
                    Res.exp = 2047;
                    Res.frac = 1;
                    return Res;
                }
                else 
                    return A;
            }
            if (A.is_inf()) 
                return A;
            if (B.is_inf()) 
                return B;
            if (A.is_NaN())
                return A;
            if (B.is_NaN()) 
                return B;
            if (Res.exp >= 2047) {
                Res.exp = 2047;
                Res.frac = 0;
                return Res;
            }
        }
        //
        ull debug = t - (1ULL << 52);
        Res.frac = debug;
        return Res;
    }

    FP64 SubFP64(FP64 A, FP64 B) {
        if (B.is_null())
            return A;
        FP64 Res;
        if (A.is_null()) {
            Res = B;
            Res.sign = 1;
            return Res;
        }
        FP64 newA = A, newB = B;
        bool sw = false;
        if (B.exp > A.exp) {
            FP64 tmp = newA;
            newA = newB;
            newB = tmp;
            sw = true;
        }
        Res.exp = newA.exp;
        ull diff = newA.exp - newB.exp;
        if (diff >= 55) {
            Res.frac = newA.frac;
            if (sw) {
                Res.sign++;
            }
            return Res;
        }
        long long sh;
        ull round, t;
        if (newA.exp == newB.exp) {
            sh = 1ULL << 52;
            if (newA.frac >= newB.frac) {
                t = newA.frac - newB.frac;
                Res.sign = 0;
            }
            else {
                t = newB.frac - newA.frac;
                Res.sign = 1;
            }
        }
        else {
            if (diff > 52)
                sh = 0;
            else sh = (1ULL << (52 - diff));
            t = newA.frac - (newB.frac >> diff) - sh + (1ULL << 52);
            bool f = false;
        }
        round = ((1ULL << 52) + newB.frac) & ((1ULL << diff) - 1);
        if (diff == 0)
            round = 0;
        //
        if (newA.is_subnorm() || newB.is_subnorm()) {
            if (newA.is_subnorm() && newB.is_subnorm()) {
                if (newA.frac >= newB.frac) {
                    t = newA.frac - newB.frac;
                    Res.sign = newA.sign;
                }
                else {
                    t = newB.frac - newA.frac;
                    Res.sign = newB.sign;
                }
            }
            else {
                if (newA.is_subnorm()) {
                    diff = B.exp - 1;
                    t = (1LL << 52) + newB.frac - (newA.frac >> diff);
                    round = newA.frac & ((1LL << diff) - 1);
                    Res.sign = newB.sign;
                    Res.exp = newB.exp;
                }
                else {
                    diff = newA.exp - 1;
                    t = (1LL << 52) + newA.frac - (newB.frac >> diff);
                    round = newB.frac & ((1LL << diff) - 1);
                    Res.sign = newA.sign;
                    Res.exp = newA.exp;
                }
            }
        }
        //
        if (round) {
            t--;
            round = (1ULL << diff) - round;
        }
        while (t < (1LL << 52LL) && Res.exp > 1) {
            t = t << 1ULL;
            Res.exp--;
            t += (round >> (diff - 1)) & 1ULL;
            round = round & ((1ULL << (diff - 1)) - 1);
            diff--;
        }
        if (!t) {
            Res = FP64();
            return Res;
        }
        if ((t >= (1ULL << 53)) && Res.exp > 0) {
            round += ((t & 1ULL) << diff);
            t = t >> 1LL;
            Res.exp++;
            diff++;
        }
        if ((round >> (diff - 1)) & 1ULL) {
            round = round & ((1ULL << (diff - 1)) - 1);
            if (round > 0) {
                t++;
            }
            else {
                t += (t & 1ULL);
            }
            if ((t >= (1ULL << 53)) && (Res.exp >= 1)) {
                t = t >> 1ULL;
                Res.exp++;
            }
        }

        if (t < (1ULL << 52) && Res.exp == 1) {
            Res.exp = 0;
        }

        if (Res.exp >= 2047 || A.exp == 2047 || B.exp == 2047) {
            if (newA.is_inf() && newB.is_inf()) {
                if (newA.sign == newB.sign) {
                    Res.exp = 2047;
                    Res.frac = 1;
                    return Res;
                }
                else
                    return A;
            }
            if (newA.is_inf())
                return newA;
            if (newB.is_inf())
                return newB;
            if (newA.is_NaN())
                return newA;
            if (newB.is_NaN())
                return newB;
            if (Res.exp >= 2047) {
                Res.exp = 2047;
                Res.frac = 0;
                return Res;
            }
        }

        ull debug = t - (1ULL << 52);
        Res.frac = debug;
        if (sw)
            Res.sign++;
        return Res;
    }

    //FP64 fmaFP64(FP64 A, FP64 B, FP64 C) {
    //    FP64 Res;
    //    ull m1 = A.frac + (1ULL << 52);
    //    ull m2 = B.frac + (1ULL << 52);
    //    ull m3 = C.frac + (1ULL << 52);
    //    ull MMSb = 0xFFFFFFFF00000000; // маска на старшие 32 бита
    //    ull MLSb = 0xFFFFFFFF; // маска на младшие 32 бита
    //    ull p = 53;
    //    ll ec = C.exp;
    //    ll eab = A.exp + B.exp - 1023;
    //    //ll eab = A.exp + B.exp - 1023 + 2;
    //    ll na = 1, nb = 1, nc = 1;
    //    if (A.is_subnorm()) {
    //        //na = (12 - __lzcnt64(A.GetFrac()));
    //        na = (12 - __builtin_clzll(A.GetFrac()));
    //        m1 -= (1ULL << 52);
    //    }
    //    if (B.is_subnorm()) {
    //        //nb = (12 - __lzcnt64(B.GetFrac()));
    //        nb = (12 - __builtin_clzll(B.GetFrac()));
    //        m2 -= (1ULL << 52);
    //    }
    //    if (C.is_subnorm()) {
    //        //nc = (12 - __lzcnt64(C.GetFrac()));
    //        nc = (12 - __builtin_clzll(C.GetFrac()));
    //        m3 -= (1ULL << 52);
    //    }
    //    ll d = ec - eab - 1 - nc + na + nb;
    //    ull a1, a2, a3, a4, a5, a6;
    //    ull b1, b2, b3, b4, b5, b6;
    //    ull ab1 = 0, ab2 = 0, ab3 = 0, ab4 = 0, ab5 = 0, ab6 = 0;
    //    // число хранится так: ab3::ab4::ab5::ab6::ab1::ab2; ab3 - старший
    //    a1 = a2 = a3 = a4 = 0;
    //    a6 = (m1 & 0x00000000ffffffff);
    //    a5 = (m1 >> 32ULL);
    //    b1 = b2 = b3 = b4 = 0;
    //    b6 = (m2 & 0x00000000ffffffff);
    //    b5 = (m2 >> 32ULL);
    //    //
    //    ab6 = a6 * b6;
    //    ab5 += (ab6 >> 32);
    //    ab6 = (ab6 & MLSb);
    //    // заполнили ab6
    //    ////ab5 += (a5 * b5); // всё хорошо, так как изначально a5 и b5 - 20-битные, переполнения не должно быть
    //    ab5 += a5 * b6;
    //    ab5 += a6 * b5;
    //    ab4 += (ab5 >> 32);
    //    ab5 = (ab5 & MLSb);
    //    // заполнили ab5
    //    ab4 += a5 * b5;
    //    ab3 += (ab4 >> 32);
    //    ab4 = (ab4 & MLSb);
    //    // заполнили ab4 и ab3
    //    ull razrab; // не может быть 0, если она 0, то фактически a * b == 0;
    //    if (ab3)
    //        razrab = 4;
    //    else {
    //        if (ab4)
    //            razrab = 3;
    //        else {
    //            if (ab5)
    //                razrab = 2;
    //            else razrab = 1;
    //        }
    //    }
    //    // сдвигаем на p+3 = 56 разрядов вправо

    //    ull orig_ab3 = ab3;
    //    ull orig_ab4 = ab4;
    //    ull orig_ab5 = ab5;
    //    ull orig_ab6 = ab6;
    //    ull orig_ab1 = ab1;
    //    ull orig_ab2 = ab2;
    //    // сначало сдвинем на 32 бита вправо, потом ещё на 24
    //    ab3 = 0;
    //    ab4 = orig_ab3 & MLSb;
    //    ab5 = orig_ab4 & MLSb;
    //    ab6 = orig_ab5 & MLSb;
    //    ab1 = orig_ab6 & MLSb;  // ab2 = 0
    //    //
    //    orig_ab3 = ab3;
    //    orig_ab4 = ab4;
    //    orig_ab5 = ab5;
    //    orig_ab6 = ab6;
    //    orig_ab1 = ab1;
    //    orig_ab2 = ab2;
    //    // 0xFFFFFF - маска на 24 единицы
    //    ab4 = ab4 >> 24ULL;
    //    orig_ab5 = ab5 + ((orig_ab4 & 0xFFFFFF) << 32ULL);
    //    ab5 = orig_ab5 >> 24ULL;
    //    orig_ab6 = ab6 + ((orig_ab5 & 0xFFFFFF) << 32ULL);
    //    ab6 = orig_ab6 >> 24ULL;
    //    orig_ab1 = ab1 + ((orig_ab6 & 0xFFFFFF) << 32ULL);
    //    ab1 = orig_ab1 >> 24ULL;
    //    orig_ab2 = (orig_ab1 & 0xFFFFFF) << 32ULL;
    //    ab2 = orig_ab2 >> 24ULL; // на самом деле ab2 должен оставаться 0
    //    // завершили сдвиг числа на 56 битов вправо

    //    // два нуля сзади добавляются автоматически (они уже есть)

    //    ll dstrih;
    //    if (d <= -105) {
    //        dstrih = 3 * 53 + 4;
    //        Res.exp = eab;
    //    }
    //    if (d > -105 && d <= 2) {
    //        dstrih = 56 - d;
    //        Res.exp = eab;
    //    }
    //    if (d > 2 && d <= 55) {
    //        dstrih = 56 - d;
    //        Res.exp = ec;
    //    }
    //    if (d >= 56) {
    //        dstrih = 0;
    //        Res.exp = ec;
    //    }

    //    // сдвигаем c вправо на dstrih бит

    //    ull number = dstrih / 32;
    //    ull position = dstrih % 32;
    //    // выгодно было бы превратить представление 192-битного числа в массив из 6 чисел для избежания if

    //    ull c[6] = { 0 };
    //    //c0::c1::c2::c3::c4::c5 - хранится
    //    ull available_bits = 32 - position;
    //    ull temp, round = 0;

    //    /*if (dstrih != 163 && !C.is_subnorm()) {
    //        c[number] = (((~0ULL) << (53 - available_bits)) & C.frac); // старшие available_bits мантисы C
    //        c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & C.frac) & (((~0ULL) << (53 - available_bits - 32)));
    //        if (53 - available_bits - 32 > 0)
    //            c[number + 2] = ((1ULL << (53 - available_bits - 32)) - 1) & C.frac;
    //    }
    //    else {
    //        round = C.frac;
    //    }*/
    //    // можно без if написать
    //    ull rfc = 64 - __builtin_clzll(m3); // разрядность мантисы C до разложения на 32-битные числа
    //    if (dstrih != 163) {
    //        if (!C.is_subnorm()) {
    //            c[number] = (((~0ULL) << (53 - available_bits)) & m3); // старшие available_bits мантисы C
    //            if (53 - available_bits - 32 >= 0) {// =?
    //                c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & m3) & (((~0ULL) << (53 - available_bits - 32)));
    //                c[number + 2] = ((1ULL << (53 - available_bits - 32)) - 1) & m3;
    //            }
    //            else {
    //                c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & m3) & m3;
    //            }
    //        }
    //        else {
    //            if (rfc >= available_bits) {
    //                c[number] = (((~0ULL) << (rfc - available_bits)) & m3); // старшие available_bits мантисы C
    //                if (rfc - available_bits >= 32) {
    //                    c[number + 1] = (((1ULL << (rfc - available_bits)) - 1) & m3) & (((~0ULL) << (rfc - available_bits - 32)));
    //                    c[number + 2] = ((1ULL << (rfc - available_bits - 32)) - 1) & m3;
    //                }
    //                else {
    //                    c[number + 1] = (((1ULL << (rfc - available_bits)) - 1) & m3) & m3;
    //                }
    //            }
    //            else {
    //                c[number] = ((1ULL << rfc) - 1) & m3; // не уверен, что идет в ячейку number
    //            }
    //        }
    //    }
    //    else { // если dstrih == 163
    //        if (!C.is_subnorm()) {
    //            round += (m3 & ((1ULL << 52) - 1));
    //            // последние 52бита уходят в round (биты округления)
    //            c[5] = 1; // так как мантиса нормального числа представляется из себя выражение 2^52 + x, то старший бит всегда будет равен 1
    //        }
    //        else {
    //            round += (m3 & ((1ULL << (rfc - 1)) - 1));
    //            // c полностью уйдет в 0
    //        }
    //    }
    //    // для начала разберемся со сложением
    //    ull s[6] = { 0 }; // сумма (a*b) + c
    //    //// лидирующий ненулевой бит c должен находиться в числе c[number] на позиции position
    //    s[5] = ab2 + c[5];
    //    s[4] = ab1 + c[4] + ((s[5] & MMSb) >> 32);
    //    s[3] = ab6 + c[3] + ((s[4] & MMSb) >> 32);
    //    s[2] = ab5 + c[2] + ((s[3] & MMSb) >> 32);
    //    s[1] = ab4 + c[1] + ((s[2] & MMSb) >> 32);
    //    s[0] = ab3 + c[0] + ((s[1] & MMSb) >> 32);
    //    // теперь необходимо сдвинуть сумму s влево так, чтобы лидирующий бит был 1, найдем величину этого сдвига
    //    ll dstrih2;
    //    if (d <= 2) {
    //        ll l; // кол-во ведущих нулей в 109 последних битах числа s
    //        // получается надо начинать смотреть с s[2]
    //        if (s[2] != 0 && (s[2] & ((1ULL << 13) - 1))) {
    //            // можем смотреть в первых 13 битах
    //            ll temp = s[2] & ((1ULL << 13) - 1);
    //            l = 64 - __builtin_clzll(temp);
    //        }
    //        else {
    //            if (s[3] != 0)
    //                l = 64 - __builtin_clzll(s[3]);
    //            else {
    //                if (s[4] != 0)
    //                    l = 64 - __builtin_clzll(s[4]);
    //                else {
    //                    if (s[5] != 0)
    //                        l = 64 - __builtin_clzll(s[5]);
    //                    // иначе число c == 0
    //                }
    //            }
    //        }
    //        //рассчитали l
    //        if (eab - l + 2 >= 1) { // в таком случае результат fma будет нормальным числом
    //            dstrih2 = 55 + l;
    //            Res.exp = eab - l + 2;
    //        }
    //        else { // иначе результат будет сабнормалом
    //            Res.exp = 1;
    //            ll leadbits = 1 - (eab - l + 2); // столько лидирующих нулей получит мантиса результата
    //            dstrih2 = 57 - 1 + eab;
    //        }
    //    }
    //    else { // d > 2
    //        dstrih2 = dstrih;
    //        // в таком случае сдвиг влево отменяет изначальный сдвиг вправо, но после сдвига ведущая 1 может сменить позицию: при слоежении оказаться левее на 1 бит, при вычитании оказаться правее на 1 бит?
    //    }
    //    // я так понимаю после сложения a*b и c необходимо сначала свдинуть результать на dstrih2 влево, затем округлить взяв старшие 53бита, остальные добавить в round (они буддут составлять 110бит), сдвинув их при это на разрядность round уже имеющегося. Затем вычесть из мантисы мнимую единицу, то есть (1ULL << 52) и получившуюся мантису округлить по round.
    //    // наверно можно не двигать всё число, а найти лидирующий бит и свдинуть его и 52 бита за ним, а так же round
    //    number = dstrih2 / 32;
    //    position = dstrih2 % 32;
    //    // сначало сдвигаем цельными фрагментами
    //    ull temps[6] = { 0 };
    //    for (int i = 0; i < 6; ++i) {
    //        if (i + number < 6)
    //            temps[i] = (s[i + number] & 0xFFFFFFFF);
    //    }
    //    memcpy(s, temps, sizeof(temps));
    //    if (position) {
    //        ull carry = 0;
    //        for (int i = 5; i >= 0; --i) {
    //            ull current = (s[i] & 0xFFFFFFFF); // на всякий случай
    //            ull shifted = (current << position) + carry;
    //            carry = (current << position) >> 32ULL;
    //            s[i] = (shifted & 0xFFFFFFFF);
    //        }
    //    }
    //    // поиск лидирующей 1
    //    ull newfrac = 0;
    //    for (int i = 0; i < 6; i++) {
    //        if (s[i]) {
    //            ull curnumb = i;
    //            ull cnt = 0;
    //            while (cnt < 53 && curnumb < 6) {
    //                available_bits = 64 - __builtin_clzll(s[curnumb]);
    //                if (available_bits + cnt > 53)
    //                    available_bits = 53 - cnt;
    //                newfrac += (s[curnumb] & ((1ULL << available_bits) - 1));
    //                curnumb++;
    //                cnt += available_bits;
    //            }
    //        }
    //    }
    //    // проверка (пока без округления)
    //    Res.sign = 0;
    //    Res.frac = newfrac;
    //    return Res;
    //}

    //FP64 DivFP64(FP64 a, FP64 b) {
    //    if (a.is_null())
    //        return FP64();
    //    int expa = a.GetExp();
    //    int expb = b.GetExp();
    //    ull mb = b.GetFrac();
    //    int shift;
    //    FP64 bs = b;
    //    int tmp = (__lzcnt64(bs.GetFrac()) - 12);
    //    if (!b.is_subnorm()) {
    //        shift = 1022 - expb;
    //    }
    //    else {
    //        shift = 1022;
    //        bs.SetFrac(bs.GetFrac() << tmp);
    //    }
    //    double bsd = (double)b;
    //    bs.SetExp(bs.GetExp() + shift);
    //    bsd = set_exponent(b, expb + shift);
    //    FP64 c1 = (double)48 / 17; // 24 // потом заменить на константу
    //    FP64 c2 = (double)32 / 17; // 8 //
    //    double c1d = (double)48 / 17;
    //    double c2d = (double)32 / 17;
    //    FP64 x0 = FP64(fma(-(double)c2, (double)bs, (double)c1));// c1 - c2 *bs
    //    double x0d = c1d - c2d * bsd;
    //    FP64 x;
    //    double xd;
    //    FP64 bag = FP64(fma(-(double)x0, (double)bs, (double)2));// 2 - x0*bs
    //    double bagd = 2 - x0d * bsd;
    //    x = x0 * bag;
    //    xd = x0d * bagd;
    //    x *= FP64(fma(-(double)x, (double)bs, (double)2));
    //    xd *= (2 - xd * bsd);
    //    x *= FP64(fma(-(double)x, (double)bs, (double)2));
    //    xd *= (2 - xd * bsd);
    //    x *= FP64(fma(-(double)x, (double)bs, (double)2));
    //    xd *= (2 - xd * bsd);
    //    FP64 xh = x;
    //    double xhd = xd;
    //    FP64 tfma = FP64(fma(double(xh), double(bs), (double)-1));
    //    double tfmad = fma(xhd, bsd, (double)-1);
    //    FP64 xl = -xh * tfma;
    //    double xld = -xhd * tfmad;
    //    FP64 as = a;
    //    double asd = (double)a;
    //    as.SetExp(as.GetExp() + shift + tmp);
    //    set_exponent(a, get_exponent(a) + shift);
    //    FP64 temp = as * xl;
    //    double tempd = asd * xld;
    //    FP64 res = FP64(fma(double(as), double(xh), double(temp)));
    //    double resd = fma(asd, xhd, tempd);
    //    return res;
    //}
    FP64 DivFP64(FP64 a, FP64 b) {
        FP64 tmp;
        if (a.exp >= 2047 || b.exp >= 2047 || b.is_null() || a.is_null()) {
            if (a.exp >= 2047) { 
                if (a.GetFrac() != 0) {
                    tmp = a;
                    tmp.sign = a.sign ^ b.sign;
                    return tmp;
                }
                else if (b.GetExp() == 2047) {
                    tmp.exp = 2047;
                    tmp.frac = 1;
                    return tmp; 
                }
                else {
                    tmp = a;
                    tmp.sign = a.sign ^ b.sign;
                    return tmp;
                }
            }
            if (b.exp == 2047) { 
                if (b.GetFrac() != 0) {
                    tmp = b;
                    tmp.sign = a.sign ^ b.sign;
                    return tmp;
                    return tmp;
                }
                else return FP64();
            }

            if (b.is_null()) {
                if (a.is_null()) {
                    tmp.exp = 2047;
                    tmp.frac = 1;
                    tmp.sign = a.sign ^ b.sign;
                    return tmp;
                }
                else {
                    tmp.exp = 2047;
                    tmp.frac = 0;
                    tmp.sign = a.sign ^ b.sign;
                    return tmp;
                }
            }

            if (a.is_null()) {
                tmp.exp = 0;
                tmp.frac = 0;
                tmp.sign = a.sign ^ b.sign;
                return tmp;
            }
        }
        if (b.sign) {
            a.sign++;
            b.sign++;
        }
        int expa = a.GetExp();
        int expb = b.GetExp();
        int shift_a, shift_b;
        FP64 bs = b;
        if (!b.is_subnorm()) {
            shift_a = 1022 - expb;
            bs.SetExp(1022);
        }
        else {
            int sh = __lzcnt64(bs.GetFrac());
            //int sh = __builtin_clzll(bs.GetFrac());
            bs.SetFrac(bs.GetFrac() << (sh - 11));
            if (sh == 63) {
                bs.SetExp(1023);
                shift_a = 1022 + sh - 11;
            }
            else {
                bs.SetExp(1022);
                shift_a = 1021 + sh - 11;
            }


        }
        FP64 c1;
        c1.SetExp(1024);
        c1.SetFrac(1854423375976087);
        FP64 c2;
        c2.SetExp(1023);
        c2.SetFrac(3973764377091614);
        FP64 x0 = FP64(fma(-(double)c2, (double)bs, (double)c1));// c1 - c2 *bs
        FP64 x;
        FP64 bag = FP64(fma(-(double)x0, (double)bs, (double)2));// 2 - x0*bs
        x = x0 * bag;
        x *= FP64(fma(-(double)x, (double)bs, (double)2));
        x *= FP64(fma(-(double)x, (double)bs, (double)2));
        x *= FP64(fma(-(double)x, (double)bs, (double)2));
        x *= FP64(fma(-(double)x, (double)bs, (double)2));
        x *= FP64(fma(-(double)x, (double)bs, (double)2));
        FP64 xh = x;
        FP64 tfma = FP64(fma(double(xh), double(bs), (double)-1));
        FP64 xl = -xh * tfma;
        FP64 as = a;
        if (a.is_subnorm()) {
            shift_a -= (__lzcnt64(a.GetFrac()) - 12);
            as.SetFrac(as.GetFrac() << (__lzcnt64(as.GetFrac()) - 11));
        }
        as.SetExp(as.GetExp() + shift_a);
        FP64 temp = as * xl;
        FP64 res = FP64(fma(double(as), double(xh), double(temp)));

        return res;
    }
};


void print(ull a1, ull a2, ull a3, ull a4, ull a5, ull a6) {
    printf("%.8x", uint32(a1));
    printf("%.8x", uint32(a2));
    printf("%.8x", uint32(a3));
    printf("%.8x", uint32(a4));
    printf("%.8x", uint32(a5));
    printf("%.8x", uint32(a6));
    printf("\n");
}


FP64 fmaFP64(FP64 A, FP64 B, FP64 C) {
    if (A.is_null() || B.is_null())
        return C;
    if (C.is_null())
        return A * B;
    FP64 Res;
    ull m1 = A.GetFrac() + (1ULL << 52);
    ull m2 = B.GetFrac() + (1ULL << 52);
    ull m3 = C.GetFrac() + (1ULL << 52);
    ull MMSb = 0xFFFFFFFF00000000; // маска на старшие 32 бита
    ull MLSb = 0xFFFFFFFF; // маска на младшие 32 бита
    ull p = 53;
    ll ec = C.GetExp();
    ll eab = A.GetExp() + B.GetExp() - 1023;
    //ll eab = A.exp + B.exp - 1023 + 2;
    ll na = 1, nb = 1, nc = 1;
    if (A.is_subnorm()) {
        na = (12 - __lzcnt64(A.GetFrac()));
        //na = (12 - __builtin_clzll(A.GetFrac()));
        m1 -= (1ULL << 52);
    }
    if (B.is_subnorm()) {
        nb = (12 - __lzcnt64(B.GetFrac()));
        //nb = (12 - __builtin_clzll(B.GetFrac()));
        m2 -= (1ULL << 52);
    }
    if (C.is_subnorm()) {
        nc = (12 - __lzcnt64(C.GetFrac()));
        //nc = (12 - __builtin_clzll(C.GetFrac()));
        m3 -= (1ULL << 52);
    }
    ll d = ec - eab - 1 - nc + na + nb; // -3
    ull a1, a2, a3, a4, a5, a6;
    ull b1, b2, b3, b4, b5, b6;
    ull ab1 = 0, ab2 = 0, ab3 = 0, ab4 = 0, ab5 = 0, ab6 = 0;
    // число хранится так: ab3::ab4::ab5::ab6::ab1::ab2; ab3 - старший
    a1 = a2 = a3 = a4 = 0;
    a6 = (m1 & 0x00000000ffffffff);
    a5 = (m1 >> 32ULL);
    ull firstbit = 0;
    b1 = b2 = b3 = b4 = 0;
    b6 = (m2 & 0x00000000ffffffff);
    b5 = (m2 >> 32ULL);
    //
    ab6 = a6 * b6;
    ab5 += (ab6 >> 32);
    ab6 = (ab6 & MLSb);
    // заполнили ab6
    ////ab5 += (a5 * b5); // всё хорошо, так как изначально a5 и b5 - 20-битные, переполнения не должно быть
    ab5 += a5 * b6;
    ab5 += a6 * b5;
    ab4 += (ab5 >> 32);
    ab5 = (ab5 & MLSb);
    // заполнили ab5
    ab4 += a5 * b5;
    ab3 += (ab4 >> 32);
    ab4 = (ab4 & MLSb);

    cout << "ab:" << endl;
    print(ab3, ab4, ab5, ab6, ab1, ab2);
    // заполнили ab4 и ab3
    ull razrab; // не может быть 0, если она 0, то фактически a * b == 0;
    if (ab3)
        razrab = 4;
    else {
        if (ab4)
            razrab = 3;
        else {
            if (ab5)
                razrab = 2;
            else razrab = 1;
        }
    }
    // для начала надо сдвинуть всё число влево, чтобы первым битом была 1
    firstbit = 0;
    if (razrab == 4) {
        firstbit = __lzcnt64(ab3) - 32; // 32 - (64 - __lzcnt64(ab3)) показывает то, насколько надо сдвинуть влево
        //cout << firstbit << endl;
        if (firstbit) {
            ab3 = ab3 << firstbit;
            ull mask = 0xFFFFFFFFULL << (32 - firstbit); // маска на старшие firstbit числа типа uint32_t
            ab3 += ((ab4 & mask) >> (32 - firstbit));
            //ab3 = ab3 & MLSb;
            ab4 = ab4 << firstbit;
            ab4 = ab4 & MLSb; // обнуляем старшие биты ab4
            ab4 += ((ab5 & mask) >> (32 - firstbit));
            ab5 = ab5 << firstbit;
            ab5 = ab5 & MLSb;
            ab5 += ((ab6 & mask) >> (32 - firstbit));
            ab6 = ab6 << firstbit;
            ab6 = ab6 & MLSb;
        }
    }
    else if (razrab == 3) {
        //firstbit = __lzcnt64(ab4) - 32;
        ab3 = ab4;
        ab4 = ab5;
        ab5 = ab6;
        ab6 = 0;
        firstbit = __lzcnt64(ab3) - 32;
        if (firstbit) {
            ab3 = ab3 << firstbit;
            ull mask = 0xFFFFFFFFULL << (32 - firstbit);
            ab3 += ((ab4 & mask) >> (32 - firstbit));
            ab4 = ab4 << firstbit;
            ab4 = ab4 & MLSb; 
            ab4 += ((ab5 & mask) >> (32 - firstbit));
            ab5 = ab5 << firstbit;
            ab5 = ab5 & MLSb;
        }
    }
    else if (razrab == 2) {
        ab3 = ab5;
        ab4 = ab6;
        ab5 = ab6 = 0;
        firstbit = __lzcnt64(ab3) - 32;
        if (firstbit) {
            ab3 = ab3 << firstbit;
            ull mask = 0xFFFFFFFFULL << (32 - firstbit);
            ab3 += ((ab4 & mask) >> (32 - firstbit));
            ab4 = ab4 << firstbit;
            ab4 = ab4 & MLSb;
        }
    }
    else {
        ab3 = ab6;
        ab4 = ab5 = ab6 = 0;
        firstbit = __lzcnt64(ab3) - 32;
        if (firstbit) {
            ab3 = ab3 << firstbit;
        }
    }
    cout << "after the left shift ab: " << endl;
    print(ab3, ab4, ab5, ab6, ab1, ab2);

    // сдвигаем на p+3 = 56 разрядов вправо

    ull orig_ab3 = ab3;
    ull orig_ab4 = ab4;
    ull orig_ab5 = ab5;
    ull orig_ab6 = ab6;
    ull orig_ab1 = ab1;
    ull orig_ab2 = ab2;
    // сначало сдвинем на 32 бита вправо, потом ещё на 24
    ab3 = 0;
    ab4 = orig_ab3 & MLSb;
    ab5 = orig_ab4 & MLSb;
    ab6 = orig_ab5 & MLSb;
    ab1 = orig_ab6 & MLSb;  // ab2 = 0
    //
    orig_ab3 = ab3;
    orig_ab4 = ab4;
    orig_ab5 = ab5;
    orig_ab6 = ab6;
    orig_ab1 = ab1;
    orig_ab2 = ab2;
    // 0xFFFFFF - маска на 24 единицы
    ab4 = ab4 >> 24ULL;
    orig_ab5 = ab5 + ((orig_ab4 & 0xFFFFFF) << 32ULL);
    ab5 = orig_ab5 >> 24ULL;
    orig_ab6 = ab6 + ((orig_ab5 & 0xFFFFFF) << 32ULL);
    ab6 = orig_ab6 >> 24ULL;
    orig_ab1 = ab1 + ((orig_ab6 & 0xFFFFFF) << 32ULL);
    ab1 = orig_ab1 >> 24ULL;
    orig_ab2 = (orig_ab1 & 0xFFFFFF) << 32ULL;
    ab2 = orig_ab2 >> 24ULL; // на самом деле ab2 должен оставаться 0
    // завершили сдвиг числа на 56 битов вправо

    // два нуля сзади добавляются автоматически (они уже есть)

    cout << "ab>>32:" << endl;
    print(orig_ab3, orig_ab4, orig_ab5, orig_ab6, orig_ab1, orig_ab2);
    cout << "ab>>56:" << endl;
    print(ab3, ab4, ab5, ab6, ab1, ab2);

    ll dstrih;
    if (d <= -105) { // и в этом
        dstrih = 3 * 53 + 4;
        Res.SetExp(eab);
    }
    if (d > -105 && d <= 2) { // в этом случае проблема //
        dstrih = 56 - d; // 59
        Res.SetExp(eab);
    }
    if (d > 2 && d <= 55) {
        dstrih = 56 - d;
        Res.SetExp(ec);
    }
    if (d >= 56) {
        dstrih = 0;
        Res.SetExp(ec);
    }

    // сдвигаем c вправо на dstrih бит

    ull number = dstrih / 32; // 2
    ull position = dstrih % 32; // 3
    // выгодно было бы превратить представление 192-битного числа в массив из 6 чисел для избежания if

    ull c[6] = { 0 };
    //c0::c1::c2::c3::c4::c5 - хранится
    //c6 = (m3 & 0x00000000ffffffff);
    //c5 = (m3 >> 32ULL);
    c[1] = (m3 & 0x00000000ffffffff); 
    c[0] = (m3 >> 32ULL);
    ull razrc;
    if (c[0])
        razrc = 2;
    else razrc = 1;
    firstbit = 0;
    if (razrc == 2) {
        firstbit = __lzcnt64(c[0]) - 32;
        if (firstbit) {
            c[0] = c[0] << firstbit;
            ull mask = 0xFFFFFFFFULL << (32 - firstbit); 
            c[0] += ((c[1] & mask) >> (32 - firstbit));
            c[1] = c[1] << firstbit;
            c[1] = c[1] & MLSb; 
        }
    }
    else {
        c[0] = c[1];
        c[1] = 0;
        firstbit = __lzcnt64(c[0]) - 32;
        if (firstbit) {
            c[0] = c[0] << firstbit;
        }
    }
    cout << "before the shift c: " << endl;
    print(c[0], c[1], c[2], c[3], c[4], c[5]);
    ull available_bits = 32 - position; // 5
    ull temp, round = 0;

    /*if (dstrih != 163 && !C.is_subnorm()) {
        c[number] = (((~0ULL) << (53 - available_bits)) & C.frac); // старшие available_bits мантисы C
        c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & C.frac) & (((~0ULL) << (53 - available_bits - 32)));
        if (53 - available_bits - 32 > 0)
            c[number + 2] = ((1ULL << (53 - available_bits - 32)) - 1) & C.frac;
    }
    else {
        round = C.frac;
    }*/
    //dstrih += 23;
    // можно без if написать
    ull newc[6] = { 0 };
    ull rfc = 64 - __lzcnt64(m3);
    //ull rfc = 64 - __builtin_clzll(m3); // разрядность мантисы C до разложения на 32-битные числа
    if (dstrih != 163) {
        if (!C.is_subnorm()) {
            //c[number] = (((~0ULL) << (53 - available_bits)) & m3) >> (53 - available_bits); // старшие available_bits мантисы C
            //if (c[number] >= ((1ULL << 32) - 1))
            //    cout << "bad1" << endl;
            //if (53 - available_bits - 32 >= 0) {// =?
            //    c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & m3) >> (53 - available_bits - 32);// >> (53 - available_bits - 32);
            //    //c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & m3) >>
            //    c[number + 2] = (((1ULL << (53 - available_bits - 32)) - 1) & m3) << (53 - available_bits - 32); // возможно 32 - 53 + available_bits
            //}
            //else {
            //    c[number + 1] = (((1ULL << (53 - available_bits)) - 1) & m3) << (32ULL - 53 + available_bits);
            //}        
            for (int i = 5 - number; i >= 0; i--) 
                c[i + number] = c[i];
            for (int i = 0; i < number; i++) 
                c[i] = 0;
            //c[number] = c[number] & MLSb;
            ull orig_c[6];
            for (int i = 0; i < 6; i++)
                orig_c[i] = c[i];
            c[number] = c[number] >> position;
            if (number + 1 < 6) {
                orig_c[number + 1] = c[number + 1] + ((orig_c[number] & ((1ULL << position) - 1)) << 32ULL);
                c[number + 1] = orig_c[number + 1] >> position;
            }
            if (number + 2 < 6) { // не уверен на счёт этих действий
                orig_c[number + 2] = c[number + 2] + ((orig_c[number + 1] & ((1ULL << position) - 1)) << 32ULL);
                c[number + 2] = orig_c[number + 2] >> position;
            }
            //orig_ab1 = ab1 + ((orig_ab6 & 0xFFFFFF) << 32ULL);
            //ab1 = orig_ab1 >> 24ULL;
            //orig_ab2 = (orig_ab1 & 0xFFFFFF) << 32ULL;
            //ab2 = orig_ab2 >> 24ULL; // на самом деле ab2 должен оставаться 0
            //for (uint32_t i = number; i < 6; ++i) {
            //    uint32_t original_index = i - number;
            //    if (original_index >= 6) continue;

            //    // Текущий 32-битный блок (учитываем только младшие 32 бита)
            //    uint64_t current = c[original_index] & 0xFFFFFFFF;
            //    // Сдвигаем текущий блок на r бит вправо
            //    uint64_t shifted = current >> position;

            //    // Если есть следующий блок, добавляем его младшие r бит
            //    if (original_index + 1 < 6) {
            //        uint64_t next = c[original_index + 1] & 0xFFFFFFFF;
            //        shifted |= (next & ((1ULL << position) - 1)) << (32 - position);
            //    }

            //    // Сохраняем результат, обрезая до 32 бит
            //    newc[i] = shifted & 0xFFFFFFFF;
            //}
        }
        else {
            if (rfc >= available_bits) {
                c[number] = (((~0ULL) << (rfc - available_bits)) & m3) >> (rfc - available_bits); // старшие available_bits мантисы C // добавил сдвиг вправо
                //if (c[number] >)
                if (rfc - available_bits >= 32) {
                    c[number + 1] = (((1ULL << (rfc - available_bits)) - 1) & m3) >> (rfc - available_bits - 32);
                    c[number + 2] = (((1ULL << (rfc - available_bits - 32)) - 1) & m3) << (rfc - available_bits - 32);
                }
                else {
                    c[number + 1] = (((1ULL << (rfc - available_bits)) - 1) & m3); // сдвинуть влево
                }
            }
            else {
                c[number] = ((1ULL << rfc) - 1) & m3; // не уверен, что идет в ячейку number
            }
        }
    }
    else { // если dstrih == 163
        if (!C.is_subnorm()) {
            round += (m3 & ((1ULL << 52) - 1));
            // последние 52бита уходят в round (биты округления)
            c[5] = 1; // так как мантиса нормального числа представляется из себя выражение 2^52 + x, то старший бит всегда будет равен 1
        }
        else {
            round += (m3 & ((1ULL << (rfc - 1)) - 1));
            // c полностью уйдет в 0
        }
    }
    //
   /* ull c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
    c1 = c2 = c3 = c4 = 0;
    c6 = (m3 & 0x00000000ffffffff);
    c5 = (m3 >> 32ULL);*/
    //
    //cout << "before the shift c: " << endl;
    //print(c3, c4, c5, c6, c1, c2);
    // cout << "before the shift: " << endl;
     //cout << m3 << endl;
    cout << "c:" << endl;
    print(c[0], c[1], c[2], c[3], c[4], c[5]);
    //cout << c[0] << endl;
    //cout << ab3 << endl;
    // для начала разберемся со сложением
    ull s[6] = { 0 }; // сумма (a*b) + c
    //// лидирующий ненулевой бит c должен находиться в числе c[number] на позиции position
    s[5] = ab2 + c[5];
    s[4] = ab1 + c[4] + ((s[5] & MMSb) >> 32);
    s[3] = ab6 + c[3] + ((s[4] & MMSb) >> 32);
    s[2] = ab5 + c[2] + ((s[3] & MMSb) >> 32);
    s[1] = ab4 + c[1] + ((s[2] & MMSb) >> 32);
    s[0] = ab3 + c[0] + ((s[1] & MMSb) >> 32);

    cout << "s:" << endl;
    print(s[0], s[1], s[2], s[3], s[4], s[5]);


    //cout << s[0] << endl;
    // теперь необходимо сдвинуть сумму s влево так, чтобы лидирующий бит был 1, найдем величину этого сдвига
    ll dstrih2; // 55
    if (d <= 2) {
        ll l = 0; // кол-во ведущих нулей в 109 младших битах числа s
        // получается надо начинать смотреть с s[2]
       // cout << s[0] << endl;
      //  cout << (s[2] & ((1ULL << 13) - 1)) << endl;
        //if (s[2] & ((1ULL << 13) - 1)) {
        //    // можем смотреть в первых 13 битах
        //    ll temp = s[2] & ((1ULL << 13) - 1);
        //    l = __lzcnt64(temp) - 51;
        //    //l = __builtin_clzll(temp) - 51; // возможно поменять местами
        //    //cout << l << endl;
        //}
        //else {
        //    l += 13;
        //    if (s[3] != 0) {
        //        l += __lzcnt64(s[3]) - 32;
        //        //l += __builtin_clzll(s[3]) - 32; //
        //        //cout << __builtin_clzll(s[3]) - 32 << endl;
        //    }
        //    else {
        //        l += 32;
        //        if (s[4] != 0) {
        //            l += __lzcnt64(s[4]) - 32;
        //            //l += __builtin_clzll(s[4]) - 32;
        //        }
        //        else {
        //            l += 32;
        //            if (s[5] != 0) {
        //                l += __lzcnt64(s[5]) - 32;
        //                //l += __builtin_clzll(s[5]) - 32;
        //            }
        //            // иначе число c == 0
        //        }
        //    }
        //}
        
        //136 младших бит
        if (s[1] & ((1ULL << 8) - 1)) {
                ll temp = s[1] & ((1ULL << 8) - 1);
                l = __lzcnt64(temp) - 56;         
        }
        else {
            l += 8;
            if (s[2] != 0) {
                l += __lzcnt64(s[2]) - 32;
            }
            else {
                l += 32;
                if (s[3] != 0) {
                    l += __lzcnt64(s[3]) - 32;
                    //l += __builtin_clzll(s[4]) - 32;
                }
                else {
                    l += 32;
                    if (s[4] != 0) {
                        l += __lzcnt64(s[4]) - 32;
                        //l += __builtin_clzll(s[5]) - 32;
                    }
                    else {
                        l += 32;
                        if (s[5] != 0) {
                            l += __lzcnt64(s[5]) - 32;
                        }
                    }
                }
            }
        }
        //рассчитали l
        if (eab - l + 2 >= 1) { // в таком случае результат fma будет нормальным числом // здесь устанавливается неправильная экспонента
            dstrih2 = 55 + l;
            Res.SetExp(eab - l + 2);
            //cout << l << endl;
            //cout << eab << endl;
        }
        else { // иначе результат будет сабнормалом
            Res.SetExp(1);
            ll leadbits = 1 - (eab - l + 2); // столько лидирующих нулей получит мантиса результата
            dstrih2 = 57 - 1 + eab;
        }
    }
    else { // d > 2
        dstrih2 = dstrih;
        // в таком случае сдвиг влево отменяет изначальный сдвиг вправо, но после сдвига ведущая 1 может сменить позицию: при слоежении оказаться левее на 1 бит, при вычитании оказаться правее на 1 бит?
    }
    // я так понимаю после сложения a*b и c необходимо сначала свдинуть результать на dstrih2 влево, затем округлить взяв старшие 53бита, остальные добавить в round (они буддут составлять 110бит), сдвинув их при это на разрядность round уже имеющегося. Затем вычесть из мантисы мнимую единицу, то есть (1ULL << 52) и получившуюся мантису округлить по round.
    // наверно можно не двигать всё число, а найти лидирующий бит и свдинуть его и 52 бита за ним, а так же round
    number = dstrih2 / 32;
    position = dstrih2 % 32;
    // сначало сдвигаем цельными фрагментами
    ull temps[6] = { 0 };
    for (int i = 0; i < 6; ++i) {
        if (i + number < 6)
            temps[i] = (s[i + number] & 0xFFFFFFFF);
    }
    memcpy(s, temps, sizeof(temps));
    cout << "after full shift: " << endl;
    print(s[0], s[1], s[2], s[3], s[4], s[5]);
    ull sn[6];
    //memcpy(sn, s, sizeof(s));
    //position = 23;// 23
    if (position) {
        for (int i = 0; i < 5; i++) {
            sn[i] = (s[i] << position) & ((1ULL << 32) - 1);
            sn[i] += (s[i + 1] >> (32 - position));
        }
        sn[5] = 0;
        //for (int i = 0; i < 6; i++)
        //    temps[i] = s[i];
        //for (int i = 0; i < 6; i++) {
        //    // Биты, переносимые из следующего блока (carry), или 0 для последнего блока
        //    uint64_t carry = (i < 5) ? ((uint64_t)temps[i + 1] >> (32 - position)) : 0;
        //    // Сдвиг текущего блока влево на m бит и добавление carry
        //    s[i] = (((uint64_t)temps[i] << position) | carry);
        //}
        /*orig_c[number + 1] = c[number + 1] + ((orig_c[number] & ((1ULL << position) - 1)) << 32ULL);
        c[number + 1] = orig_c[number + 1] >> position;*/
        //ull orig_s[6];
        //memcpy(orig_s, s, sizeof(s));
        //ull mask = 0xFFFFFFFFULL << (32 - position);
        ///*for (int i = 5; i >= 1; i--) {
        //    orig_s[i - 1] = s[i - 1] + ((orig_s[i] & mask) >> 32ULL);
        //    s[i - 1] = orig_s[i - 1] << position;
        //}*/
        //ull carry = 0;
        //for (int i = 5; i >= 0; --i) {
        //    ull current = (s[i] & 0xFFFFFFFF); // на всякий случай
        //    ull shifted = (current << position) + carry;
        //    carry = (current << position) >> 32ULL;
        //    s[i] = (shifted & 0xFFFFFFFF);
        //}
    }
    /*if (d <= 2 && ((s[0] & (1ULL << 31)) == 0)) {
        Res.SetExp(Res.GetExp() - 1);
        ull carry = 0;
        for (int i = 5; i >= 0; --i) {
            ull current = (s[i] & 0xFFFFFFFF);
            ull shifted = (current << 1) + carry;
            carry = (current << 1) >> 32ULL;
            s[i] = (shifted & 0xFFFFFFFF);
        }
    }*/
    cout << (s[0] & (1ULL << 31)) << endl;
    cout << "after shift s: " << endl; // корректно происходит
    print(sn[0], sn[1], sn[2], sn[3], sn[4], sn[5]);
    // поиск лидирующей 1
    ull newfrac = 0;
    //for (int i = 0; i < 6; i++) { // возможно здесь проблема
    //    if (s[i]) {
    //        ull curnumb = i;
    //        ull cnt = 0;
    //        while (cnt < 53 && curnumb < 6) {
    //            available_bits = 64 - __lzcnt64(s[curnumb]);
    //            if (available_bits > 32) //
    //                available_bits = 32; //
    //            if (available_bits + cnt > 53)
    //                available_bits = 53 - cnt;
    //            newfrac += (s[curnumb] & ((1ULL << available_bits) - 1));
    //            curnumb++;
    //            cnt += available_bits;
    //        }
    //    }
    //}
    cout << (s[0] & (1ULL << 31)) << endl;
    newfrac += ((s[0] & MLSb) << 21ULL);
    ull mask = 0xFFFFFFFFULL << (32 - 21);
    newfrac += ((s[1] & mask) >> (32 - 21));
    //ull curnumb;
    //ull len;
    //for (int i = 0; i < 6; ++i) {
    //    if (s[i]) {
    //        curnumb = i;
    //        break;
    //    }
    //}
    //ull razrMSb = 64 - __lzcnt64(s[curnumb] & MLSb);
    ////ull razrMSb = 64 - __builtin_clzll(s[curnumb] & MLSb);
    //if (razrMSb >= 21) {
    //    if (curnumb < 5)
    //        len = 2;
    //    else len = 1;
    //}
    //else {
    //    if (curnumb < 4)
    //        len = 3;
    //    if (curnumb == 4)
    //        len = 2;
    //    if (curnumb == 5)
    //        len = 1;
    //}
    //for (int i = 0; i < len; ++i) { // лишний
    //    if (i == 0) {
    //        //newfrac += ((s[curnumb] & ((1ULL << razrMSb) - 1)) << (53 - razrMSb)); // надо учитывать len
    //        newfrac += ((s[curnumb] & ((1ULL << razrMSb) - 1)) << 32);
    //        curnumb++;
    //        if (curnumb > 5)
    //            break;
    //    }
    //    if (i == 1) {
    //        //newfrac += ((s[curnumb] & MLSb) << (53 - razrMSb - 32));
    //        newfrac += (s[curnumb] & MLSb);
    //        curnumb++;
    //        if (curnumb > 5)
    //            break;
    //    }
    //    if (i == 2) { // отсюда надо взять 53 - 32 - razrMSb страрших бит
    //        newfrac += ((s[curnumb] >> (53ULL - 32 - razrMSb)));
    //    }
    //}
    cout << "my frac from s: " << newfrac << endl;
    //newfrac = 4503604187917403; // для конкретного примера
    //if (newfrac >= (1ULL << 53) || newfrac < (1ULL << 52))
     //   cout << "bred" << endl;
    if (newfrac >= (1ULL << 53)) {
        newfrac >>= 1ULL;
        Res.SetExp(Res.GetExp() + 1);
    }
    if (newfrac < (1ULL << 52)) {
        newfrac <<= 1ULL;
        Res.SetExp(Res.GetExp() - 1);
    }
    //cout << newfrac << endl;
    //ull debug = newfrac - (1ULL << 52);
   /* while (newfrac < (1LL << 52LL) && Res.GetExp() > 1) {
        newfrac = newfrac << 1ULL;
        Res.SetExp(Res.GetExp() - 1);
    }*/
    ull debug = newfrac - (1ULL << 52);
    // проверка (пока без округления)
    Res.SetSign(0);
    Res.SetFrac(debug);
    return Res;
}



