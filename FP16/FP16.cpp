#include <stdfloat>
#include <iostream>

using float16_t = std::float16_t;

//распределение бит для FP16 1:5:10
const int manLength = 10; //биты для мантиссы 0-1023
const int expLength = 5; //биты для экспоненты 0-31
const int signLength = 1; //биты для знаков 0-1
const int shiftExp = (1 << (expLength - 1)) - 1;  //смещение мантиссы на 15

class FP16 {
private:
    unsigned int man : manLength; //мантисса       %1024
    unsigned int exp : expLength; //порядок         %32
    unsigned int sign : signLength; //знак           %2

    //нормал: 2^(exp-15) * 1.man, exp!=0
    //сабнормал: 2^-14 * 0.man, exp=0

public:
    FP16(int _sign = 0, int _exp = 0, int _man = 0) :sign(_sign), exp(_exp), man(_man) {}

    bool IsSubnormal() const {
        if ((exp == 0) && (man != 0)) return true;
        else return false;
    }

    bool IsInf() const {
        if (exp == ((1 << expLength) - 1)) return true;
        else return false;
    }

    bool IsNull() const noexcept {
        if ((exp == 0) && (man == 0)) return true;
        return false;
    }

    void printFP16_bites() {
        std::cout << sign << "_";
        for (char i = 0; i < expLength; ++i) {
            if ((exp & (1 << (expLength - i - 1))) == 0) std::cout << 0;
            else std::cout << 1;
        }
        std::cout << "_";
        for (char i = 0; i < manLength; ++i) {
            if ((man & (1 << (manLength - i - 1))) == 0) std::cout << 0;
            else std::cout << 1;
        }
    }

    //
    FP16& operator=(const FP16& N) {
        if (this != &N) {
            sign = N.sign;
            exp = N.exp;
            man = N.man;
        }
        return *this;
    }


    //
    FP16 operator+(FP16 N2) {
        if ((IsSubnormal() + N2.IsSubnormal()) == 0) {
            return AddFP16_N_N(*this, N2);
        }
        else {
            if ((IsSubnormal() + N2.IsSubnormal()) == 1) return AddFP16_N_S(*this, N2);
            else return AddFP16_S_S(*this, N2);
        }
    }

    //
    FP16 operator-(FP16 N2) {
        FP16 temp(N2);
        ++temp.sign;
        return (*this + temp);
    }
    //
    FP16 operator-() {
        FP16 temp;
        temp = *this;
        ++temp.sign;
        return (temp);
    }


};





void print_float16_bites(float16_t f) {
    uint8_t bytes[2];
    bool bites[16];

    for (int i = 0; i < sizeof(f); i++) bytes[i] = *((uint8_t*)(&f) + i);

    for (int i = 0; i < 16; ++i) bites[15 - i] = (bytes[i / 8] >> (i % 8)) & 1;

    for (int i = 0; i < 16; ++i) {
        if (i == 1 || i == 6) std::cout << "_";
        std::cout << bites[i];
    }
}


unsigned short float16_bites_to_short(float16_t f) {
    unsigned short temp = 0;
    uint8_t bytes[2];
    bool bites[16];

    for (int i = 0; i < sizeof(f); i++) bytes[i] = *((uint8_t*)(&f) + i);

    for (int i = 0; i < 16; ++i) bites[15 - i] = (bytes[i / 8] >> (i % 8)) & 1;

    for (int i = 0; i < 16; ++i) temp |= (bites[15 - i] << i);

    return temp;
}






using namespace std;
int main() {
    float A = 1;
    float B = 0;
    float16_t a = A;
    float16_t b = B;



    cout << a + b << endl;
    return 0;
}



