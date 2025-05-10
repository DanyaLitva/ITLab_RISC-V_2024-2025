#include <iostream>
#include <cstdint>
#include <cstring>
#include <cmath>
#pragma once
const int manLength = 10;
const int expLength = 5;
const int signLength = 1;
const int shiftExp = (1 << (expLength - 1)) - 1;


class FP16 {
private:
    unsigned int man : manLength;
    unsigned int exp : expLength;
    unsigned int sign : signLength;

public:
    FP16(int _sign = 0, int _exp = 0, int _man = 0):sign(_sign), exp(_exp), man(_man) {}
    FP16(const FP16& right) :sign(right.sign), exp(right.exp), man(right.man) {}
    bool IsSubnormal() const noexcept;
    bool IsInf() const noexcept;
    bool IsNull() const noexcept;
    bool IsNan() const noexcept;
    void printFP16_bites() const noexcept;
    //operator float16_t() const noexcept;
    FP16& operator=(const FP16& N) noexcept;
    void set_from_uint(uint16_t temp) noexcept;
    uint16_t get_int() const noexcept;
    //FP16(uint16_t temp = 0) noexcept;
    FP16 operator+(FP16 right) const noexcept;
    FP16 operator-(FP16 right) const noexcept;
    FP16 operator*(FP16 right) const noexcept;
    friend uint16_t mul_half(uint16_t a, uint16_t b) noexcept;
    friend uint16_t usub_half(uint16_t temp) noexcept;
    friend FP16 fma4(FP16 a, FP16 b, FP16 c) noexcept;
    friend uint16_t fma_half(uint16_t a, uint16_t b, uint16_t c) noexcept;
    friend std::ostream& operator<<(std::ostream& os, FP16 num) noexcept;
    FP16 operator/(FP16 right) const noexcept;
    bool operator == (const FP16& right) const noexcept;
    bool operator != (const FP16& right) const noexcept;
    bool operator < (const FP16& right) const noexcept;
    bool operator > (const FP16& right) const noexcept;
    bool operator >= (const FP16& right) const noexcept;
    bool operator <= (const FP16& right) const noexcept;
    FP16& operator+=(const FP16& right) noexcept;
    FP16& operator-=(const FP16& right) noexcept;
    FP16& operator*=(const FP16& right) noexcept;
    FP16& operator/=(const FP16& right) noexcept;
    FP16 operator-() noexcept;
    FP16& operator=(float N) noexcept;
    FP16& operator=(double N) noexcept;
    FP16(float num) noexcept{
        *this = num;
    }
    FP16(double num) noexcept{
        *this = num;
    }
    FP16& fma(const FP16& a, const FP16& b){
        *this = fma4(*this,a,b);
        return *this;
    }
    operator double() const noexcept;
    operator float() const noexcept;
};


float fp16_to_float(uint16_t h) noexcept;
std::ostream& operator<<(std::ostream& os, FP16 num) noexcept;
