#include <iostream>
#include <cmath>

#if !__GNUC__ && !__clang__
int32_t __builtin_clz(uint32_t num) { \
	if (num == 0) return 32; \
	int32_t res = 0; \
	while (!(num & 0x8000'0000)) { \
		num <<= 1; \
		++res; \
	} \
	return res; \
} 

int32_t __builtin_clzll(uint64_t num) { \
	if (num == 0) return 64; \
	int32_t res = 0; \
	while (!(num & 0x8000'0000'0000'0000)) { \
		num <<= 1; \
		++res; \
	} \
	return res; \
} 
#endif

class BF16 {
    private:
        uint16_t data;
    
        uint16_t add(uint16_t l, uint16_t r) const noexcept;
    
        uint16_t usub(uint16_t r) const noexcept;
    
        uint16_t sub(uint16_t l, uint16_t r) const noexcept;
    
        uint16_t mul(uint16_t l, uint16_t r) const noexcept;
    
        uint16_t fma(uint16_t a, uint16_t b, uint16_t c) const noexcept;
    
        uint16_t div(uint16_t l, uint16_t r) const noexcept;
    
    public:
        BF16() = default;
        BF16(uint16_t t) noexcept;
        BF16(float f) noexcept;
        BF16(double d) noexcept;
        BF16(const BF16& fp) noexcept;
        BF16& operator=(const BF16& fp) noexcept;
        BF16& operator=(uint16_t t) noexcept;
        BF16& operator=(float f) noexcept;
        BF16& operator=(double f) noexcept;
        operator float() const noexcept;
        operator double() const noexcept;
    
        BF16 operator- () noexcept;
        BF16& operator+= (const BF16& fp) noexcept;
        BF16& operator-= (const BF16& fp) noexcept;
        BF16& operator*= (const BF16& fp) noexcept;
        BF16& operator/= (const BF16& fp) noexcept;
        BF16 operator+ (const BF16& fp) const noexcept;
        BF16 operator- (const BF16& fp) const noexcept;
        BF16 operator* (const BF16& fp) const noexcept;
        BF16 operator/ (const BF16& fp) const noexcept;
        BF16 fma (const BF16& a, const BF16& b) const noexcept;
        friend std::ostream& operator << (std::ostream& stream, const BF16& fp);
        friend std::istream& operator >> (std::istream& stream, BF16& fp);
    
        bool operator == (const BF16& fp) const noexcept;
        bool operator < (const BF16& fp) const noexcept;
        bool operator > (const BF16& fp) const noexcept;
        bool operator >= (const BF16& fp) const noexcept;
        bool operator <= (const BF16& fp) const noexcept;
        bool operator != (const BF16& fp) const noexcept;
    };