#include <iostream>

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

class FP32 {
private:
	uint32_t data;

	uint32_t add(uint32_t l, uint32_t r) const noexcept;

	uint32_t sub(uint32_t l, uint32_t r) const noexcept;

	uint32_t usub(uint32_t a) const noexcept;

	uint32_t mul(uint32_t l, uint32_t r) const noexcept;

	uint32_t fma(uint32_t a, uint32_t b, uint32_t c) const noexcept;

	uint32_t div(uint32_t l, uint32_t r) const noexcept;

public:
	FP32() = default;
	FP32(uint32_t t);
	FP32(float f) noexcept;
	FP32(double f) noexcept;
	FP32(const FP32& fp);
	FP32& operator=(const FP32& fp) noexcept;
	FP32& operator=(uint32_t t) noexcept;
	FP32& operator=(float f) noexcept;
	FP32& operator=(double f) noexcept;
	operator float() const noexcept;
	operator double() const noexcept;
	FP32 operator- () noexcept;
	FP32& operator+= (const FP32& fp) noexcept;
	FP32& operator-= (const FP32& fp) noexcept;
	FP32& operator*= (const FP32& fp) noexcept;
	FP32& operator/= (const FP32& fp) noexcept;
	FP32 operator+ (const FP32& fp) const noexcept;
	FP32 operator- (const FP32& fp) const noexcept;
	FP32 operator* (const FP32& fp) const noexcept;
	FP32 operator/ (const FP32& fp) const noexcept;
	FP32 fma (const FP32& a, const FP32& b) const noexcept;

	friend std::ostream& operator << (std::ostream& stream, const FP32& fp);
	friend std::istream& operator >> (std::istream& stream, FP32& fp);

	bool operator == (const FP32& fp) const noexcept;
	bool operator < (const FP32& fp) const noexcept;
	bool operator > (const FP32& fp) const noexcept;
	bool operator >= (const FP32& fp) const noexcept;
	bool operator <= (const FP32& fp) const noexcept;
	bool operator != (const FP32& fp) const noexcept;
};