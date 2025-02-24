#include <iostream>
class FP16 {
private:
	uint16_t data;

	uint16_t add(uint16_t l, uint16_t r) const noexcept {
		
	}

	uint16_t sub(uint16_t l, uint16_t r) const noexcept {
		
	}

	uint16_t mul(uint16_t l, uint16_t r) const noexcept {
		
	}

	uint16_t fma(uint16_t a, uint16_t b, uint16_t c) const noexcept {

	}

	uint16_t div(uint16_t l, uint16_t r) const noexcept {

	}

public:
	FP16() = default;
	FP16(uint16_t t) noexcept : data(t) {}
	FP16(float f) noexcept {
		
	}
	FP16(double f) noexcept {
		
	}
	FP16(const FP16& fp) noexcept : data(fp.data) {}
	FP16& operator=(const FP16& fp) {
		data = fp.data;
		return *this;
	}
	FP16& operator=(uint16_t t) {
		data = t;
		return *this;
	}
	FP16& operator=(float f) {
		data = FP16(f).data;
		return *this;
	}
    FP16& operator=(double f) {
		data = FP16(f).data;
		return *this;
	}
	operator float() const noexcept {
		
	}
	operator double() const noexcept {
		
	}

	FP16 operator- () {
		uint16_t tmp = (data & 0x7FFF) + (~data & 0x8000);
		return FP16(tmp);
	}
	FP16& operator+= (const FP16& fp) {
		data = add(data, fp.data);
		return *this;
	}
	FP16& operator-= (const FP16& fp) {
		data = sub(data, fp.data);
		return *this;
	}
	FP16& operator*= (const FP16& fp) {
		data = mul(data, fp.data);
		return *this;
	}
	FP16& operator/= (const FP16& fp) {
		data = div(data, fp.data);
		return *this;
	}
	FP16 operator+ (const FP16& fp) const noexcept {
		FP16 tmp = *this;
		tmp += fp;
		return tmp;
	}
	FP16 operator- (const FP16& fp) const noexcept {
		FP16 tmp = *this;
		tmp -= fp;
		return tmp;
	}
	FP16 operator* (const FP16& fp) const noexcept {
		FP16 tmp = *this;
		tmp *= fp;
		return tmp;
	}
	FP16 operator/ (const FP16& fp) const noexcept {
		FP16 tmp = *this;
		tmp /= fp;
		return tmp;
	}
	FP16 fma (const FP16& a, const FP16& b) const noexcept {
		FP16 tmp = *this;
		tmp = fma(data, a.data, b.data);
		return tmp;
	}

	friend std::ostream& operator << (std::ostream& stream, const FP16& fp) {
		float tmp = fp;
		stream << tmp;
		return stream;
	}
	friend std::istream& operator >> (std::istream& stream, FP16& fp) {
		float tmp;
		stream >> tmp;
		fp = tmp;
		return stream;
	}

	bool operator == (const FP16& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		return data == fp.data;
	};
	bool operator < (const FP16& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP16 tmp = *this - fp;
		return tmp.data > 0x80000000ul;
	}
	bool operator > (const FP16& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP16 tmp = *this - fp;
		return (tmp.data < 0x80000000ul) && (tmp.data != 0ul);
	}
	bool operator >= (const FP16& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP16 tmp = *this - fp;
		return (tmp.data <= 0x80000000ul);
	}
	bool operator <= (const FP16& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP16 tmp = *this - fp;
		return (tmp.data >= 0x80000000ul) && (tmp.data == 0ul);
	}
	bool operator != (const FP16& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return true;
		return data != fp.data;
	}
};