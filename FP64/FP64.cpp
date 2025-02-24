#include <iostream>
class FP64 {
private:
	uint64_t data;

	uint64_t add(uint64_t l, uint64_t r) const noexcept {
		
	}

	uint64_t sub(uint64_t l, uint64_t r) const noexcept {
		
	}

	uint64_t mul(uint64_t l, uint64_t r) const noexcept {
		
	}

	uint64_t fma(uint64_t a, uint64_t b, uint64_t c) const noexcept {

	}

	uint64_t div(uint64_t l, uint64_t r) const noexcept {

	}

public:
	FP64() = default;
	FP64(uint64_t t) noexcept : data(t) {}
	FP64(float f) noexcept {
        double df = f;
		data = reinterpret_cast<uint64_t*>(&df)[0];
	}
    FP64(double f) noexcept {
		data = reinterpret_cast<uint64_t*>(&f)[0];
	}
	FP64(const FP64& fp) noexcept : data(fp.data) {}
	FP64& operator=(const FP64& fp) {
		data = fp.data;
		return *this;
	}
	FP64& operator=(uint64_t t) {
		data = t;
		return *this;
	}
	FP64& operator=(float f) {
		data = FP64(f).data;
		return *this;
	}
    FP64& operator=(double f) {
		data = FP64(f).data;
		return *this;
	}
    operator float() const noexcept {
		uint64_t tmp = data;
		return float(reinterpret_cast<double*>(&tmp)[0]);
	}
	operator double() const noexcept {
		uint64_t tmp = data;
		return reinterpret_cast<double*>(&tmp)[0];
	}

	FP64 operator- () {
		uint64_t tmp = (data & 0x7FFF'FFFF'FFFF'FFFF) + (~data & 0x8000'0000'0000'0000);
		return FP64(tmp);
	}
	FP64& operator+= (const FP64& fp) {
		data = add(data, fp.data);
		return *this;
	}
	FP64& operator-= (const FP64& fp) {
		data = sub(data, fp.data);
		return *this;
	}
	FP64& operator*= (const FP64& fp) {
		data = mul(data, fp.data);
		return *this;
	}
	FP64& operator/= (const FP64& fp) {
		data = div(data, fp.data);
		return *this;
	}
	FP64 operator+ (const FP64& fp) const noexcept {
		FP64 tmp = *this;
		tmp += fp;
		return tmp;
	}
	FP64 operator- (const FP64& fp) const noexcept {
		FP64 tmp = *this;
		tmp -= fp;
		return tmp;
	}
	FP64 operator* (const FP64& fp) const noexcept {
		FP64 tmp = *this;
		tmp *= fp;
		return tmp;
	}
	FP64 operator/ (const FP64& fp) const noexcept {
		FP64 tmp = *this;
		tmp /= fp;
		return tmp;
	}
	FP64 fma (const FP64& a, const FP64& b) const noexcept {
		FP64 tmp = *this;
		tmp = fma(data, a.data, b.data);
		return tmp;
	}

	friend std::ostream& operator << (std::ostream& stream, const FP64& fp) {
		double tmp = fp;
		stream << tmp;
		return stream;
	}
	friend std::istream& operator >> (std::istream& stream, FP64& fp) {
		double tmp;
		stream >> tmp;
		fp = tmp;
		return stream;
	}

	bool operator == (const FP64& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		return data == fp.data;
	};
	bool operator < (const FP64& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP64 tmp = *this - fp;
		return tmp.data > 0x80000000ul;
	}
	bool operator > (const FP64& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP64 tmp = *this - fp;
		return (tmp.data < 0x80000000ul) && (tmp.data != 0ul);
	}
	bool operator >= (const FP64& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP64 tmp = *this - fp;
		return (tmp.data <= 0x80000000ul);
	}
	bool operator <= (const FP64& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return false;
		FP64 tmp = *this - fp;
		return (tmp.data >= 0x80000000ul) && (tmp.data == 0ul);
	}
	bool operator != (const FP64& fp) const noexcept {
		if (((data & 0x7FFFFFFF) > 0x7F800000) || ((fp.data & 0x7FFFFFFF) > 0x7F800000)) return true;
		return data != fp.data;
	}
};