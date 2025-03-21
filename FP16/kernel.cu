
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <cuda_fp16.h>
#include <iostream>
#include <ctime>
cudaError_t addWithCuda(int* c, const int* a, const int* b, unsigned int size);

__global__ void addKernel(int* c, const int* a, const int* b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

using float16_t = half;

const int manLength = 10;
const int expLength = 5;
const int signLength = 1;
const int shiftExp = (1 << (expLength - 1)) - 1;

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

class FP16 {
private:
    unsigned int man : manLength;
    unsigned int exp : expLength;
    unsigned int sign : signLength;

public:
    FP16(int _sign = 0, int _exp = 0, int _man = 0) :sign(_sign), exp(_exp), man(_man) {}

    FP16(float16_t f) {
        unsigned short bits = float16_bites_to_short(f);
        sign = (bits >> 15) & 0x1;
        exp = (bits >> manLength) & ((1 << expLength) - 1);
        man = bits & ((1 << manLength) - 1);
    }

    bool IsSubnormal() { return (exp == 0) && (man != 0); }
    bool IsInf() { return exp == ((1 << expLength) - 1) && (man == 0); }
    bool IsNull() { return (exp == 0) && (man == 0); }
    bool IsNan() { return ((exp == ((1 << expLength) - 1)) && (man != 0)); }
    void printFP16_bites() {
        std::cout << sign << "_";
        for (char i = 0; i < expLength; ++i) std::cout << ((exp >> (expLength - i - 1)) & 1);
        std::cout << "_";
        for (char i = 0; i < manLength; ++i) std::cout << ((man >> (manLength - i - 1)) & 1);
    }

    operator float16_t() const {
        unsigned short bits = (sign << 15) | (exp << manLength) | man;
        return *reinterpret_cast<float16_t*>(&bits);
    }

    FP16& operator=(FP16& N) {
        if (this != &N) {
            sign = N.sign;
            exp = N.exp;
            man = N.man;
        }
        return *this;
    }

    uint16_t get_int() {
        uint16_t temp = (sign << (manLength + expLength)) + (exp << manLength) + man;
        return temp;
    }

    FP16 operator+(FP16 right) {
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

        if (temp < (1 << manLength) && max_exp == 1) {
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

        if (temp < (1 << manLength) && (max_exp == 1)) {
            temp += (1 << manLength);
            max_exp = 0;
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

    FP16 operator-(FP16 right) {
        right.sign += 1;
        return *this + right;
    }

    FP16 operator*(FP16 right) {
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

    friend FP16 fma(FP16 a, FP16 b, FP16 c);
};

//a*b + c
FP16 fma(FP16 a, FP16 b, FP16 c) {
    FP16 Res;
    uint32_t temp_man;


    return Res;
}



void print_bites(FP16 a) { a.printFP16_bites(); }
void print_bites(float16_t a) { print_float16_bites(a); }

void Test_Add();
void Test_Sub();
void Test_Mul();
void TimeTest();

using namespace std;
int main() {
    FP16 A(0, 0, 1);
    FP16 B(1, 17, 1);

    float16_t a = A;
    float16_t b = B;
    float16_t c;
    cout << "A: "; print_bites(A); cout << ", " << __half2float(A) << endl;
    cout << "B: "; print_bites(B); cout << ", " << __half2float(B) << endl;
    cout << "A + B: "; print_bites(A + B); cout << " - My" << ", " << __half2float(A + B) << endl;
    cout << "a + b: "; print_bites(a + b); cout << ", " << __half2float(a + b) << endl;
    cout << "A - B: "; print_bites(A - B); cout << " - My" << ", " << __half2float(A - B) << endl;
    cout << "a - b: "; print_bites(a - b); cout << ", " << __half2float(a - b) << endl;
    cout << "A * B: "; print_bites(A * B); cout << ", " << __half2float(A * B) << " - My" << endl;
    cout << "a * b: "; print_bites(a * b); cout << ", " << __half2float(a * b) << endl;
    
    Test_Mul();
    

    return 0;
}


void Test_Add() {
    cout << "Addition test:" << endl;
    float16_t a, b, c;
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength) - 1; ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength) - 1; ++exp2) {
                        for (size_t man2 = 0; man2 < (1 << manLength); ++man2) {
                            FP16 A(sign1, exp1, man1);
                            FP16 B(sign2, exp2, man2);
                            FP16 C = A + B;
                            a = A; b = B;
                            c = a + b;
                            if (abs(float16_bites_to_short(C) - float16_bites_to_short(c)) >= 1) {
                                cout << "A: ";
                                print_bites(A);
                                cout << endl;
                                cout << "B: ";
                                print_bites(B);
                                cout << endl;
                                cout << "A + B: ";
                                print_bites(A + B);
                                cout << " - My " << endl;
                                cout << "a + b: ";
                                print_bites(a + b);
                                cout << endl;
                                cout << sign1 << " " << exp1 << " " << man1 << " " << endl;
                                cout << sign2 << " " << exp2 << " " << man2 << " " << endl;
                                //return;
                            }
                        }
                    }
                }
            }
            cout << "Exp = " << exp1 << endl;
        }
    }
}

void Test_Sub() {
    cout << "Subrtact test:" << endl;
    float16_t a, b, c;
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength) - 1; ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength) - 1; ++exp2) {
                        for (size_t man2 = 0; man2 < (1 << manLength); ++man2) {
                            FP16 A(sign1, exp1, man1);
                            FP16 B(sign2, exp2, man2);
                            FP16 C = A - B;
                            a = A; b = B;
                            c = a - b;
                            if (abs(float16_bites_to_short(C) - float16_bites_to_short(c)) >= 1) {
                                cout << "A: ";
                                print_bites(A);
                                cout << endl;
                                cout << "B: ";
                                print_bites(B);
                                cout << endl;
                                cout << "A - B: ";
                                print_bites(A - B);
                                cout << " - My " << endl;
                                cout << "a - b: ";
                                print_bites(a - b);
                                cout << endl;
                                cout << sign1 << " " << exp1 << " " << man1 << " " << endl;
                                cout << sign2 << " " << exp2 << " " << man2 << " " << endl;
                                //return;

                            }
                        }
                    }
                }
            }
            cout << "Exp = " << exp1 << endl;
        }
    }
}


void Test_Mul() {
    cout << "Mult test:" << endl;
    float16_t a, b, c;
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength) - 1; ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength) - 1; ++exp2) {
                        for (size_t man2 = 0; man2 < (1 << manLength); ++man2) {
                            FP16 A(sign1, exp1, man1);
                            FP16 B(sign2, exp2, man2);
                            FP16 C = A * B;
                            a = A; b = B;
                            c = a * b;
                            if (abs(float16_bites_to_short(C) - float16_bites_to_short(c)) >= 1) {
                                cout << "A: ";
                                print_bites(A);
                                cout << endl;
                                cout << "B: ";
                                print_bites(B);
                                cout << endl;
                                cout << "A * B: ";
                                print_bites(A * B);
                                cout << " - My " << endl;
                                cout << "a * b: ";
                                print_bites(a * b);
                                cout << endl;
                                cout << sign1 << " " << exp1 << " " << man1 << " " << endl;
                                cout << sign2 << " " << exp2 << " " << man2 << " " << endl;
                                //return;                                
                            }
                        }
                    }
                }
            }
            cout << "Exp = " << exp1 << endl;
        }
    }
}

void TimeTest() {
    FP16 A;
    FP16 B;
    FP16 C;
    int start,end,t;
    cout << endl << "Addition Time Test:" << endl;
    start = clock();
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength) - 1; ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength) - 1; ++exp2) {
                        for (size_t man2 = 0; man2 < (1 << manLength); ++man2) {
                            A = FP16(sign1, exp1, man1);
                            B = FP16(sign2, exp2, man2);
                            C = A + B;
                        }
                    }
                }
            }
        }
    }
    end = clock();
    t = (end - start) / CLOCKS_PER_SEC;

    cout << "Time: " << t << " seconds" << endl;

    cout << endl << "Subtract Time Test:" << endl;

    start = clock();
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength) - 1; ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength) - 1; ++exp2) {
                        for (size_t man2 = 0; man2 < (1 << manLength); ++man2) {
                            A = FP16(sign1, exp1, man1);
                            B = FP16(sign2, exp2, man2);
                            C = A - B;
                        }
                    }
                }
            }
        }
    }
    end = clock();
    t = (end - start) / CLOCKS_PER_SEC;

    cout << "Time: " << t << " seconds" << endl;

    cout << endl << "Multiplication Time Test:" << endl;

    start = clock();
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength) - 1; ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength) - 1; ++exp2) {
                        for (size_t man2 = 0; man2 < (1 << manLength); ++man2) {
                            A = FP16(sign1, exp1, man1);
                            B = FP16(sign2, exp2, man2);
                            C = A * B;
                        }
                    }
                }
            }
        }
    }
    end = clock();
    t = (end - start) / CLOCKS_PER_SEC;

    cout << "Time: " << t << " seconds" << endl;

}


// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int* c, const int* a, const int* b, unsigned int size)
{
    int* dev_a = 0;
    int* dev_b = 0;
    int* dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel << <1, size >> > (dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);

    return cudaStatus;
}
