#pragma once
// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// Copyright (c) Сысоев А.В.
//
//

#ifndef __TDynamicMatrix_H__
#define __TDynamicMatrix_H__
#define _CRT_SECURE_NO_WARNINGS

//#include <cmath>
#include <vector>
#include <math.h>
#include <iostream>
#include <random>
#include <iomanip>
#include <stdio.h>
#include <omp.h>
//#include <boost/numeric/ublas/blas.hpp>

using namespace std;

const long int MAX_VECTOR_SIZE = 32786 * 32786;
const int MAX_MATRIX_SIZE = 32786;

// Динамический вектор - 
// шаблонный вектор на динамической памяти
template<typename T>
class TDynamicVector
{
public:
    long long sz;
    T* pMem;
    TDynamicVector(size_t size = 1) : sz(size)
    {
        pMem = new T[sz]();// {}; // У типа T д.б. констуктор по умолчанию
    }
    TDynamicVector(T* arr, size_t s) : sz(s)
    {
        if (sz == 0)
            throw out_of_range("Vector size should be greater than zero");
        if (sz > MAX_VECTOR_SIZE)
            throw out_of_range("Too large vector");
        if (arr == nullptr)
            throw invalid_argument("TDynamicVector ctor requires non-nullptr arg");
        pMem = new T[sz]();
//        memcpy(pMem, arr, sz + sizeof(T));
        std::copy(arr, arr + sz, pMem);
    }
    TDynamicVector(const TDynamicVector& v)
    {
        sz = v.sz;
        pMem = new T[sz]();
        std::copy(v.pMem, v.pMem + sz, pMem);

    }
    TDynamicVector(TDynamicVector&& v) noexcept
    {
        pMem = nullptr;
        swap(*this, v);
    }
    ~TDynamicVector()
    {
        delete[] pMem;
    }
    TDynamicVector& operator=(const TDynamicVector& v)
    {
        if (this == &v) return *this;
        if (sz != v.sz) {
            T* p = new T[v.sz]();
            delete[] pMem;
            sz = v.sz;
            pMem = p;
        }
        std::copy(v.pMem, v.pMem + sz, pMem);
        return *this;
    }
    TDynamicVector& operator=(TDynamicVector&& v) noexcept
    {
        swap(*this, v);
        return *this;
    }

    long long size() const noexcept { return sz; }

    // индексация
    T& operator[](size_t ind)
    {
        return pMem[ind];
    }
    const T& operator[](size_t ind) const
    {
        return pMem[ind];
    }
    // индексация с контролем
    T& at(size_t ind)
    {
        if (ind >= sz) throw out_of_range("Out of range");
        return pMem[ind];
    }
    const T& at(size_t ind) const
    {
        if (ind >= sz) throw out_of_range("Out of range");
        return pMem[ind];
    }

    // сравнение
    bool operator==(const TDynamicVector& v) const noexcept
    {
        if (sz != v.sz) return false;
        for (size_t i = 0; i < sz; ++i)
            if (pMem[i] != v.pMem[i])
                return false;

        return true;
    }
    bool operator!=(const TDynamicVector& v) const noexcept
    {
        return !(*this == v);
    }

    // скалярные операции
    TDynamicVector operator+(T val)
    {
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] + val;
        return tmp;
    }
    TDynamicVector operator-(T val)
    {
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] - val;
        return tmp;
    }
    TDynamicVector operator*(T val)
    {
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] * val;
        return tmp;
    }

    // векторные операции
    TDynamicVector operator+(const TDynamicVector& v) const
    {
        if (sz != v.sz) throw length_error("Incompatible sizes");
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] + v.pMem[i];
        return tmp;
    }
    TDynamicVector operator-(const TDynamicVector& v) const
    {
        if (sz != v.sz) throw length_error("Incompatible sizes");
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] - v.pMem[i];
        return tmp;
    }
    T operator*(const TDynamicVector& v) const
    {
        if (sz != v.sz) throw length_error("Incompatible sizes");
        T tmp = T();
        for (size_t i = 0; i < sz; i++)
            tmp += pMem[i] * v.pMem[i];
        return tmp;
    }

    friend void swap(TDynamicVector& lhs, TDynamicVector& rhs) noexcept
    {
        std::swap(lhs.sz, rhs.sz);
        std::swap(lhs.pMem, rhs.pMem);
    }

    // ввод/вывод
    friend istream& operator>>(istream& istr, TDynamicVector& v)
    {
        for (size_t i = 0; i < v.sz; i++)
            istr >> v.pMem[i]; // требуется оператор>> для типа T
        return istr;
    }
    friend ostream& operator<<(ostream& ostr, const TDynamicVector& v)
    {
        for (size_t i = 0; i < v.sz; i++)
            ostr << std::setw(15) << v.pMem[i]; // требуется оператор<< для типа T
        return ostr;
    }

    double norm1() {
        double res = 0.0;
        for (size_t i = 0; i < sz; ++i)
            res += abs(double(pMem[i]));
        return res;
    }

    double norm2() {
        double res = 0.0;
        for (size_t i = 0; i < sz; ++i)
            res += double(pMem[i]) * double(pMem[i]);
        return std::sqrt(res);
    }

    void generate() {
        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_real_distribution<double> coef_gen(-1.0, 1.0);

        for (size_t i = 0; i < sz; ++i) {
            pMem[i] = coef_gen(e);
        }
    }
};


// Динамическая матрица - 
// шаблонная матрица на динамической памяти
template<typename T>
class TDynamicMatrix
{
public:
    T* pMem;
    long long sz;
    TDynamicMatrix(size_t s = 1)
    {
        pMem = new T[s * s]();
        sz = s;
    }
    ~TDynamicMatrix() {
        delete[] pMem;
    }
    TDynamicMatrix(const TDynamicMatrix& m) {
        pMem = new T[m.sz * m.sz]();
        sz = m.sz;
        std::copy(m.pMem, m.pMem + sz * sz, pMem);
    }
    TDynamicMatrix(TDynamicMatrix&& m) {
        pMem = m.pMem;
        sz = m.sz;
        m.pMem = nullptr;
    }


    TDynamicMatrix<T>& operator=(const TDynamicMatrix& m) {
        if (this == &m) return *this;
        if (sz != m.sz) {
            T* p = new T[m.sz * m.sz]();
            delete[] pMem;
            sz = m.sz;
            pMem = p;
        }
        std::copy(m.pMem, m.pMem + sz * sz, pMem);
        return *this;
    }
    TDynamicMatrix<T>& operator=(TDynamicMatrix&& m) {
        if (this != &m) {
            std::swap(pMem, m.pMem);
            sz = m.sz;
        }
        return *this;
    }

    // сравнение
    bool operator==(const TDynamicMatrix& m) const noexcept
    {
        if (m.sz != sz) return false;
        for (size_t i = 0; i < sz * sz; ++i)
            if (pMem[i] != m.pMem[i])
                return false;
        return true;
    }

    bool operator!=(const TDynamicMatrix& m) const noexcept // my func
    {
        return !(*this == m);
    }

    // матрично-скалярные операции
    TDynamicMatrix<T> operator*(const T& val)
    {
        TDynamicMatrix<T> tmp(sz);
        for (size_t i = 0; i < sz * sz; ++i)
            tmp.pMem[i] = pMem[i] * val;
        return tmp;
    }

    long long size() const noexcept {
        return sz;
    }

    T& operator()(size_t i, size_t j) noexcept {
        return pMem[j * sz + i];
    }

    const T& operator()(size_t i, size_t j) const noexcept {
        return pMem[j * sz + i];
    }

    // матрично-векторные операции
    TDynamicVector<T> operator*(const TDynamicVector<T>& v)
    {
        if (sz != v.size()) throw logic_error("ERR");
        T sum;
        TDynamicVector<T> tmp(sz);
        for (size_t i = 0; i < sz; ++i) {
            sum = 0.0;
            for (size_t j = 0; j < sz; ++j) {
                sum += pMem[j * sz + i] * v[j];
            }
            tmp[i] = sum;
        }
        return tmp;
    }

    // матрично-матричные операции
    TDynamicMatrix operator+(const TDynamicMatrix& m)
    {
        if (sz != m.sz) throw logic_error("ERR");
        TDynamicMatrix tmp(sz);
        for (size_t i = 0; i < sz * sz; i++)
            tmp.pMem[i] = pMem[i] + m.pMem[i];
        return tmp;
    }

    TDynamicMatrix operator-(const TDynamicMatrix& m)
    {
        if (sz != m.sz) throw logic_error("ERR");
        TDynamicMatrix tmp(sz);
        for (size_t i = 0; i < sz * sz; i++)
            tmp.pMem[i] = pMem[i] - m.pMem[i];
        return tmp;
    }
    TDynamicMatrix operator*(const TDynamicMatrix& m) //square matrixes
    {
        if (sz != m.sz) throw length_error("incompatible sizes");
        TDynamicMatrix res(sz);

        for (long long k = 0; k < sz; ++k) {
#pragma omp parallel for
            for (long long i = 0; i < sz; ++i) {
                for (long long j = 0; j < sz; ++j) {
                    res(i, j) += (*this)(i, k) * m(k, j);
                }
            }
        }

        return res;
    }

    double norm1() const
    {
        double res = 0;
        double sum;
        for (size_t j = 0; j < sz; ++j) {
            sum = 0;
            for (size_t i = 0; i < sz; ++i)
                sum += std::abs(double((*this)(i, j)));
            if (sum > res) res = sum;
        }
        return res;
    }

    double norminf() const
    {
        double res = 0;
        double sum;
        for (size_t i = 0; i < sz; ++i) {
            sum = 0;
            for (size_t j = 0; j < sz; ++j)
                sum += std::abs(double((*this)(i, j)));
            if (sum > res) res = sum;
        }
        return res;
    }

    double normmax() const {
        double res = 0;
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                res = std::max(double(std::abs((*this)(i, j))), res);
            }
        }
        return res;
    }

    double normmin() const
    {
        double res = 0;
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                res = std::min(double(std::abs((*this)(i, j))), res);
            }
        }
        return res;
    }

    void LUdecomposition(TDynamicMatrix& L, TDynamicMatrix& U) const {
        L = TDynamicMatrix<T>(sz);
        U = TDynamicMatrix<T>(sz);
        for (size_t i = 0; i < sz; ++i) {
            L(i, i) = 1.0;
        }

        T sum;
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                sum = 0;
                if (i <= j) {
                    for (size_t k = 0; k <= i; ++k) {
                       sum += L(i, k) * U(k, j);
                    }
                    U(i, j) = (*this)(i, j) - sum;
                }
                else {
                    for (size_t k = 0; k <= j; ++k) {
                        sum += L(i, k) * U(k, j);
                    }
                    L(i, j) = ((*this)(i, j) - sum) / U(j, j);

                }
            }
        }
    }

    void LUdecompositionOptimized(TDynamicMatrix& L, TDynamicMatrix& U) const {
        L = TDynamicMatrix<T>(sz);
        U = TDynamicMatrix<T>(sz);
        for (size_t i = 0; i < sz; ++i) {
            L.pMem[i * sz + i] = 1.0;
        }

        T sum;
        for (size_t i = 0; i < sz; ++i) {

#pragma omp parallel for
            for (long long j = 0; j < sz; ++j) {
                sum = 0;
                if (i <= j) {

//#pragma omp simd for
                   for (size_t k = 0; k <= i; ++k) {
                        sum += L.pMem[k * sz + i] * U.pMem[j * sz + k];
                   }

                    U.pMem[j * sz + i] = pMem[j * sz + i] - sum;
                }
                else {

//#pragma omp simd for
                    for (size_t k = 0; k <= j; ++k) {
                        sum += L.pMem[k * sz + i] * U.pMem[j * sz + k];
                    }

                    L.pMem[j * sz + i] = (pMem[j * sz + i] - sum) / U.pMem[j * sz + j];
                }
            }
        }
    }

    // ввод/вывод
    friend istream& operator>>(istream& istr, TDynamicMatrix& m)
    {
        for (size_t i = 0; i < m.sz; ++i) {
            for (size_t j = 0; j < m.sz; ++j)
                istr >> m(i, j);
        }
        return istr;
    }
    friend ostream& operator<<(ostream& ostr, const TDynamicMatrix& m)
    {
        for (size_t i = 0; i < m.sz; ++i) {
            for (size_t j = 0; j < m.sz; ++j)
                ostr << std::setw(15) << double(m(i, j));
            std::cout << std::endl;
        }
        return ostr;
    }
    void simpleScan(const vector<double>& vec) {
        for (long long i = 0; i < sz * sz; ++i) pMem[i] = vec[i];
    }

    void LUdecompress(TDynamicMatrix& L, TDynamicMatrix& U) const {
        L = TDynamicMatrix<T>(sz);
        U = TDynamicMatrix<T>(sz);
        for (size_t i = 0; i < sz; ++i) {
            L(i, i) = 1;
        }
        for (size_t j = 0; j < sz; ++j) {
            for (size_t i = j + 1; i < sz; ++i) {
                L(i, j) = (*this)(i, j);
            }
        }
        for (size_t j = 0; j < sz; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                U(i, j) = (*this)(i, j);
            }
        }
    }

    void LUdecompositionV3(TDynamicMatrix& L, TDynamicMatrix& U) const {
        long long stepr = 2; // 64 for 80Kbyte L1 cash 
        T sum;
        long long i, j, k, r, nr;

        // MATRIX DEFENITION
        TDynamicMatrix<T> A = *this;
        for (i = 0; i < sz * sz; ++i) {
            L.pMem[i] = 0.0;
            U.pMem[i] = 0.0;
        }
        for (i = 0; i < sz; ++i) {
            L.pMem[i * sz + i] = 1.0;
        }

        // MAIN CYCLE (as recursive algorithm)
        for (r = 0; r < sz; r += stepr) {
//            cout << "r: " << r << endl << A << endl << L << endl << U << endl;
//            cout << "r: " << r << endl;
            nr = sz - r - stepr;
            if (nr < 0) break;

            // Gauss for L11 and U11
            for (i = 0; i < stepr; ++i) {
//#pragma omp parallel for
                for (j = 0; j < stepr; ++j) {
                    sum = 0;
                    if (i <= j) {
                        //#pragma omp simd for
                        for (k = 0; k <= i; ++k) {
                            sum += L.pMem[(k + r) * sz + i + r] * U.pMem[(j + r) * sz + k + r];
                        }

                        U.pMem[(j + r) * sz + i + r] = A.pMem[(j + r) * sz + i + r] - sum;
                    }
                    else {
                        //#pragma omp simd for
                        for (k = 0; k <= j; ++k) {
                            sum += L.pMem[(k + r) * sz + i + r] * U.pMem[(j + r) * sz + k + r];
                        }

                        L.pMem[(j + r) * sz + i + r] = (A.pMem[(j + r) * sz + i + r] - sum) / U.pMem[(j + r) * sz + j + r];
                    }
                }
            }

            // SOLVING NR OF TRIANGLE SYSTEMS
            for (j = 0; j < nr; ++j) {
                // L - lower triangle system
//#pragma omp parallel for
                for (i = 0; i < stepr; ++i) {
                    U.pMem[(j + r + stepr) * sz + i + r] = A.pMem[(j + r + stepr) * sz + i + r];
                    for (k = i - 1; k >= 0; --k) {
                        U.pMem[(j + r + stepr) * sz + i + r] -= U.pMem[(j + r + stepr) * sz + k + r] * L.pMem[(k + r) * sz + i + r];
                    }
                }
                // U - upper triangle system
//#pragma omp parallel for
                for (i = stepr - 1; i >= 0; --i) {
                    L.pMem[(i + r) * sz + j + r + stepr] = A.pMem[(i + r) * sz + j + r + stepr];
                    for (k = i + 1; k < stepr; ++k) {
                        L.pMem[(i + r) * sz + j + r + stepr] -= L.pMem[(k + r) * sz + j + r + stepr] * U.pMem[(i + r) * sz + k + r];
                    }
                    L.pMem[(i + r) * sz + j + r + stepr] = L.pMem[(i + r) * sz + j + r + stepr] / U.pMem[(i + r) * sz + i + r];
                }
            }

            // REDUCTION matrix computing
//#pragma omp parallel
//#pragma omp parallel
            for (j = 0; j < nr; ++j) {
//#pragma omp for
                for (i = 0; i < nr; ++i) {
#pragma omp simd
                    for (k = 0; k < stepr; ++k) {
                        A.pMem[(j + r + stepr) * sz + i + r + stepr] -= L.pMem[(k + r) * sz + i + r + stepr] * U.pMem[(j + r + stepr) * sz + k + r];
                    }
                }
            }
        }

        r = sz % stepr;
        nr = sz - r;
        // Gauss for a bottom-left matrix
        for (i = 0; i < r; ++i) {
//#pragma omp parallel for // CREATES AN ERROR WHEN TURNED ON!!!!!!!!!!!!!!! WHY????????????????
            for (j = 0; j < r; ++j) {
                sum = 0;
                if (i <= j) {
                    //#pragma omp simd for
                    for (k = 0; k <= i; ++k) {
                        sum += L.pMem[(k + nr) * sz + i + nr] * U.pMem[(j + nr) * sz + k + nr];
                    }

                    U.pMem[(j + nr) * sz + i + nr] = A.pMem[(j + nr) * sz + i + nr] - sum;
                }
                else {
                    //#pragma omp simd for
                    for (k = 0; k <= j; ++k) {
                        sum += L.pMem[(k + nr) * sz + i + nr] * U.pMem[(j + nr) * sz + k + nr];
                    }

                    L.pMem[(j + nr) * sz + i + nr] = (A.pMem[(j + nr) * sz + i + nr] - sum) / U.pMem[(j + nr) * sz + j + nr];
                }
            }
        }
    }

    TDynamicVector<T> solver(const TDynamicVector<T>& b, int version) const {
        TDynamicMatrix<T> A(sz);
        A = *this;
        if (version == 2) LUdecompositionV2(A);
        TDynamicMatrix<T> L(sz), U(sz);
        if (version == 0) LUdecomposition(L, U);
        if (version == 1) LUdecompositionOptimized(L, U);
        if (version == 2) A.LUdecompress(L, U); // can reduce memory usage!
        if (version == 3) A.LUdecompositionV3(L, U);
        TDynamicVector<T> y(sz);
        // Ly = b, L - upper triangle
        for (size_t i = 0; i < sz; ++i) {
            y[i] = b[i];
            for (size_t j = i; j > 0; --j) {
                y[i] -= y[j - 1] * L(i, j - 1);
            }
        }
        TDynamicVector<T> x(sz);
        // Ux = y, U - bottom triangle
        for (size_t i = sz; i > 0; --i) {
            x[i - 1] = y[i - 1];
            for (size_t j = i; j < sz; ++j) {
                x[i - 1] -= x[j] * U(i - 1, j);
            }
            x[i - 1] = x[i - 1] / U(i - 1, i - 1);
        }
        return x;
    }

    TDynamicVector<T> solver(const TDynamicMatrix<T>& L, const TDynamicMatrix<T>& U, const TDynamicVector<T>& b) const {
        TDynamicVector<T> y(sz);
        // Ly = b, L - upper triangle
        for (size_t i = 0; i < sz; ++i) {
            y[i] = b[i];
            for (size_t j = i; j > 0; --j) {
                y[i] -= y[j - 1] * L(i, j - 1);
            }
        }
        TDynamicVector<T> x(sz);
        // Ux = y, U - bottom triangle
        for (size_t i = sz; i > 0; --i) {
            x[i - 1] = y[i - 1];
            for (size_t j = i; j < sz; ++j) {
                x[i - 1] -= x[j] * U(i - 1, j);
            }
            x[i - 1] = x[i - 1] / U(i - 1, i - 1);
        }
        return x;
    }

    void generate() {
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                (*this)(i, j) = 0.0;
            }
            (*this)(i, i) = 1.0;
        }

        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_int_distribution<size_t> index_gen(0, sz - 1);
        std::uniform_real_distribution<double> coef_gen(-10.0 / (sz * sz), 10.0 / (sz * sz));
        size_t count = sz * sz;
        size_t indexI, indexJ;
        double coef;

#pragma omp parallel for
        for (long long _ = 0; _ < count; ++_) {
            indexI = index_gen(e);
            indexJ = index_gen(e);
            coef = coef_gen(e);

            for (size_t i = 0; i < sz; ++i) {
                pMem[indexI * sz + i] += pMem[indexJ * sz + i] * T(coef);
            }
        }
    }

    void simple_generate() {
        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_real_distribution<double> coef_gen(-1.0, 1.0);
        double coef;

        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                (*this)(i, j) = coef_gen(e);
            }
        }
    }

    void writeMatrix(FILE* file) { // max 32786 (8gb), 16384 (2gb) i think 10000 is suitable

        fwrite(pMem, sizeof(double), sz * sz, file);
    }

    void readDMatrix(FILE* file, size_t readsize, size_t offset) {
//        if (offset + sz > readsize) return;
        fread(pMem, sizeof(double), sz * sz, file);

        
//        double* input = (double*)malloc(sizeof(double) * readsize * readsize);
//        fread(input, sizeof(double), readsize * readsize, file);
//        for (size_t i = 0; i < sz; ++i) {
//            memcpy(pMem, input + (i * readsize + offset) * sizeof(double), sz * sizeof(double));
//        }
//        free(input);
    }
};

template <typename T>
void printSimpleMatrix(T* A, size_t m, size_t n) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << std::setw(15) << double(A[j * m + i]);
        }
        std::cout << std::endl;
    }
}

template <typename T>
void LUdecompositionV2(TDynamicMatrix<T>& A) {
    size_t r = 128ull;
    size_t nr = A.sz - r;
    long long i, j, k;
    
    // CHECK ON SIZE
    if (A.sz <= r) {
        r = A.sz;
        T* L = (T*)malloc(sizeof(T) * r * r);
        T* U = (T*)malloc(sizeof(T) * r * r);

        // GAUSS METHOD FOR L U
        for (i = 0; i < r * r; ++i) {
            L[i] = 0.0;
            U[i] = 0.0;
        }
        for (i = 0; i < r; ++i) {
            L[i * r + i] = 1.0;
        }
        T sum;
        for (i = 0; i < r; ++i) {
#pragma omp parallel for
            for (j = 0; j < r; ++j) {
                sum = 0;
                if (i <= j) {
                    //#pragma omp simd for
                    for (k = 0; k <= i; ++k) {
                        sum += L[k * r + i] * U[j * r + k];
                    }

                    U[j * r + i] = A.pMem[j * A.sz + i] - sum;
                }
                else {
                    //#pragma omp simd for
                    for (k = 0; k <= j; ++k) {
                        sum += L[k * r + i] * U[j * r + k];
                    }

                    L[j * r + i] = (A.pMem[j * A.sz + i] - sum) / U[j * r + j];
                }
            }
        }

        // Compress
        for (j = 0; j < r; ++j) {
            for (i = j + 1; i < r; ++i) {
                A.pMem[j * A.sz + i] = L[r * j + i];
            }
        }
        for (j = 0; j < r; ++j) {
            for (i = 0; i <= j; ++i) {
                A.pMem[j * A.sz + i] = U[r * j + i];
            }
        }

        free(L);
        free(U);
        return;
    }

    // ALLOCATE MEMORY
    T* L11 = (T*)malloc(sizeof(T) * r * r);
    T* U11 = (T*)malloc(sizeof(T) * r * r);
    T* L21 = (T*)malloc(sizeof(T) * r * nr);
    T* U12 = (T*)malloc(sizeof(T) * nr * r);
//    T* L22 = (T*)malloc(sizeof(T) * nr * nr);
//    T* U22 = (T*)malloc(sizeof(T) * nr * nr);
    TDynamicMatrix<T> A22(nr);

    // GAUSS METHOD FOR L11 U11
    for (i = 0; i < r * r; ++i) {
        L11[i] = 0.0;
        U11[i] = 0.0;
    }
    for (i = 0; i < r; ++i) {
        L11[i * r + i] = 1.0;
    }
    T sum;
    for (i = 0; i < r; ++i) {
#pragma omp parallel for
        for (j = 0; j < r; ++j) {
            sum = 0;
            if (i <= j) {
                //#pragma omp simd for
                for (k = 0; k <= i; ++k) {
                    sum += L11[k * r + i] * U11[j * r + k];
                }

                U11[j * r + i] = A.pMem[j * A.sz + i] - sum;
            }
            else {
                //#pragma omp simd for
                for (k = 0; k <= j; ++k) {
                    sum += L11[k * r + i] * U11[j * r + k];
                }

                L11[j * r + i] = (A.pMem[j * A.sz + i] - sum) / U11[j * r + j];
            }
        }
    }

//    cout << "L11 and U11" << endl;
//    printSimpleMatrix(L11, r, r);
//    cout << endl;
//    printSimpleMatrix(U11, r, r);
//    cout << endl;

    // SOLVING NR OF TRIANGLE SYSTEMS
    for (j = 0; j < nr; ++j) {
        // L - lower triangle system
        for (i = 0; i < r; ++i) {
            U12[j * r + i] = A.pMem[(j + r) * A.sz + i];
            for (k = i - 1; k >= 0; --k) {
                U12[j * r + i] -= U12[j * r + k] * L11[k * r + i];
            }
        }

        // U - upper triangle system
        for (i = r - 1; i >= 0; --i) {
            L21[i * r + j] = A.pMem[i * A.sz + j + r];
            for (k = i + 1; k < r; ++k) {
                L21[i * r + j] -= L21[k * r + j] * U11[i * r + k];
            }
            L21[i * r + j] = L21[i * r + j] / U11[i * r + i];
        }
    }

//    cout << "L21 and U12" << endl;
//    printSimpleMatrix(L21, nr, r);
//    cout << endl;
//    printSimpleMatrix(U12, r, nr);
//    cout << endl;

    // A22 - reduced matrix calculating
    // GEMM usage required!
    for (i = 0; i < nr; ++i) {
        for (k = 0; k < r; ++k) {
            for (j = 0; j < nr; ++j) {
                A22.pMem[j * nr + i] += L21[k * nr + i] * U12[j * r + k];
            }
        }
    }
//    cout << "L21*U12" << endl;
//    printSimpleMatrix(A22.pMem, nr, nr);
//    cout << endl;
    for (j = 0; j < nr; ++j) {
        for (i = 0; i < nr; ++i) {
            A22.pMem[j * nr + i] = A.pMem[(r + j) * A.sz + i + r] - A22.pMem[j * nr + i];
        }
    }
//    cout << "A22 reducted" << endl;
//    printSimpleMatrix(A22.pMem, nr, nr);
//    cout << endl;

    // RECURSIVE COMPUTATION OF A22 = L22*U22
    //...
    LUdecompositionV2(A22);

    // LUcompress
    for (j = 0; j < A.sz; ++j) {
        for (i = 0; i < A.sz; ++i) {
            A.pMem[j * A.sz + i] = 0.0;
        }
    }
    for (j = 0; j < r; ++j) {
        for (i = j + 1; i < r; ++i) {
            A.pMem[j * A.sz + i] = L11[r * j + i];
        }
    }
    for (j = 0; j < r; ++j) {
        for (i = 0; i <= j; ++i) {
            A.pMem[j * A.sz + i] = U11[r * j + i];
        }
    }
    for (j = 0; j < nr; ++j) {
        for (i = 0; i < nr; ++i) {
            A.pMem[(r + j) * A.sz + i + r] = A22.pMem[j * nr + i];
        }
    }
    for (j = 0; j < r; ++j) {
        for (i = 0; i < nr; ++i) {
            A.pMem[j * A.sz + i + r] = L21[nr * j + i];
        }
    }
    for (j = 0; j < nr; ++j) {
        for (i = 0; i < r; ++i) {
            A.pMem[(r + j) * A.sz + i] = U12[r * j + i];
        }
    }
//    cout << "A result" << endl;
//    printSimpleMatrix(A.pMem, A.sz, A.sz);
//    cout << endl;

    free(L11);
    free(U11);
    free(L21);
    free(U12);
//    free(L22);
//    free(U22);
}

#endif
