// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// Copyright (c) Сысоев А.В.
//
//

#ifndef __TDynamicMatrix_H__
#define __TDynamicMatrix_H__

#include <iostream>
#include <cassert>
#include <random>
#include <iomanip>
#include <vector>
using namespace std;

const int MAX_VECTOR_SIZE = 100000000;
const int MAX_MATRIX_SIZE = 10000;

// Динамический вектор - 
// шаблонный вектор на динамической памяти
template<typename T>
class TDynamicVector
{
protected:
    size_t sz;
    //    T* pMem;
    vector<T> pMem;
public:
    TDynamicVector(size_t size = 1) : sz(size)
    {
        if (sz == 0)
            throw std::out_of_range("Vector size should be greater than zero");
        if (sz > MAX_VECTOR_SIZE)
            throw std::out_of_range("Too large vector");
        //        pMem = new T[sz]();// {}; // У типа T д.б. констуктор по умолчанию
        pMem.resize(sz);
    }
    TDynamicVector(T* arr, size_t s) : sz(s)
    {
        if (sz == 0)
            throw std::out_of_range("Vector size should be greater than zero");
        if (sz > MAX_VECTOR_SIZE)
            throw std::out_of_range("Too large vector");
        if (arr == nullptr)
            throw std::invalid_argument("TDynamicVector ctor requires non-nullptr arg");
        //        pMem = new T[sz];
        pMem.resize(sz);
        std::copy(arr, arr + sz, pMem);
    }
    TDynamicVector(const TDynamicVector& v)
    {
        sz = v.sz;
        //        pMem = new T[sz];
        pMem.resize(sz);
        std::copy(v.pMem.begin(), v.pMem.end(), pMem.begin());
    }
    TDynamicVector(TDynamicVector&& v) noexcept
    {
        //        pMem = nullptr;
        swap(*this, v);
    }
    ~TDynamicVector()
    {
        //        delete[] pMem;
    }
    TDynamicVector& operator=(const TDynamicVector& v)
    {
        if (this == &v) return *this;
        //        if (sz != v.sz) {
        //            T* p = new T[v.sz];
        //            delete[] pMem;
        //            sz = v.sz;
        //            pMem = p;
        //        }
        //        std::copy(v.pMem, v.pMem + sz, pMem);
        pMem = v.pMem;
        sz = v.sz;
        return *this;
    }
    TDynamicVector& operator=(TDynamicVector&& v) noexcept
    {
        swap(*this, v);
        return *this;
    }

    size_t size() const noexcept { return sz; }

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
        if (ind >= sz) throw std::out_of_range("Out of range");
        return pMem[ind];
    }
    const T& at(size_t ind) const
    {
        if (ind >= sz) throw std::out_of_range("Out of range");
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
        if (sz != v.sz) throw std::length_error("Incompatible sizes");
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] + v.pMem[i];
        return tmp;
    }
    TDynamicVector operator-(const TDynamicVector& v) const
    {
        if (sz != v.sz) throw std::length_error("Incompatible sizes");
        TDynamicVector tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = (pMem[i] - v.pMem[i]);
        return tmp;
    }
    T operator*(const TDynamicVector& v) const
    {
        if (sz != v.sz) throw std::length_error("Incompatible sizes");
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

    void generate() {
        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_real_distribution<double> coef_gen(0.,1.);

        for (size_t i = 0; i < sz; ++i) {
            pMem[i] = coef_gen(e) * (1 - (2* (rand()%2)));
        }
    }

    template <typename type2>
    operator TDynamicVector<type2>() {
        TDynamicVector<type2> V(size());
        for (size_t i = 0; i < size(); ++i) {
            V[i] = type2(((*this)[i]));
        }
        return V;
    }
};


// Динамическая матрица - 
// шаблонная матрица на динамической памяти
template<typename T>
class TDynamicMatrix : private TDynamicVector<TDynamicVector<T>>
{
    using TDynamicVector<TDynamicVector<T>>::pMem;
    using TDynamicVector<TDynamicVector<T>>::sz;
public:
    TDynamicMatrix(size_t s = 1) : TDynamicVector<TDynamicVector<T>>(s)
    {
        if (sz > MAX_MATRIX_SIZE)
            throw std::out_of_range("Too large matrix");
        for (size_t i = 0; i < sz; i++)
            pMem[i] = TDynamicVector<T>(sz);
    }

    using TDynamicVector<TDynamicVector<T>>::size;
    using TDynamicVector<TDynamicVector<T>>::operator[];
    using TDynamicVector<TDynamicVector<T>>::at;

    // сравнение
    bool operator==(const TDynamicMatrix& m) const noexcept
    {
        return TDynamicVector<TDynamicVector<T>>::operator==(m);
    }

    bool operator!=(const TDynamicMatrix& m) const noexcept // my func
    {
        return TDynamicVector<TDynamicVector<T>>::operator!=(m);
    }

    // матрично-скалярные операции
    TDynamicMatrix<T> operator*(const T& val)
    {
        TDynamicMatrix<T> tmp(sz);
        for (size_t i = 0; i < sz; ++i)
            tmp.pMem[i] = pMem[i] * val;
        return tmp;
    }

    // матрично-векторные операции
    TDynamicVector<T> operator*(const TDynamicVector<T>& v)
    {
        TDynamicVector<T> tmp(sz);
        for (size_t i = 0; i < sz; ++i)
            tmp[i] = pMem[i] * v;
        return tmp;
    }

    // матрично-матричные операции
    TDynamicMatrix operator+(const TDynamicMatrix& m)
    {
        TDynamicMatrix tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] + m.pMem[i];
        return tmp;
    }

    TDynamicMatrix operator-(const TDynamicMatrix& m)
    {
        TDynamicMatrix tmp(sz);
        for (size_t i = 0; i < sz; i++)
            tmp.pMem[i] = pMem[i] - m.pMem[i];
        return tmp;
    }
    TDynamicMatrix operator*(const TDynamicMatrix& m) //square matrixes
    {
        if (sz != m.sz) throw std::length_error("incompatible sizes");
        TDynamicMatrix res(sz);

        for (size_t i = 0; i < sz; ++i) {
            for (size_t k = 0; k < sz; ++k) {
                for (size_t j = 0; j < sz; ++j) {
                    res[i][j] += pMem[i][k] * m[k][j];
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
                sum += std::abs(double(pMem[i][j]));
            if (sum > res) res = sum;
        }
        return res;
    }

    //void LUdecomposition(TDynamicMatrix& L, TDynamicMatrix& U) const {
    //    L = TDynamicMatrix<T>(sz);
    //    U = TDynamicMatrix<T>(sz);
    //    for (size_t i = 0; i < sz; ++i) {
    //        /*
    //        for (size_t j = 0; j < sz; ++j) {
    //            U[i][j] = 0;
    //            L[i][j] = 0;
    //        }
    //        */
    //        L[i][i] = 1;
    //    }

    //    T sum;
    //    for (size_t i = 0; i < sz; ++i) {
    //        for (size_t j = 0; j < sz; ++j) {
    //            sum = 0;
    //            if (i <= j) {
    //                for (size_t k = 0; k <= i; ++k) {
    //                    sum += L[i][k] * U[k][j];
    //                }
    //                U[i][j] = pMem[i][j] - sum;
    //            }
    //            else {
    //                for (size_t k = 0; k <= j; ++k) {
    //                    sum += L[i][k] * U[k][j];
    //                }
    //                L[i][j] = (pMem[i][j] - sum) / U[j][j];
    //            }
    //        }
    //    }
    //}

    // ввод/вывод
    friend istream& operator>>(istream& istr, TDynamicMatrix& m)
    {
        for (size_t i = 0; i < m.sz; ++i) {
            istr >> m.pMem[i];
        }
        return istr;
    }
    friend ostream& operator<<(ostream& ostr, const TDynamicMatrix& m)
    {
        for (size_t i = 0; i < m.sz; ++i)
            ostr << m.pMem[i] << endl;
        return ostr;
    }

    //TDynamicVector<T> solver(TDynamicVector<T>& b) const {
    //    TDynamicMatrix<T> L(sz), U(sz);
    //    LUdecomposition(L, U);
    //    TDynamicVector<T> y(sz);
    //    //        cout << endl << L << endl << U << endl;
    //            // Ly = b, L - upper triangle
    //    for (size_t i = 0; i < sz; ++i) {
    //        y[i] = b[i];
    //        for (size_t j = i; j > 0; --j) {
    //            y[i] -= y[j - 1] * L[i][j - 1];
    //        }
    //    }
    //    //        cout << y << endl << L * y << endl << endl;
    //    TDynamicVector<T> x(sz);
    //    // Ux = y, U - bottom triangle
    //    for (size_t i = sz; i > 0; --i) {
    //        x[i - 1] = y[i - 1];
    //        for (size_t j = i; j < sz; ++j) {
    //            x[i - 1] -= x[j] * U[i - 1][j];
    //        }
    //        x[i - 1] = x[i - 1] / U[i - 1][i - 1];
    //    }
    //    return x;
    //}

    void generate() {
        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_real_distribution<double> coef_gen(1.,2.);
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                pMem[i][j] = T(coef_gen(e) * (1 - (2*(rand()%2))));
            }
//            pMem[i][i]*=2.;
        }
    }

    template <typename type2>
    operator TDynamicMatrix<type2>(){
        TDynamicMatrix<type2> M(size());
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = 0; j < size(); ++j) {
                M[i][j] = type2((*this)[i][j]);
            }
        }
        return M;
    }


    void generateGoodMatrix() {
        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_real_distribution<double> coef_gen(5.0, 100.0);  // Диапазон для диагональных элементов
        std::uniform_real_distribution<double> small_gen(-4., 4.); // Малые значения для недиагональных элементов


        // Заполнение диагонали положительными числами > 1
        for (size_t i = 0; i < sz; ++i) {
            pMem[i][i] = T(coef_gen(e)); // Гарантированно ненулевые элементы (>1)
        }

        // Добавление небольших недиагональных элементов
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                if (i != j) {
                    pMem[i][j] = T(small_gen(e)); // Малые значения для обеспечения доминирования диагонали
                }
            }
        }
        
    }

    
    void generateGoodMatrix2() {
        std::random_device r;
        std::default_random_engine e(r());
        std::uniform_real_distribution<double> coef_gen(0.1,10.);
        std::uniform_real_distribution<double> small_gen(1.001,1.1);
        double sum = 0;
        for (size_t i = 0; i < sz; ++i) {
            sum = 0;
                for (size_t j = 0; j < sz; ++j) {
                    if (i != j) {
                        sum+= abs(pMem[i][j] = T(coef_gen(e) * (1 - (2 * (rand()%2)))));
                    }
                }
                pMem[i][i]=T(small_gen(e) * sum * (1 - (2 * (rand()%2))));
            
            while(pMem[i][i] >10000.){
                for (size_t j = 0; j < sz; ++j) {
                    pMem[i][j]/=2.;
                }
            }
            
        }
    }

};


/*
рекурсивный определитель матрицы разложением
//Возвращает матрицу matrix без row-ой строки и col-того столбца, результат в newMatrix
template <typename type>
void getMatrixWithoutRowAndCol(TDynamicMatrix<type>& matrix, int row, int col, TDynamicMatrix<type>& newMatrix) {
    size_t size = matrix.size();
    int offsetRow = 0; //Смещение индекса строки в матрице
    int offsetCol = 0; //Смещение индекса столбца в матрице
    for (int i = 0; i < size - 1; i++) {
        //Пропустить row-ую строку
        if (i == row) {
            offsetRow = 1; //Как только встретили строку, которую надо пропустить, делаем смещение для исходной матрицы
        }

        offsetCol = 0; //Обнулить смещение столбца
        for (int j = 0; j < size - 1; j++) {
            //Пропустить col-ый столбец
            if (j == col) {
                offsetCol = 1; //Встретили нужный столбец, проускаем его смещением
            }

            newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
        }
    }
}


template <typename type>
type matrixDet(TDynamicMatrix<type> matrix) {
    type det = 0;
    int degree = 1; // (-1)^(1+j) из формулы определителя
    size_t size = matrix.size();

    //Условие выхода из рекурсии
    if (size == 1) {
        return matrix[0][0];
    }
    //Условие выхода из рекурсии
    else if (size == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    else {
        //Матрица без строки и столбца
        TDynamicMatrix<type> newMatrix(size-1);

        //Раскладываем по 0-ой строке, цикл бежит по столбцам
        for (int j = 0; j < size; j++) {
            //Удалить из матрицы i-ю строку и j-ый столбец
            //Результат в newMatrix
            getMatrixWithoutRowAndCol(matrix, 0, j, newMatrix);

            //Рекурсивный вызов
            //По формуле: сумма по j, (-1)^(1+j) * matrix[0][j] * minor_j (это и есть сумма из формулы)
            //где minor_j - дополнительный минор элемента matrix[0][j]
            // (напомню, что минор это определитель матрицы без 0-ой строки и j-го столбца)
            det = det + (degree * matrix[0][j] * matrixDet(newMatrix));
            //"Накручиваем" степень множителя
            degree = -degree;
        }

        //Чистим память на каждом шаге рекурсии(важно!)
        //delete newMatrix;
    }

    return det;
}
*/
template <typename type>
type Det(TDynamicMatrix<type> matrix)
{
    type det = 0;
    int sign = 1;
    size_t size = matrix.size();
    // Base Case
    if (size == 1) {
        det = matrix[0][0];
    }
    else if (size == 2) {
        det = (matrix[0][0] * matrix[1][1])
            - (matrix[0][1] * matrix[1][0]);
    }

    // Perform the Laplace Expansion
    else {
        for (int i = 0; i < size; i++) {

            // Stores the cofactor matrix
            TDynamicMatrix<type> cofactor(size-1);
            int sub_i = 0, sub_j = 0;
            for (int j = 1; j < size; j++) {
                for (int k = 0; k < size; k++) {
                    if (k == i) {
                        continue;
                    }
                    cofactor[sub_i][sub_j] = matrix[j][k];
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }

            // Update the determinant value
            det += sign * matrix[0][i] * Det(cofactor);
                sign = -sign;
        }
    }

    // Return the final determinant value
    return det;
}
template <typename type>
bool IsDiagDiffFromZero(TDynamicMatrix<type> M, type ref) {
    for (size_t i = 0; i < M.size(); ++i) {
        if (M[i][i] <= ref) return false;
    }
    return true;
}

template <typename type>
TDynamicVector<type> DiagMatrix(TDynamicMatrix<type> M) {
    TDynamicVector<type> Diag(M.size());
    for (size_t i = 0; i < M.size(); ++i) Diag[i] = M[i][i];
    return Diag;
}

template <typename type>
bool CloseSol(TDynamicMatrix<type> A, TDynamicVector<type> x, TDynamicVector<type> b, type ref) {
    for (size_t i = 0; i < A.size(); ++i) {
        if ((A * x - b)[i] >= ref){
            double t = (A * x - b)[i];
            t = t;
            x[i]=x[i];
            b[i] = b[i];
            return false;
        }
    }
    return true;
}

template <typename type>
bool IsMatrixDiffFromInf(TDynamicMatrix<type> M, type ref) {
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M.size(); ++j) {
            if (M[i][j] >= ref) return false;
        }
    }
    return true;
}

template <typename type>
type MinVal(TDynamicMatrix<type> M) {
    type temp = M[0][0];
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M.size(); ++j) {
            if (M[i][j] < temp) temp = M[i][j];
        }
    }
    return temp;
}

template <typename type>
type MaxVal(TDynamicMatrix<type> M) {
    type temp = M[0][0];
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M.size(); ++j) {
            if (M[i][j] > temp) temp = M[i][j];
        }
    }
    return temp;
}

template <typename type>
type MinVal(TDynamicVector<type> M) {
    type temp = M[0];
    for (size_t i = 0; i < M.size(); ++i) {
            if (M[i] < temp) temp = M[i];
    }
    return temp;
}

template <typename type>
type MaxVal(TDynamicVector<type> M) {
    type temp = M[0];
    for (size_t i = 0; i < M.size(); ++i) {
            if (M[i] > temp) temp = M[i];
    }
    return temp;
}
#endif
