#include <iostream>
#include "tmatrix.h"
#include <string>
#include <chrono>
#include "FP.h"

//auto start = std::chrono::steady_clock::now();
//
//auto end = std::chrono::steady_clock::now();
//    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
//cout << "double time: "<<elapsed/1000. <<" seconds"<< endl;

double fp64perc =1.e-13;
float fp32perc =0.0001;
float fp16perc =0.5;
size_t count4pre = 10;
size_t count_repeat = 1;
size_t MSize = 512;
size_t count_it = 10;


template<typename type>
TDynamicVector<type> Gauss_Seidel_accurate(TDynamicMatrix<type> A, TDynamicVector<type> x, TDynamicVector<type> b, type ref = fp64perc) {
    TDynamicVector<type> temp(b.size());

    //initial approximation
    temp = x;

    size_t count_it = 0;
    while(true) {
//    while(!CloseSol(A, x, b, ref)) {
                auto start = std::chrono::steady_clock::now();
        if(CloseSol(A, x, b, ref)) break;
        for(size_t repeat = 0; repeat<count_repeat; repeat++){
            for (size_t i = 0; i < b.size(); ++i) {
                x[i] = b[i];
                
                for (size_t j = 0; j < i; ++j) {
                    x[i] -= A[i][j] * x[j];
                }
                
                for (size_t j = i + 1; j < b.size(); ++j) {
                    x[i] -= A[i][j] * temp[j];
                }
                
                x[i] *= (type(1.) / A[i][i]);
            }
            temp = x;
            count_it++;

//            CloseSol(A, x, b, ref);
            auto end = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
            cout << "it time: "<<elapsed/1000. <<" seconds"<< endl;
        }
       

    }
    cout<<"count it: "<<count_it<<endl;
    return x;
}


template<typename type>
TDynamicVector<type> Gauss_Seidel_4pre(TDynamicMatrix<type> A, TDynamicVector<type> x, TDynamicVector<type> b, type ref = fp64perc) {
    TDynamicVector<type> temp(b.size());

    //initial approximation
    temp = x;

    size_t count_it = 0;
    for(size_t i = 0; i<count4pre;++i){
    //while(!CloseSol(A, x, b, ref)) {
//        auto start = std::chrono::steady_clock::now();
        for (size_t i = 0; i < b.size(); ++i) {
            x[i] = b[i];

            for (size_t j = 0; j < i; ++j) {
                x[i] -= A[i][j] * x[j];
            }

            for (size_t j = i + 1; j < b.size(); ++j) {
                x[i] -= A[i][j] * temp[j];
            }

            x[i] *= (type(1.) / A[i][i]);
        }
        temp = x;
        count_it++;
//        auto start = std::chrono::steady_clock::now();
//        CloseSol(A, x, b, ref);
//        auto end = std::chrono::steady_clock::now();
//            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
//        cout << "it time: "<<elapsed/1000. <<" seconds"<< endl;
    }
    cout<<"count it: "<<count_it<<endl;
    return x;
}


TDynamicVector<double> Gauss_Seidel_double(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicVector<double> x_fp64(b.size());
    for (size_t i = 0; i < x_fp64.size(); ++i) x_fp64[i] = 1.;

    cout<<"double ";
    auto start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A, x_fp64, b, ref);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "double time: "<<elapsed/1000. <<" seconds"<< endl;
    
    return x_fp64;
}


TDynamicVector<double> Gauss_Seidel_float_double(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<float> A_fp32(A.size());
    TDynamicVector<float> x_fp32(b.size());
    TDynamicVector<float> b_fp32(b.size());
    //initial approximation
    for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = 1.f;

    A_fp32 = TDynamicMatrix<float>(A);
    b_fp32 = TDynamicVector<float>(b);
    cout<<endl<<"float ";
    auto start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, fp32perc);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "float time: "<<elapsed/1000. <<" seconds"<< endl;


    TDynamicVector<double> x_fp64(b.size());

//        start = std::chrono::steady_clock::now();
    x_fp64 = TDynamicVector<double>(x_fp32);
//        end = std::chrono::steady_clock::now();
//        elapsed = end - start;
//        cout<<"float to double:"<<elapsed<<endl;
    
    cout<<"double ";
    start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A, x_fp64, b, ref);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "double time: "<<elapsed/1000. <<" seconds"<< endl;
    
    return x_fp64;
}

TDynamicVector<double> Gauss_Seidel_fp16_fp32(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP16> A_fp16(A.size());
    TDynamicVector<FP16> x_fp16(b.size());
    TDynamicVector<FP16> b_fp16(b.size());

    //initial approximation
    FP16 temp(0, 15, 0);
    for (size_t i = 0; i < x_fp16.size(); ++i) x_fp16[i] = temp;

    A_fp16 = TDynamicMatrix<FP16>(A);
    b_fp16 = TDynamicVector<FP16>(b);
    cout<<"FP16 ";
    auto start = std::chrono::steady_clock::now();
    x_fp16 = Gauss_Seidel_accurate(A_fp16, x_fp16, b_fp16, FP16(fp16perc));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "FP16 time: "<<elapsed/1000. <<" seconds"<< endl;
    
    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());

//    start = std::chrono::steady_clock::now();
    x_fp32 = TDynamicVector<FP32>(x_fp16);
//    end = std::chrono::steady_clock::now();
//    elapsed = end - start;
//    cout<<"FP16 to FP32:"<<elapsed<<endl;
    
    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    cout<<"FP32 ";
    start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(fp32perc));
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    
    return x_fp32;
}

TDynamicVector<double> Gauss_Seidel_fp32(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());

    for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = 1.;
    
    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    cout<<"FP32 ";
    auto start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(ref));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    return x_fp32;
}




using namespace std;

int main() {
    TDynamicMatrix<double> A(MSize);
    TDynamicVector<double> b(A.size());
    TDynamicVector<double> x(A.size());
    vector<double> VectorTime = { 0,0,0,0,0 };
    vector<int> VectorTimeCount = { 0,0,0,0,0 };
    cout<<"Matrix size: "<<MSize<<"*"<<MSize<<endl;
    cout<<"number of preiterations: "<<count4pre<<endl;
    cout<<"number of iterations before checking the accuracy: "<<count_repeat<<endl;
    
    TDynamicMatrix<double> minA(MSize);
    TDynamicMatrix<double> maxA(MSize);
    TDynamicVector<double> minb(A.size());
    TDynamicVector<double> maxb(A.size());
    double minT = 10000000.;
    double maxT = 0;
    
    for(size_t temp_it = 0; temp_it < count_it;temp_it++){
        A.generateGoodMatrix2();
        b.generate();
        cout<<endl<<endl<<"Iteration #"<<temp_it+1<<endl;
        cout<<"About Matrix:"<<endl;
        cout<<"A min: "<<MinVal(A)<<endl;
        cout<<"A max: "<<MaxVal(A)<<endl;
        cout<<"b min: "<<MinVal(b)<<endl;
        cout<<"b max: "<<MaxVal(b)<<endl<<endl;
        
        auto start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_double(A, b, fp64perc);
        auto end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
        cout << "general double time: "<<elapsed/1000. <<" seconds"<< endl;
        cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
        cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
        

        
        start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_float_double(A, b, fp64perc);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
        cout << "general float+double time: "<<elapsed/1000. <<" seconds"<< endl;
        //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        //    cout<<"x max: "<<MaxVal(x)<<endl;
        cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
        cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
        
        start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_fp32(A, b, fp32perc);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
        cout << "general FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
        //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        //    cout<<"x max: "<<MaxVal(x)<<endl;
        cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
        cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
        
        if(elapsed > maxT) {
            maxT = elapsed;
                maxA = A;
                maxb = b;
            }
            if(elapsed < minT) {
            minT = elapsed;
                minA = A;
                minb = b;
            }
        
        start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_fp16_fp32(A, b, fp32perc);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
        cout << "general FP16+FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
        //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        //    cout<<"x max: "<<MaxVal(x)<<endl;
        cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
        cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    }
        
    
    
    
    
    
    A = minA;
    b = minb;
    
    cout<<endl<<endl<<"Best Iteration"<<endl;
    cout<<"About Matrix:"<<endl;
    cout<<"A min: "<<MinVal(A)<<endl;
    cout<<"A max: "<<MaxVal(A)<<endl;
    cout<<"b min: "<<MinVal(b)<<endl;
    cout<<"b max: "<<MaxVal(b)<<endl<<endl;
    
    auto start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_double(A, b, fp64perc);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general double time: "<<elapsed/1000. <<" seconds"<< endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_float_double(A, b, fp64perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general float+double time: "<<elapsed/1000. <<" seconds"<< endl;
    //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
    //    cout<<"x max: "<<MaxVal(x)<<endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_fp32(A, b, fp32perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
    //    cout<<"x max: "<<MaxVal(x)<<endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_fp16_fp32(A, b, fp32perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general FP16+FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
    //    cout<<"x max: "<<MaxVal(x)<<endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    
    
    A = maxA;
    b = maxb;
    cout<<endl<<endl<<"Worst Iteration"<<endl;
    cout<<"About Matrix:"<<endl;
    cout<<"A min: "<<MinVal(A)<<endl;
    cout<<"A max: "<<MaxVal(A)<<endl;
    cout<<"b min: "<<MinVal(b)<<endl;
    cout<<"b max: "<<MaxVal(b)<<endl<<endl;
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_double(A, b, fp64perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general double time: "<<elapsed/1000. <<" seconds"<< endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
 
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_float_double(A, b, fp64perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general float+double time: "<<elapsed/1000. <<" seconds"<< endl;
    //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
    //    cout<<"x max: "<<MaxVal(x)<<endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_fp32(A, b, fp32perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
    //    cout<<"x max: "<<MaxVal(x)<<endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    
    start = std::chrono::steady_clock::now();
    x = Gauss_Seidel_fp16_fp32(A, b, fp32perc);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    cout << "general FP16+FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
    //    cout<<"x max: "<<MaxVal(x)<<endl;
    cout<<"(Ax-b) min: "<<MinVal((A*x-b))<<endl;
    cout<<"(Ax-b) max: "<<MaxVal((A*x-b))<<endl<<endl;
    
    
    return 0;
}

//выдать количество итераций для precondition

//спектральное разложение сами задаем хар числа
//отношение макс к минимуму хар чисел
//отношение макс к минимуму хар чисел
//
//задать генератор характер чисел от 0.1 до 5 например
