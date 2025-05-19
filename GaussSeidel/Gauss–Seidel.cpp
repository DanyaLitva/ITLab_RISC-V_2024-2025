#include <iostream>
#include "tmatrix.h"
#include <string>
#include <chrono>
#include "FP.h"
#include "FP64.h"
#define TYPEMATRIX FP32
#define CloseSol CloseSol2

////for CloseSol
//double fp64perc = 1.e-13;
//float fp32perc = 0.0005;
//float fp16perc = 0.05;

//for CloseSol2
double fp64perc =1.e-24;
float fp32perc = 0.000005;
float fp16perc = 0.5;

size_t count4pre = 1;
size_t count_repeat = 1;
size_t MSize = 2000;
size_t count_it = 1;
//true = triangle
bool Generate = true;
int Mode = -1;
bool WriteMatrix = false;

template<typename type>
TDynamicVector<type> Gauss_Seidel_accurate(TDynamicMatrix<type> A, TDynamicVector<type> x, TDynamicVector<type> b, type ref = fp64perc) {
    TDynamicVector<type> temp(b.size());
    
    //initial approximation
    temp = x;
    auto start = std::chrono::steady_clock::now();
    size_t count_it = 0;
    while(true) {
        //    while(!CloseSol(A, x, b, ref)) {
        start = std::chrono::steady_clock::now();
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
                
                x[i] /= A[i][i];
            }
            if (temp == x) {            
                    cout << "CYCLE after " << count_it << " iterations" << endl;
                    if (Mode == 0) {
                        auto end = std::chrono::steady_clock::now();
                        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                        cout << "last it time: " << elapsed / 1000. << " seconds" << endl;
                    }
                return x;
            }

            temp = x;
            count_it++;
        }
    }
    if (Mode == 0) {
        cout << "count it: " << count_it << endl;
        auto end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << "last it time: " << elapsed / 1000. << " seconds" << endl;
    }
    return x;
}



TDynamicVector<double> Gauss_Seidel_double(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicVector<double> x_fp64(b.size());
    for (size_t i = 0; i < x_fp64.size(); ++i) x_fp64[i] = 1.;

    if(Mode == 0) cout<<"double ";
    auto start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A, x_fp64, b, ref);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    if(Mode == 0) cout << "double time: "<<elapsed/1000. <<" seconds"<< endl;
    
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
    if (Mode == 0)cout<<endl<<"float ";
    auto start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, fp32perc);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    if (Mode == 0) cout << "float time: "<<elapsed/1000. <<" seconds"<< endl;


    TDynamicVector<double> x_fp64(b.size());

//        start = std::chrono::steady_clock::now();
    x_fp64 = TDynamicVector<double>(x_fp32);
//        end = std::chrono::steady_clock::now();
//        elapsed = end - start;
//        cout<<"float to double:"<<elapsed<<endl;
    
    if (Mode==0)cout<<"double ";
    start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A, x_fp64, b, ref);
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    if (Mode==0)cout << "double time: "<<elapsed/1000. <<" seconds"<< endl;
    
    return x_fp64;
}

TDynamicVector<double> Gauss_Seidel_fp64(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP64> A_fp64(A.size());
    TDynamicVector<FP64> x_fp64(b.size());
    TDynamicVector<FP64> b_fp64(b.size());

    for (size_t i = 0; i < x_fp64.size(); ++i) x_fp64[i] = 1.;

    A_fp64 = TDynamicMatrix<FP64>(A);
    b_fp64 = TDynamicVector<FP64>(b);
    if (Mode == 0) cout << "FP64 ";
    auto start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A_fp64, x_fp64, b_fp64, FP64(ref));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    if (Mode == 0)  cout << "FP64 time: " << elapsed / 1000. << " seconds" << endl;
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
    if (Mode==0)cout<<"FP16 ";
    auto start = std::chrono::steady_clock::now();
    x_fp16 = Gauss_Seidel_accurate(A_fp16, x_fp16, b_fp16, FP16(fp16perc));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    
    
    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());


    bool flag = false;
    for (size_t i = 0; i < x_fp16.size(); ++i) {
        if (x_fp16[i].IsNan() || x_fp16[i].IsInf()) {
            flag = true;
            break;
        }
    }
    if (flag) {
        for (size_t i = 0; i < x_fp16.size(); ++i) x_fp32[i] = 1.;
        cout << "FP16 Inf, ";
    }
    else x_fp32 = TDynamicVector<FP32>(x_fp16);

    if (Mode==0)    cout << "FP16 time: " << elapsed / 1000. << " seconds" << endl;

    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    if (Mode==0)    cout<<"FP32 ";
    start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(fp32perc));
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    if (Mode==0)   cout << "FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    
    return x_fp32;
}

TDynamicVector<double> Gauss_Seidel_fp32_fp64(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());

    //initial approximation
    FP32 temp(1.f);
    for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = temp;

    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    if (Mode == 0)cout << "FP32 ";
    auto start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(fp32perc));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


    TDynamicMatrix<FP64> A_fp64(A.size());
    TDynamicVector<FP64> x_fp64(b.size());
    TDynamicVector<FP64> b_fp64(b.size());



    x_fp64 = TDynamicVector<FP64>(x_fp32);

    if (Mode == 0)    cout << "FP32 time: " << elapsed / 1000. << " seconds" << endl;

    A_fp64 = TDynamicMatrix<FP64>(A);
    b_fp64 = TDynamicVector<FP64>(b);
    if (Mode == 0)    cout << "FP64 ";
    start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A_fp64, x_fp64, b_fp64, FP64(fp64perc));
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
    if (Mode == 0)   cout << "FP64 time: " << elapsed / 1000. << " seconds" << endl;

    return x_fp64;
}

TDynamicVector<double> Gauss_Seidel_fp16_fp32_fp64(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP16> A_fp16(A.size());
    TDynamicVector<FP16> x_fp16(b.size());
    TDynamicVector<FP16> b_fp16(b.size());

    //initial approximation
    FP16 temp(0, 15, 0);
    for (size_t i = 0; i < x_fp16.size(); ++i) x_fp16[i] = temp;

    A_fp16 = TDynamicMatrix<FP16>(A);
    b_fp16 = TDynamicVector<FP16>(b);
    if (Mode == 0)cout << "FP16 ";
    auto start = std::chrono::steady_clock::now();
    x_fp16 = Gauss_Seidel_accurate(A_fp16, x_fp16, b_fp16, FP16(fp16perc));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());

    ////initial approximation
    //FP32 temp(1.f);
    //for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = temp;

    x_fp32 = TDynamicVector<FP32>(x_fp16);

    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    if (Mode == 0)cout << "FP32 ";
    start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(fp32perc));
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


    TDynamicMatrix<FP64> A_fp64(A.size());
    TDynamicVector<FP64> x_fp64(b.size());
    TDynamicVector<FP64> b_fp64(b.size());



    x_fp64 = TDynamicVector<FP64>(x_fp32);

    if (Mode == 0)    cout << "FP32 time: " << elapsed / 1000. << " seconds" << endl;

    A_fp64 = TDynamicMatrix<FP64>(A);
    b_fp64 = TDynamicVector<FP64>(b);
    if (Mode == 0)    cout << "FP64 ";
    start = std::chrono::steady_clock::now();
    x_fp64 = Gauss_Seidel_accurate(A_fp64, x_fp64, b_fp64, FP64(fp64perc));
    end = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
    if (Mode == 0)   cout << "FP64 time: " << elapsed / 1000. << " seconds" << endl;

    return x_fp64;
}

TDynamicVector<double> Gauss_Seidel_fp32(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());

    for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = 1.;
    
    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    if (Mode==0) cout<<"FP32 ";
    auto start = std::chrono::steady_clock::now();
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(ref));
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds > (end - start).count();
    if (Mode==0)  cout << "FP32 time: "<<elapsed/1000. <<" seconds"<< endl;
    return x_fp32;
}




using namespace std;

int main(int argc, char* argv[]) {
    if (Mode == 0) {
        TDynamicMatrix<double> A(MSize);
        TDynamicVector<double> b(A.size());
        TDynamicVector<double> x(A.size());
        TDynamicMatrix<FP32> A_temp(MSize);
        TDynamicVector<FP32> b_temp(A.size());
        TDynamicVector<FP32> x_temp(A.size());
        cout << "Matrix size: " << MSize << "*" << MSize << endl;
        //cout<<"number of preiterations: "<<count4pre<<endl;
        cout << "number of iterations before checking the accuracy: " << count_repeat << endl;

        TDynamicMatrix<double> minA(MSize);
        TDynamicMatrix<double> maxA(MSize);
        TDynamicVector<double> minb(A.size());
        TDynamicVector<double> maxb(A.size());
        double minT = 10000000.;
        double maxT = 0;

        for (size_t temp_it = 0; temp_it < count_it; temp_it++) {
            if(Generate)     A.generateGoodMatrix();
            else A.generateGoodMatrix2();
            x.generate();
            A_temp = TDynamicMatrix<TYPEMATRIX>(A);
            x_temp = TDynamicVector<TYPEMATRIX>(x);
            b_temp = A_temp * x_temp;
            //b.generate();
            A = TDynamicMatrix<double>(A_temp);
            b = TDynamicVector<double>(b_temp);

            cout << endl << endl << "Iteration #" << temp_it + 1 << endl;
            cout << "About Matrix:" << endl;
            cout << "A min: " << MinVal(A) << endl;
            cout << "A max: " << MaxVal(A) << endl;
            cout << "x min: " << MinVal(x) << endl;
            cout << "x max: " << MaxVal(x) << endl;
            cout << "b min: " << MinVal(b) << endl;
            cout << "b max: " << MaxVal(b) << endl << endl;
            if (WriteMatrix == true) {
                cout << A << endl << endl;
                cout << "x: " << x << endl << endl;
                cout << "b: " << b << endl << endl;
            }
            auto start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_double(A, b, fp64perc);
            auto end = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << "general double time: " << elapsed / 1000. << " seconds" << endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl;
            //cout << A*x-b << endl;


            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_float_double(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << "general float+double time: " << elapsed / 1000. << " seconds" << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;


            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp32(A, b, fp32perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << "general FP32 time: " << elapsed / 1000. << " seconds" << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;


            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp16_fp32(A, b, fp32perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            cout << "general FP16+FP32 time: " << elapsed / 1000. << " seconds" << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp64(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << "general FP64 time: " << elapsed / 1000. << " seconds" << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;
            /*cout << A * x - b << endl;*/


            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp32_fp64(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            cout << "general FP32+FP64 time: " << elapsed / 1000. << " seconds" << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp16_fp32_fp64(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            cout << "general FP16+FP32+FP64 time: " << elapsed / 1000. << " seconds" << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

        }
        return 0;
    }
    if (Mode==1) {
        TDynamicMatrix<double> A(MSize);
        TDynamicVector<double> b(A.size());
        TDynamicVector<double> x(A.size());
        TDynamicMatrix<FP32> A_temp(MSize);
        TDynamicVector<FP32> b_temp(A.size());
        TDynamicVector<FP32> x_temp(A.size());
        cout << MSize << endl;
        //cout<<"number of preiterations: "<<count4pre<<endl;
        //cout << "number of iterations before checking the accuracy: " << count_repeat << endl;

        TDynamicMatrix<double> minA(MSize);
        TDynamicMatrix<double> maxA(MSize);
        TDynamicVector<double> minb(A.size());
        TDynamicVector<double> maxb(A.size());
        double minT = 10000000.;
        double maxT = 0;

        if (Generate)     A.generateGoodMatrix();
        else A.generateGoodMatrix2();
        x.generate();
        A_temp = TDynamicMatrix<FP32>(A);
        x_temp = TDynamicVector<FP32>(x);
        b_temp = A_temp * x_temp;
        //b.generate();
        A = TDynamicMatrix<double>(A_temp);
        b = TDynamicVector<double>(b_temp);

        //cout << endl << endl << "Iteration #" << temp_it + 1 << endl;
        //cout << "About Matrix:" << endl;
        //cout << "A min: " << MinVal(A) << endl;
        //cout << "A max: " << MaxVal(A) << endl;
        //cout << "x min: " << MinVal(x) << endl;
        //cout << "x max: " << MaxVal(x) << endl;
        //cout << "b min: " << MinVal(b) << endl;
        //cout << "b max: " << MaxVal(b) << endl << endl;
        if (WriteMatrix == true) {
            cout << A << endl << endl;
            cout << "x: " << x << endl << endl;
            cout << "b: " << b << endl << endl;
        }
        auto start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_double(A, b, fp64perc);
        auto end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << elapsed / 1000. << endl;
        //cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
        //cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl;

        //cout << x << endl;

        start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_float_double(A, b, fp64perc);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << elapsed / 1000. << endl;
        //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        //    cout<<"x max: "<<MaxVal(x)<<endl;
        //cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
        //cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

        //start = std::chrono::steady_clock::now();
        //x = Gauss_Seidel_fp32(A, b, fp32perc);
        //end = std::chrono::steady_clock::now();
        //elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        //cout <<  elapsed / 1000. << endl;
        ////    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        ////    cout<<"x max: "<<MaxVal(x)<<endl;
        ///*cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
        //cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;*/



        //start = std::chrono::steady_clock::now();
        //x = Gauss_Seidel_fp16_fp32(A, b, fp32perc);
        //end = std::chrono::steady_clock::now();
        //elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        //cout << elapsed / 1000.  << endl;
        ////    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        ////    cout<<"x max: "<<MaxVal(x)<<endl;
        ////cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
        ////cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

        start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_fp64(A, b, fp64perc);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << elapsed / 1000. << endl;
        //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        //    cout<<"x max: "<<MaxVal(x)<<endl;
        /*cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
        cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;*/



        //start = std::chrono::steady_clock::now();
        //x = Gauss_Seidel_fp32_fp64(A, b, fp64perc);
        //end = std::chrono::steady_clock::now();
        //elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        //cout << elapsed / 1000. << endl;
        ////    cout<<endl<<"x min: "<<MinVal(x)<<endl;
        ////    cout<<"x max: "<<MaxVal(x)<<endl;
        ////cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
        ////cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;


        start = std::chrono::steady_clock::now();
        x = Gauss_Seidel_fp32_fp64(A, b, fp64perc);
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        cout << elapsed / 1000. << endl;
        return 0;
    }
    if (Mode == -1) {
        TDynamicMatrix<double> A_double(MSize);
        TDynamicVector<double> b_double(MSize);
        TDynamicVector<double> x_double(MSize);
        TDynamicMatrix<float> A_float(MSize);
        TDynamicVector<float> b_float(MSize);
        TDynamicVector<float> x_float(MSize);
        TDynamicMatrix<FP32> A_FP32(MSize);
        TDynamicVector<FP32> b_FP32(MSize);
        TDynamicVector<FP32> x_FP32(MSize);
        TDynamicMatrix<FP16> A_FP16(MSize);
        TDynamicVector<FP16> b_FP16(MSize);
        TDynamicVector<FP16> x_FP16(MSize);
        TDynamicMatrix<FP64> A_FP64(MSize);
        TDynamicVector<FP64> b_FP64(MSize);
        TDynamicVector<FP64> x_FP64(MSize);
        A_double.generateGoodMatrix();
        x_double.generate();
        b_double = A_double * x_double;

        A_float.generateGoodMatrix();
        x_float.generate();
        b_float = A_float * x_float;

        A_FP16.generateGoodMatrix();
        x_FP16.generate();
        b_FP16 = A_FP16 * x_FP16;

        A_FP32.generateGoodMatrix();
        x_FP32.generate();
        b_FP32 = A_FP32 * x_FP32;

        A_FP64.generateGoodMatrix();
        x_FP64.generate();
        b_FP64 = A_FP64 * x_FP64;
        auto start = std::chrono::steady_clock::now();
        x_float = A_float * x_float - b_float;
        auto end = std::chrono::steady_clock::now();
        float elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << "float time: " << elapsed / 1000 << endl;

        start = std::chrono::steady_clock::now();
        A_double* x_double - b_double;
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << "double time: " << elapsed / 1000 << endl;

        //start = std::chrono::steady_clock::now();
        //A_FP16* x_FP16 - b_FP16;
        //end = std::chrono::steady_clock::now();
        //elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        //cout << "FP16 time: " << elapsed / 1000 << endl;

        start = std::chrono::steady_clock::now();
        A_FP32* x_FP32 - b_FP32;
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << "FP32 time: " << elapsed/1000 << endl;

        start = std::chrono::steady_clock::now();
        A_FP64* x_FP64 - b_FP64;
        end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
        cout << "FP64 time: " << elapsed / 1000 << endl;



        return 0;
    }

    if (Mode == 2) {
        for (size_t MSize = 100; MSize < 2100; MSize += 100) {
            TDynamicMatrix<double> A(MSize);
            TDynamicVector<double> b(A.size());
            TDynamicVector<double> x(A.size());
            TDynamicMatrix<FP32> A_temp(MSize);
            TDynamicVector<FP32> b_temp(A.size());
            TDynamicVector<FP32> x_temp(A.size());
            cout << MSize << endl;
            //cout<<"number of preiterations: "<<count4pre<<endl;
            //cout << "number of iterations before checking the accuracy: " << count_repeat << endl;

            TDynamicMatrix<double> minA(MSize);
            TDynamicMatrix<double> maxA(MSize);
            TDynamicVector<double> minb(A.size());
            TDynamicVector<double> maxb(A.size());
            double minT = 10000000.;
            double maxT = 0;

            if (Generate)     A.generateGoodMatrix();
            else A.generateGoodMatrix2();
            x.generate();
            A_temp = TDynamicMatrix<FP32>(A);
            x_temp = TDynamicVector<FP32>(x);
            b_temp = A_temp * x_temp;
            //b.generate();
            A = TDynamicMatrix<double>(A_temp);
            b = TDynamicVector<double>(b_temp);

            //cout << endl << endl << "Iteration #" << temp_it + 1 << endl;
            //cout << "About Matrix:" << endl;
            //cout << "A min: " << MinVal(A) << endl;
            //cout << "A max: " << MaxVal(A) << endl;
            //cout << "x min: " << MinVal(x) << endl;
            //cout << "x max: " << MaxVal(x) << endl;
            //cout << "b min: " << MinVal(b) << endl;
            //cout << "b max: " << MaxVal(b) << endl << endl;
            if (WriteMatrix == true) {
                cout << A << endl << endl;
                cout << "x: " << x << endl << endl;
                cout << "b: " << b << endl << endl;
            }
            auto start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_double(A, b, fp64perc);
            auto end = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << elapsed / 1000. << endl;
            //cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            //cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl;

            //cout << x << endl;

            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_float_double(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << elapsed / 1000. << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            //cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            //cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

            //start = std::chrono::steady_clock::now();
            //x = Gauss_Seidel_fp32(A, b, fp32perc);
            //end = std::chrono::steady_clock::now();
            //elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            //cout <<  elapsed / 1000. << endl;
            ////    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            ////    cout<<"x max: "<<MaxVal(x)<<endl;
            ///*cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            //cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;*/



            //start = std::chrono::steady_clock::now();
            //x = Gauss_Seidel_fp16_fp32(A, b, fp32perc);
            //end = std::chrono::steady_clock::now();
            //elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            //cout << elapsed / 1000.  << endl;
            ////    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            ////    cout<<"x max: "<<MaxVal(x)<<endl;
            ////cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            ////cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;

            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp64(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end - start).count();
            cout << elapsed / 1000. << endl;
            //    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            //    cout<<"x max: "<<MaxVal(x)<<endl;
            /*cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;*/



            //start = std::chrono::steady_clock::now();
            //x = Gauss_Seidel_fp32_fp64(A, b, fp64perc);
            //end = std::chrono::steady_clock::now();
            //elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            //cout << elapsed / 1000. << endl;
            ////    cout<<endl<<"x min: "<<MinVal(x)<<endl;
            ////    cout<<"x max: "<<MaxVal(x)<<endl;
            ////cout << "(Ax-b) min: " << MinVal((A * x - b)) << endl;
            ////cout << "(Ax-b) max: " << MaxVal((A * x - b)) << endl << endl;


            start = std::chrono::steady_clock::now();
            x = Gauss_Seidel_fp32_fp64(A, b, fp64perc);
            end = std::chrono::steady_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            cout << elapsed / 1000. << endl;
        }
            return 0;
    }
    return 0;
}

//выдать количество итераций для precondition

//спектральное разложение сами задаем хар числа
//отношение макс к минимуму хар чисел
//отношение макс к минимуму хар чисел
//
//задать генератор характер чисел от 0.1 до 5 например
