#include <iostream>
#include "tmatrix.h"
#include <string>
#include <time.h>
#include "FP/FP.h"




template<typename type>
TDynamicVector<type> Gauss_Seidel_accurate(TDynamicMatrix<type> A, TDynamicVector<type> x, TDynamicVector<type> b, type ref = 1.e-8) {
	TDynamicVector<type> temp(b.size());

	//initial approximation
	temp = x;

	size_t count_it = 0;
	while(!CloseSol(A, x, b, ref)) {
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
	}
	return x;
}

template<typename type>
TDynamicVector<type> Gauss_Seidel_close(TDynamicMatrix<type> A, TDynamicVector<type> x, TDynamicVector<type> b, type ref = 1.e-8) {
	TDynamicVector<type> temp(b.size());

	//initial approximation
	temp = x;
	bool flag;
	size_t count_it = 0;
	while (!CloseSol(A, x, b, ref)) {
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

		flag = true;
		for (size_t i = 0; i < x.size(); ++i) if (abs(x[i] - temp[i]) >= ref) flag = false;
		if (flag == true) break;

		temp = x;
		count_it++;
	}
	return x;
}




TDynamicVector<double> Gauss_Seidel(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
	TDynamicMatrix<float> A_fp32(A.size());
	TDynamicVector<float> x_fp32(b.size());
	TDynamicVector<float> b_fp32(b.size());

	//initial approximation
	for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = 1.f;

	A_fp32 = TDynamicMatrix<float>(A);
	b_fp32 = TDynamicVector<float>(b);
	x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, float(0.00048828125));


	TDynamicVector<double> x_fp64(b.size());
	x_fp64 = TDynamicVector<float>(x_fp32);

	x_fp64 = Gauss_Seidel_accurate(A, x_fp64, b, ref);

	return x_fp64;
}

TDynamicVector<double> Gauss_Seidel2(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
	TDynamicMatrix<FP16> A_fp16(A.size());
	TDynamicVector<FP16> x_fp16(b.size());
	TDynamicVector<FP16> b_fp16(b.size());

	//initial approximation
	FP16 temp(0, 15, 0);
	for (size_t i = 0; i < x_fp16.size(); ++i) x_fp16[i] = temp;

    A_fp16 = TDynamicMatrix<FP16>(A);
    b_fp16 = TDynamicVector<FP16>(b);
    x_fp16 = Gauss_Seidel_accurate(A_fp16, x_fp16, b_fp16, FP16(0,10,0));

	TDynamicMatrix<FP32> A_fp32(A.size());
	TDynamicVector<FP32> x_fp32(b.size());
	TDynamicVector<FP32> b_fp32(b.size());

    x_fp32 = TDynamicVector<FP16>(x_fp16);
    
	A_fp32 = TDynamicMatrix<FP32>(A);
	b_fp32 = TDynamicVector<FP32>(b);
	x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(float(0.00048828125)));

    
	return x_fp32;
}

TDynamicVector<double> Gauss_Seidel3(TDynamicMatrix<double> A, TDynamicVector <double> b, double ref) {
    TDynamicMatrix<FP32> A_fp32(A.size());
    TDynamicVector<FP32> x_fp32(b.size());
    TDynamicVector<FP32> b_fp32(b.size());

    for (size_t i = 0; i < x_fp32.size(); ++i) x_fp32[i] = 1.;
    
    A_fp32 = TDynamicMatrix<FP32>(A);
    b_fp32 = TDynamicVector<FP32>(b);
    x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, FP32(float(0.00048828125)));
    
    return x_fp32;
}

using namespace std;

int main() {
	time_t t, t1, t2, t3, t4, t5;
	TDynamicMatrix<double> A(10);
	TDynamicVector<double> b(A.size());
	TDynamicVector<double> x(A.size());
	vector<double> VectorTime = { 0,0,0,0,0 };
	vector<int> VectorTimeCount = { 0,0,0,0,0 };
	A.generateGoodMatrix();
	b.generate();

    t = clock();
    x = Gauss_Seidel(A, b, 1.e-10);
    t -= clock();
    t1 = t *= -1;
    cout << double(t) << ", ";
    
	t = clock();
	x = Gauss_Seidel2(A, b, 1.e-10);
	t -= clock();
	t1 = t *= -1;
	cout << double(t) << ", ";
	
    t = clock();
    x = Gauss_Seidel3(A, b, 1.e-10);
    t -= clock();
    t1 = t *= -1;
    cout << double(t) << ", ";

    cout << (A * x - b) << endl;
	VectorTimeCount[min_element(VectorTime.begin(), VectorTime.end()) - VectorTime.begin()]++;

    cout<<endl;
    t = clock();
    FP16 fp16_1, fp16_2;
    fp16_1 = 1.f;
    fp16_2=1.f;
    for (size_t i = 0; i<100000; ++i){
        fp16_1 = fp16_1 + fp16_2;
    }
    t -= clock();
    t1 = t *= -1;
    cout << double(t) << ", ";
    
    t = clock();
    FP32 fp32_1, fp32_2;
    fp32_1 = 1.f;
    fp32_2=1.f;
    for (size_t i = 0; i<100000; ++i){
        fp32_1 = fp32_1 + fp32_2;
    }
    t -= clock();
    t1 = t *= -1;
    cout << double(t) << ", ";
        
    t = clock();
    BF16 bf16_1, bf16_2;
    bf16_1 = 1.f;
    bf16_2=1.f;
    for (size_t i = 0; i<100000; ++i){
        bf16_1 = bf16_1 + bf16_2;
    }
    t -= clock();
    t1 = t *= -1;
    cout << double(t) << endl;
    
    TDynamicMatrix<BF16> B(100);
    for(size_t i = 0; i<0; ++i) B[i][i]=1.f;
    t = clock();
    B*B;
    t -= clock();
    t1 = t *= -1;
    cout << double(t) << endl;

    FP16 Q;
    Q = 10.f;
    Q+=10.;
    cout<<Q<<endl;
    float q = Q;
    cout<<q<<endl;
    Q = q;
    cout<<Q<<endl;
    q = Q;
    cout<<q;
    
    
    
	return 0;
}

