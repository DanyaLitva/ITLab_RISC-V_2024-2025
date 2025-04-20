#include <iostream>
#include "tmatrix.h"
#include <string>
#include <time.h>
#include "FP/FP16.h"
//#include "FP/FP32.h"

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

			x[i] *= (1 / A[i][i]);
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

			x[i] *= (1 / A[i][i]);
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


	TDynamicMatrix<float> A_fp32(A.size());
	TDynamicVector<float> x_fp32(b.size());
	TDynamicVector<float> b_fp32(b.size());

	

	A_fp32 = TDynamicMatrix<float>(A);
	b_fp32 = TDynamicVector<float>(b);
	x_fp32 = Gauss_Seidel_accurate(A_fp32, x_fp32, b_fp32, float(0.00048828125));


	TDynamicVector<double> x_fp64(b.size());
	x_fp64 = TDynamicVector<float>(x_fp32);

	x_fp64 = Gauss_Seidel_accurate(A, x_fp64, b, ref);

	return x_fp64;
}

int main() {
	time_t t, t1, t2, t3, t4, t5;
	TDynamicMatrix<double> A(100);
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
	VectorTime[0] = double(t1);
	cout << double(t) << ", ";
	cout << (A * x - b) << endl;


	//VectorTimeCount[min_element(VectorTime.begin(), VectorTime.end()) - VectorTime.begin()]++;



	return 0;
}

