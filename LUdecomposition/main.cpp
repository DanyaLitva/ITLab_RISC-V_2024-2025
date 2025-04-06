#include <iostream>
#include "tmatrixMinimal.h"
#include "../FP/BF16.h"
#include "../FP/FP32.h"
#include <chrono>
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS

typedef double type;

double testLU(TDynamicMatrix<type>& A, bool coutflag) {
	if (coutflag) std::cout << A << std::endl;
	TDynamicMatrix<type> L(A.size()), U(A.size());
	TDynamicMatrix<type> tmp = A;
//	cout << "LU started" << endl;
//	A.LUdecompositionOptimized(L, U);
	A.LUdecompositionV3(L, U, 8);
//	cout << "LU ended" << endl;
	A = tmp;
	tmp = L * U;
	if (coutflag) std::cout << L << std::endl << U << std::endl << tmp << std::endl << A - tmp << std::endl << (A - tmp).norm1() << std::endl << endl;
	return (A - tmp).norm1();
}

void test() {
	TDynamicMatrix<type> A(1003);
	TDynamicMatrix<type> savedA;
	double maxerr = 0.0;
	double err = 0.0;
	vector<double> input = { 0.113249, -4.77064, 0.230813, -0.0496871, 1.0, -0.00777988, 0.0, 0.0, 1.0 };
	for (size_t i = 0; i < 1; ++i) {
		A.generate();
//		A.simpleScan(input);
		err = testLU(A, false);
		if (err > maxerr) {
			savedA = A;
			maxerr = err;
		}
	}
	std::cout << endl << maxerr << " " << err << endl;
//	testLU(savedA, true);
}

void testSolver() {
	cout << "TestSolver started" << endl << endl;
	size_t s = 1024; // 2048, 16384
	TDynamicMatrix<type> A(s);
	TDynamicVector<type> b(s);
	TDynamicVector<type> x(s);
	double maxerr = 0.0;
	int64_t mintime = 0x3FFF'FFFF'FFFF'FFFF;
	int64_t time;
	double err = 0.0;
	FILE* file = fopen("C:/Users/maksi/Desktop/disk_D/cpp/matrix/matrix0.bin", "rb");
	std::chrono::steady_clock::time_point start, finish;
	for (long long step = 1; step < 100; ++step) {
		cout << "step = " << step << endl;
		maxerr = 0.0;
		mintime = 0x3FFF'FFFF'FFFF'FFFF;
	for (size_t i = 0; i < 10; ++i) {
//		cout << i << endl; //
		A.generate();
//		A.readDMatrix(file, s, 0);
		b.generate();

		start = std::chrono::steady_clock::now();
		x = A.solver(b, step);
		finish = std::chrono::steady_clock::now();
		time = std::chrono::duration_cast<std::chrono::milliseconds > (finish - start).count();

		err = (b - A * x).norm2();
		if (err != err) {
			cout << "NAN APPEARED, EXITING..." << endl;
			return;
		}
		maxerr = std::max(maxerr, err);
		mintime = std::min(mintime, time);
	}

	cout << "Maximum error is: " << maxerr << endl;
	cout << "Minimal time is: " << mintime << endl;
	cout << endl;
}
	//	fclose(file);
}

int main() {

	testSolver();
//	test();
	system("pause");

	return 0;
}
