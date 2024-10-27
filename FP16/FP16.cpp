// FP16.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdio.h>
#include <math.h>


    //распределение бит для FP16 1:5:10
const int signLength = 1;
const int expLength = 5;
const int manLength = 10;

typedef struct {
    unsigned int sign : 1; //знак           %2
    unsigned int exp : 5; //порядок         %32
    unsigned int man : 10; //мантисса       %1024
    float ref; //референс для сравнения
} FP16;

int GetSignFP16(FP16 A);
int GetExpFP16(FP16 A);
int GetManFP16(FP16 A);
float GetRefFP16(FP16 N);
//operator=
void AssignFP16(FP16* N1, FP16 N2);
//конвертация FP16 в float
void ConvertFP16tof(FP16 N, float* f);
//конвертация float в FP16(не реализовано)
void ConvertftoFP16(float f, FP16* N);
//вывод битов FP16
void PrintFP16(FP16 N);
//вывод битов FP16 с разделением знака, порядка и мантиссы
void PrintFP16_ed(FP16 N);
//вывод FP16 в виде float
void PrintFP16_f(FP16 N);
//operator+ (реализован + для нормалов)
void AddFP16(FP16 N1, FP16 N2, FP16* Res);
//operator* (реализовано для нормалов)
void MulFP16(FP16 N1, FP16 N2, FP16* Res);


int main()
{
    FP16 A,B,C;
    float f1, f2, f3;

    A.sign = 0;
    A.exp = 17;
    A.man = 25;
    B.sign = 0;
    B.exp = 16;
    B.man = 700;

    AddFP16(A, B, &C);
    ConvertFP16tof(A, &f1);
    ConvertFP16tof(B, &f2);
    ConvertFP16tof(C, &f3);
    printf("A = %f\n", f1);
    PrintFP16_ed(A);
    printf("\nB = %f\n", f2);
    PrintFP16_ed(B);
    printf("\nC = A + B = %f\n", f3);
    PrintFP16_ed(C);
    printf("\nthe result when using fleets: %f", (f1 + f2));
    printf("\nerror rate: %f", (f3 - (f1 + f2)));
    printf("\n\n");

    A.sign = 0;
    A.exp = 18;
    A.man = 253;
    B.sign = 0;
    B.exp = 20; 
    B.man = 51;  

    MulFP16(A, B, &C);
    ConvertFP16tof(A, &f1);
    ConvertFP16tof(B, &f2);
    ConvertFP16tof(C, &f3);
    printf("A = %f\n", f1);
    PrintFP16_ed(A);
    printf("\nB = %f\n", f2);
    PrintFP16_ed(B);
    printf("\nC = A * B = %f\n", f3); 
    PrintFP16_ed(C);
    printf("\nthe result when using fleets: %f", (f1 * f2));
    printf("\nerror rate: %f", (f3 - (f1 * f2)));
}







int GetSignFP16(FP16 A) { return A.sign; }

int GetExpFP16(FP16 A) { return A.exp; }

int GetManFP16(FP16 A) { return A.man; }

float GetRefFP16(FP16 N) { return N.ref; }

//operator=
void AssignFP16(FP16* N1, FP16 N2) {
    N1->sign = N2.sign;
    N1->exp = N2.exp;
    N1->man = N2.man;
    N1->ref = N2.ref;
}


void ConvertFP16tof(FP16 N, float* f) {
    int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
    if (N.sign == 1) (*f *= -1);

    float temp = 1.0f;
    float temp2 = 1.0f;
    for (int i = 0; i < manLength; ++i) {
        temp2 /= 2.0f;
        if ((N.man & (1 << (manLength - i - 1))) != 0)   temp += temp2;
    }
    //printf("%d\n", N.exp - shiftExp);
    if (N.exp >= shiftExp) *f = (float)(pow(2, (N.exp - shiftExp))) * temp;
    else *f = (float)(1/(pow(2, (shiftExp - N.exp)))) * temp;
}


void PrintFP16(FP16 N) {
    int temp;
    printf("%c", (N.sign + '0'));
    for (int i = 0; i < expLength; ++i) {
        if ((N.man & (1 <<( expLength - i-1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
    for (int i = 0; i < manLength; ++i) {
        if ((N.exp & (1 << (manLength - i-1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
}

void PrintFP16_ed(FP16 N) {
    int temp;
    printf("%c_", (N.sign + '0'));
    for (int i = 0; i < expLength; ++i) {
        if ((N.exp & (1 << (expLength - i-1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
    printf("_");
    for (int i = 0; i < manLength; ++i) {
        if ((N.man & (1 << (manLength - i-1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
}

void PrintFP16_f(FP16 N) {
    float temp;
    ConvertFP16tof(N, &temp);
    printf("%f", temp);
}


void AddFP16(FP16 N1, FP16 N2, FP16* Res) {
    int diff = abs(N1.exp - N2.exp); //разница порядков
    int MaxExp;
    if (N1.exp >= N2.exp) MaxExp = N1.exp;
    else MaxExp = N2.exp;
    int temp;

    int shift;
    if ((manLength - diff) <= 0) shift = 0;
    else shift = (1 << (manLength - diff));

    temp = (N1.man >> (MaxExp - N1.exp)) + (N2.man >> (MaxExp - N2.exp)) + shift + (1 << manLength);
    if (temp >= (1 << (manLength+1))) {
       temp = (temp / 2);
       MaxExp++;
    }
    Res->exp = MaxExp;
    Res->man = (temp - (1<<manLength));

}


void MulFP16(FP16 N1, FP16 N2, FP16* Res) {
    int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
    int temp = ((N1.man * N2.man)/(1<<manLength)) + N1.man + N2.man;

    Res->exp = ((N1.exp-shiftExp) + (N2.exp-shiftExp) + shiftExp);
    int temp2 = (1 << manLength);
    while (temp >= (1<<(manLength))) {
        temp = (temp / 2);
        temp2 /= 2;
        temp -= temp2;
        Res->exp++;
    }
    Res->man = (temp);
    Res->sign = N1.sign + N2.sign;

}