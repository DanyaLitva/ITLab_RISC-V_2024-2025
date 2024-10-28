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
//конвертация float в FP16
void ConvertftoFP16(float f, FP16* N);
//вывод битов FP16
void PrintFP16(FP16 N);
//вывод битов FP16 с разделением знака, порядка и мантиссы
void PrintFP16_ed(FP16 N);
//вывод FP16 в виде float
void PrintFP16_f(FP16 N);
//вернуть значение FP16 в виде float
float GetFloat(FP16 N);
//operator+ (реализован + для нормалов)
void AddFP16(FP16 N1, FP16 N2, FP16* Res);
//operator+ v2 с разными знаками(черновой вариант)
void AddFP16_2(FP16 N1, FP16 N2, FP16* Res);
//operator* (реализовано для нормалов)
void MulFP16(FP16 N1, FP16 N2, FP16* Res);

void Test();

int main()
{
    
    FP16 A,B,C,D;
    float f1, f2, f3;
    
    A.sign = 0;
    A.exp = 17;
    A.man = 25;
    
    f1 = -17.05;
    ConvertftoFP16(f1, &A);
    
    B.sign = 0;
    B.exp = 16;
    B.man = 700;
    
    f2=-19.01;
    ConvertftoFP16(f2, &B);

    AddFP16_2(A, B, &C);
    ConvertFP16tof(A, &f1);
    ConvertFP16tof(B, &f2);
    ConvertFP16tof(C, &f3);
    printf("A = %f\n", f1);
    PrintFP16_ed(A);
    printf("\nB = %f\n", f2);
    PrintFP16_ed(B);
    printf("\nC = A + B = %f\n", f3);
    PrintFP16_ed(C);
    printf("\nthe result when using floats: %f", (f1 + f2));
    printf("\nerror rate: %f", (f3 - (f1 + f2)));
    printf("\n\n");

    


    printf("\n");
    A.sign = 0;
    A.exp = 17;
    A.man = 25;
    
    f1 = -17.05;
    ConvertftoFP16(f1, &A);
    
    B.sign = 0;
    B.exp = 16;
    B.man = 700;
    
    f2=19.01f;
    ConvertftoFP16(f2, &B);

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
    printf("\nthe result when using floats: %f", (f1 * f2));
    printf("\nerror rate: %f", (f3 - (f1 * f2)));
    ConvertftoFP16(f1*f2, &D);
    printf("\n");
    
    //Test();
    
    
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
    if (N.exp >= shiftExp) *f = pow(-1,N.sign)*(float)(pow(2, (N.exp - shiftExp))) * temp;
    else *f = pow(-1,N.sign)*(float)(1/(pow(2, (shiftExp - N.exp)))) * temp;
}


void PrintFP16(FP16 N) {
    int temp;
    printf("%c", (N.sign + '0'));
    for (int i = 0; i < expLength; ++i) {
        if ((N.exp & (1 << (expLength - i - 1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
    for (int i = 0; i < manLength; ++i) {
        if ((N.man & (1 << (manLength - i - 1))) != 0) temp = 1; else temp = 0;
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

void AddFP16_2(FP16 N1, FP16 N2, FP16* Res) {
    int diff = abs(N1.exp - N2.exp); //разница порядков
    int MaxExp;
    int sign_Min,sign_Max;
    if (N1.exp >= N2.exp) {
        MaxExp = N1.exp;
        sign_Min=N2.sign;
        sign_Max=N1.sign;
    }
    else {
        MaxExp = N2.exp;
        sign_Min=N1.sign;
        sign_Max=N2.sign;
    }
    int temp;
    
    int shift;
    if ((manLength - diff) <= 0) shift = 0;
    else shift = (1 << (manLength - diff));

    temp = pow(-1,N1.sign)*(N1.man >> (MaxExp - N1.exp)) + pow(-1,N2.sign)*(N2.man >> (MaxExp - N2.exp)) + pow(-1,sign_Min)*shift + pow(-1,sign_Max)*(1 << manLength);
    if (temp<0) Res->sign = 1;
    else Res->sign = 0;
    temp=abs(temp);
    while(temp<(1<<manLength)){
        temp*=2;
        MaxExp--;
    }
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

char bytes[sizeof(float)];
char bites[sizeof(float) * 8];
//float: 1:8:23
void ConvertftoFP16(float f, FP16* N){
    //f = -12.1171875f;
    for (int i = 0; i < sizeof(float) / sizeof(char); i++)
    {
        bytes[i] = *((char*)(&f) + i * sizeof(char));     
    }
    for (int i = 0; i < (sizeof(float)*8); ++i) {
        bites[31-i] = (bytes[i / 8] >> (i % 8)) & 1;
        //if(bites[31-i]!=0)    printf("%d = %d ", (31-i),bites[31-i]);
    }
    
    N->sign = bites[0];
    N->exp = bites[1] * (1 << (expLength-1));
    N->man = 0;
    for (int i = 0; i < 4; ++i) {
        N->exp = N->exp + (bites[5+i] * (1 << (expLength - 2 - i)));
    }
    for (int i = 0; i < 10; ++i) {
        N->man = N->man + (bites[9 + i] * (1 << (manLength - 1 - i)));
    }
    //PrintFP16_ed(*N);
}

float GetFloat(FP16 N){
    float temp;
    ConvertFP16tof(N, &temp);
    return temp;
}









/*
void Test(){
    //
    int flag = 0;
    FP16 A,B,C,D;
    float f1,f2,f3,f4,f5;
    //for(int _sign = 0; _sign<2;_sign++){
    for(int _exp = 20;_exp<32;_exp++){
        for(int _man = 0; _man<1024;_man++){
            if (flag!=1){
                //A.sign=B.sign=_sign;
                A.sign=B.sign=0;
                A.exp=B.exp=_exp;
                A.man=B.man=_man;
                
                AddFP16(A, B, &C);
                ConvertFP16tof(A, &f1);
                ConvertFP16tof(B, &f2);
                ConvertftoFP16(f1+f2, &D);
                //C - сумма fp16
                //D - сумма float
                PrintFP16_ed(A);
                printf("\n");
                if(GetFloat(C)!=GetFloat(D)) {
                    printf("\nbreak");
                    printf("\n%f %f\n",GetFloat(C),GetFloat(D));
                    flag=1;
                }
            }
        }
    }
    
}

*/
