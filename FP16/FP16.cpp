// FP16.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdio.h>
#include <math.h>


    //распределение бит для FP16 1:5:10
const int signLength = 1;
const int expLength = 5;
const int manLength = 10;
/*
typedef struct {
    unsigned int sign : 1; //знак           %2
    unsigned int exp : 5; //порядок         %32
    unsigned int man : 10; //мантисса       %1024
    float ref; //референс для сравнения
} FP16;
*/
typedef union {
    struct{
        unsigned int man : 10; //мантисса       %1024
        unsigned int exp : 5; //порядок         %32
        unsigned int sign : 1; //знак           %2
    } in;
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
//operator+ с разными знаками
void AddFP16(FP16 N1, FP16 N2, FP16* Res);
//operator* (реализовано для нормалов)
void MulFP16(FP16 N1, FP16 N2, FP16* Res);

void Test();

int main()
{
    FP16 A,B,C,D;
    float f1, f2, f3;
    
    A.in.sign = 0;
    A.in.exp = 17;
    A.in.man = 25;
    
    f1 = 1025.0f;
    ConvertftoFP16(f1, &A);
    
    B.in.sign = 0;
    B.in.exp = 16;
    B.in.man = 700;
    
    f2=-1.61f;
    ConvertftoFP16(f2, &B);

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
    printf("\nthe result when using floats: %f", (f1 + f2));
    printf("\nerror rate: %f", (f3 - (f1 + f2)));
    printf("\n\n");

    


    printf("\n");
    A.in.sign = 0;
    A.in.exp = 17;
    A.in.man = 25;
    
    f1 = 6.57f;
    ConvertftoFP16(f1, &A);
    
    B.in.sign = 0;
    B.in.exp = 16;
    B.in.man = 700;
    
    f2=13.7565555f;
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







int GetSignFP16(FP16 A) { return A.in.sign; }

int GetExpFP16(FP16 A) { return A.in.exp; }

int GetManFP16(FP16 A) { return A.in.man; }

float GetRefFP16(FP16 N) { return N.ref; }

//operator=
void AssignFP16(FP16* N1, FP16 N2) {
    N1->in.sign = N2.in.sign;
    N1->in.exp = N2.in.exp;
    N1->in.man = N2.in.man;
    N1->ref = N2.ref;
}


void ConvertFP16tof(FP16 N, float* f) {
    int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
    if (N.in.sign == 1) (*f *= -1);

    float temp = 1.0f;
    float temp2 = 1.0f;
    for (int i = 0; i < manLength; ++i) {
        temp2 /= 2.0f;
        if ((N.in.man & (1 << (manLength - i - 1))) != 0)   temp += temp2;
    }
    //printf("%d\n", N.in.exp - shiftExp);
    if (N.in.exp >= shiftExp) *f = pow(-1,N.in.sign)*(float)(pow(2, (N.in.exp - shiftExp))) * temp;
    else *f = pow(-1,N.in.sign)*(float)(1/(pow(2, (shiftExp - N.in.exp)))) * temp;
}


void PrintFP16(FP16 N) {
    int temp;
    printf("%c", (N.in.sign + '0'));
    for (int i = 0; i < expLength; ++i) {
        if ((N.in.exp & (1 << (expLength - i - 1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
    for (int i = 0; i < manLength; ++i) {
        if ((N.in.man & (1 << (manLength - i - 1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
}

void PrintFP16_ed(FP16 N) {
    int temp;
    printf("%c_", (N.in.sign + '0'));
    for (int i = 0; i < expLength; ++i) {
        if ((N.in.exp & (1 << (expLength - i-1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
    printf("_");
    for (int i = 0; i < manLength; ++i) {
        if ((N.in.man & (1 << (manLength - i-1))) != 0) temp = 1; else temp = 0;
        printf("%c", (temp + '0'));
    }
}

void PrintFP16_f(FP16 N) {
    float temp;
    ConvertFP16tof(N, &temp);
    printf("%f", temp);
}

/*
void AddFP16(FP16 N1, FP16 N2, FP16* Res) {
    int diff = abs(N1.in.exp - N2.in.exp); //разница порядков
    int MaxExp;
    if (N1.in.exp >= N2.in.exp) MaxExp = N1.in.exp;
    else MaxExp = N2.in.exp;
    int temp;

    int shift;
    if ((manLength - diff) <= 0) shift = 0;
    else shift = (1 << (manLength - diff));

    temp = (N1.in.man >> (MaxExp - N1.in.exp)) + (N2.in.man >> (MaxExp - N2.in.exp)) + shift + (1 << manLength);
    if (temp >= (1 << (manLength+1))) {
       temp = (temp / 2);
       MaxExp++;
    }
    Res->in.exp = MaxExp;
    Res->in.man = (temp - (1<<manLength));
}
*/

void AddFP16(FP16 N1, FP16 N2, FP16* Res) {
    int diff = abs(N1.in.exp - N2.in.exp); //разница порядков
    int MaxExp;
    int sign_Min,sign_Max;
    if (N1.in.exp >= N2.in.exp) {
        MaxExp = N1.in.exp;
        sign_Min=N2.in.sign;
        sign_Max=N1.in.sign;
    }
    else {
        MaxExp = N2.in.exp;
        sign_Min=N1.in.sign;
        sign_Max=N2.in.sign;
    }

    int temp;
    int shift;
    if ((manLength - diff) <= 0) shift = 0;
    else shift = (1 << (manLength - diff));
    if (manLength == diff) shift = 1;

    temp = pow(-1,N1.in.sign)*(N1.in.man >> (MaxExp - N1.in.exp)) + pow(-1,N2.in.sign)*(N2.in.man >> (MaxExp - N2.in.exp)) + pow(-1,sign_Min)*shift + pow(-1,sign_Max)*(1 << manLength);

    if (temp<0) Res->in.sign = 1;
    else Res->in.sign = 0;
    temp=abs(temp);

    while((temp<(1<<manLength)) && (MaxExp>0)){
        temp*=2;
        MaxExp--;
    }

    if ((temp >= (1 << (manLength+1)) && (MaxExp<(1<<(expLength))))) {
       temp = (temp / 2);
       MaxExp++;
    }

    Res->in.exp = MaxExp;
    Res->in.man = (temp - (1<<manLength));
}


void MulFP16(FP16 N1, FP16 N2, FP16* Res) {
    int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
    int temp2 = (N1.in.man * N2.in.man); //сохранять, тут последниий бит разложить на / и %
    int temp = ((N1.in.man * N2.in.man)/(1<<manLength)) + N1.in.man + N2.in.man  + (1<<manLength);

    Res->in.exp = ((N1.in.exp-shiftExp) + (N2.in.exp-shiftExp) + shiftExp);

    int flag = 0;
    if(temp>= (1<<(manLength+1))) {
        flag = (temp & 1);
        temp/=2; //и здесь один бит теряется
        Res->in.exp++;
    }

    //округление по отрезанным битам
    if(((temp2%(1<<manLength)) >= (1<<(manLength-1))) || flag) (temp |= 1);

    Res->in.man = (temp);
    Res->in.sign = N1.in.sign + N2.in.sign;
}
//прибавить 1024 к temp
//если больше 2048 делим все на пополам


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
    
    N->in.sign = bites[0];
    N->in.exp = bites[1] * (1 << (expLength-1));
    N->in.man = 0;
    for (int i = 0; i < 4; ++i) {
        N->in.exp = N->in.exp + (bites[5+i] * (1 << (expLength - 2 - i)));
    }
    for (int i = 0; i < 10; ++i) {
        N->in.man = N->in.man + (bites[9 + i] * (1 << (manLength - 1 - i)));
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
                A.in.sign=B.in.sign=0;
                A.in.exp=B.in.exp=_exp;
                A.in.man=B.in.man=_man;
                
                MulFP16(A, B, &C);
                ConvertFP16tof(A, &f1);
                ConvertFP16tof(B, &f2);
                ConvertftoFP16(f1*f2, &D);
                //C - сумма fp16
                //D - сумма float
                PrintFP16_ed(A);
                //printf("\n");
                
                if(GetFloat(C)!=GetFloat(D)) {
                    printf("\nbreak");
                    //printf("%f %f\n",GetFloat(C),GetFloat(D));
                    flag=1;
                }
            }
        }
    }
    
}


*/
/*
//сабнормал + нормал
если эксп позволяет, то сдвигаем мантиссу
 если эксп 0, оставляем мантиссу как есть


*/
