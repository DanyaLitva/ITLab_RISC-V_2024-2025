// FP16.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <stdio.h>
#include <math.h>


    //распределение бит для FP16 1:5:10
const int signLength = 1;
const int expLength = 5;
const int manLength = 10;
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
//проверка числа на сабнормальность
int IsSubnormal(FP16 N);
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
//сумма нормального и сабнормального
void AddFP16_N_S(FP16 N1, FP16 N2, FP16* Res);
//сумма сабнормальных
void AddFP16_S_S(FP16 N1, FP16 N2, FP16* Res);
//operator*
void MulFP16(FP16 N1, FP16 N2, FP16* Res);
//произведение нормального и субнормального
void MulFP16_N_S(FP16 N1, FP16 N2, FP16* Res);
//произведение субнормальных
void MulFP16_S_S(FP16 N1, FP16 N2, FP16* Res);

void Test();

//void GetX0(FP16 D);

int main()
{
    FP16 A,B,C,D;
    float f1=0, f2=0, f3=0;
    
    A.in.sign = 0;
    A.in.exp = 2;
    A.in.man = 512;
    
    //f1 = 0.000003;
    //ConvertftoFP16(f1, &A);
    
    B.in.sign = 0;
    B.in.exp = 0;
    B.in.man = 180;
    
    //f2= 0.000005;
    //ConvertftoFP16(f2, &B);

    AddFP16(A, B, &C);
    ConvertFP16tof(A, &f1);
    ConvertFP16tof(B, &f2);
    ConvertFP16tof(C, &f3);
    printf("A = %.15f\n", f1);
    PrintFP16_ed(A);
    printf("\nB = %.15f\n", f2);
    PrintFP16_ed(B);
    printf("\nC = A + B = %.15f\n", f3);
    PrintFP16_ed(C);
    printf("\nthe result when using floats: %.15f", (f1 + f2));
    printf("\nerror rate: %.15f", (f3 - (f1 + f2)));
    printf("\n\n");

    

    printf("\n");
    A.in.sign = 0;
    A.in.exp =13;
    A.in.man = 164;
    
    f1 = 0.55555;
    ConvertftoFP16(f1, &A);
    
    B.in.sign = 0;
    B.in.exp = 0;
    B.in.man = 511;
    
    f2 = 0.000030517578125;
    //f2 = 0.000030457973480;
    ConvertftoFP16(f2, &B);
    
    MulFP16(A, B, &C);
    ConvertFP16tof(A, &f1);
    ConvertFP16tof(B, &f2);
    ConvertFP16tof(C, &f3);
    printf("A = %.15f\n", f1);
    PrintFP16_ed(A);
    printf("\nB = %.15f\n", f2);
    PrintFP16_ed(B);
    printf("\nC = A * B = %.15f\n", f3);
    PrintFP16_ed(C);
    printf("\nthe result when using floats: %.15f", (f1 * f2));
    printf("\nerror rate: %.15f", (f3 - (f1 * f2)));
    ConvertftoFP16(f1*f2, &D);
    printf("\n");

    /*
    FP16 temp;
    temp.in.sign=0;
    temp.in.exp=15;
    temp.in.man=20;
    GetX0(temp);
     */
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


void ConvertFP16tof(FP16 A, float* f) {
    unsigned int t = 0;
    *f = 0;
    int shiftExp = (1 << (expLength - 1)) - 1;
    t |= (((unsigned int)A.in.sign) << 31);
    if (IsSubnormal(A) == 0) {
        if (((int)(A.in.exp) - shiftExp) > 0) {
            t |= 1 << 30;
            t |= (((A.in.exp & ((1 << (expLength - 1)) - 1)) << 23));
        }
        else {
            if (A.in.exp != 0) {
                for (int i = 0; i < 3; ++i) {
                    t += (1 << (27 + i));
                }
                t += ((A.in.exp & ((1 << (expLength - 1)) - 1)) << 23);
            }
        }
        t |= (A.in.man) << 13;
        *(f) = *((float*)(&t));
    }
    else {
        for (int i = 0; i < manLength; ++i) {
            if((A.in.man & (1 << (manLength - 1 - i))) != 0) *f += (float)(1/(pow(2,15+i)));
        }
    }
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

char bytes[sizeof(float)];
char bites[sizeof(float) * 8];
//float: 1:8:23
void ConvertftoFP16(float f, FP16* N){
    int shiftExp = (1<<(expLength-1)) - 1;
    for (int i = 0; i < sizeof(float) / sizeof(char); i++)
    {
        bytes[i] = *((char*)(&f) + i * sizeof(char));
    }
    for (int i = 0; i < (sizeof(float)*8); ++i) {
        bites[31-i] = (bytes[i / 8] >> (i % 8)) & 1;
    }
    N->in.sign = bites[0];
    int manf = 0;
    for (int i = 0; i < 10; ++i) {
        manf += (bites[9 + i] * (1 << (manLength - 1 - i)));
    }
    
    int expf=0;
    for (int i = 0; i < 8; ++i) {
        expf+=(bites[1 + i])*(1<<(7 - i));
    }
    expf -= ((1<<(8-1)) - 1); //смещение float
    if (expf > 16) expf=16;
    
    if (expf==-15){
        manf=manf>>1;
        manf+=(1<<(manLength-1));
    }
    
    if (expf < -15){
        manf = (manf>>(-shiftExp - expf));
        expf = -15;
    }
    
    N->in.man = manf;
    N->in.exp = (expf + shiftExp);
}



float GetFloat(FP16 N){
    float temp;
    ConvertFP16tof(N, &temp);
    return temp;
}



int IsNull(FP16 N) {
    if ((N.in.exp == 0) && (N.in.man == 0)) return 1;
    return 0;
}

int IsSubnormal(FP16 N) {
    if ((N.in.exp != 0) || (N.in.man == 0)) return 0;
    return 1;
}


void AddFP16_N_S(FP16 N1, FP16 N2, FP16* Res) {
    FP16 A, B;
    int shiftExp = (1 << (expLength - 1)) - 1;
    if (IsSubnormal(N1)) {
        AssignFP16(&A,N1);
        AssignFP16(&B,N2);
    }
    else {
        AssignFP16(&A, N2);
        AssignFP16(&B, N1);
    }
    //A - субнормальное В - нормальное
    Res->in.sign = B.in.sign;
    int diff = (B.in.exp-shiftExp) + (shiftExp - 1);
    int temp;
    if (diff != 0) temp = pow(-1, B.in.sign) * ((1 << manLength) + B.in.man) + (pow(-1, A.in.sign) * (A.in.man >> diff));
    else temp = pow(-1,B.in.sign)*((1 << manLength) + B.in.man) + (pow(-1,A.in.sign) * A.in.man);

    Res->in.exp = B.in.exp;
    if ((temp >= (1 << (manLength + 1)) && (Res->in.exp < ((1 << (expLength))-1)))) {
        temp = (temp / 2);
        Res->in.exp++;
    }
    if ((temp < (1 << manLength)) && (Res->in.exp > 0)) {
        temp *= 2;
        Res->in.exp--;
    }
    Res->in.man = (temp - (1 << manLength));
}

void AddFP16_S_S(FP16 N1, FP16 N2, FP16* Res) {
    if (N1.in.man > N2.in.man) Res->in.sign = N1.in.sign;
    else Res->in.sign = N2.in.sign;

    Res->in.exp = 0;
    int temp = (pow(-1, N1.in.sign) * N1.in.man) + (pow(-1, N2.in.sign) * N2.in.man);
    if (temp >= (1 << manLength)) {
        temp -=1024;
        Res->in.exp++;
    }
    Res->in.man = temp;
}

void AddFP16(FP16 N1, FP16 N2, FP16* Res) {
    if ((IsSubnormal(N1) + IsSubnormal(N2)) == 0) {
        if (IsNull(N1) + IsNull(N2) == 0) {
            int diff = abs(N1.in.exp - N2.in.exp); //разница порядков
            int MaxExp;
            int sign_Min, sign_Max;
            int flag = 0; //для остатка
            if (N1.in.exp >= N2.in.exp) {
                MaxExp = N1.in.exp;
                sign_Min = N2.in.sign;
                sign_Max = N1.in.sign;
            }
            else {
                MaxExp = N2.in.exp;
                sign_Min = N1.in.sign;
                sign_Max = N2.in.sign;
            }

            int temp;
            int shift;
            if ((manLength - diff) <= 0) shift = 0;
            else shift = (1 << (manLength - diff));
            if (manLength == diff) shift = 1;

            temp = pow(-1, N1.in.sign) * (N1.in.man >> (MaxExp - N1.in.exp)) + pow(-1, N2.in.sign) * (N2.in.man >> (MaxExp - N2.in.exp)) + pow(-1, sign_Min) * shift + pow(-1, sign_Max) * (1 << manLength);

            if (temp < 0) Res->in.sign = 1;
            else Res->in.sign = 0;
            temp = abs(temp);

            if ((temp < (1 << manLength)) && (MaxExp > 0)) {
                temp *= 2;
                MaxExp--;
            }

            if ((temp >= (1 << (manLength + 1)) && (MaxExp < (1 << (expLength))))) {
                flag += temp & 1;
                temp = (temp / 2);
                MaxExp++;
            }

            Res->in.exp = MaxExp;
            Res->in.man = (temp - (1 << manLength));

            //обработка последнего бита
            if (flag || N1.in.exp > N2.in.exp) {
                //проверяем меньшую мантиссу, которую приводили к порядку большей
                if ((N2.in.man % (1 << (MaxExp - N2.in.exp))) >= (1 << (MaxExp - N2.in.exp - 1)))  Res->in.man++;
            }
            else {
                if (flag || N1.in.exp < N2.in.exp) {
                    if ((N1.in.man % (1 << (MaxExp - N1.in.exp))) >= (1 << (MaxExp - N1.in.exp - 1)))  Res->in.man++;
                }
            }
        }
        else {
            Res->in.exp = N1.in.exp + N2.in.exp;
            Res->in.man = N1.in.man + N2.in.man;
            if (IsNull(N1)) Res->in.sign = N2.in.sign;
            else Res->in.sign = N1.in.sign;
        }
    }
    else {
        if ((IsSubnormal(N1) + IsSubnormal(N2)) == 1) AddFP16_N_S(N1, N2, Res);
        else AddFP16_S_S(N1, N2, Res);
    }
}



void MulFP16_N_S(FP16 N1, FP16 N2, FP16* Res) {
    FP16 A, B;
    int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
    if (IsSubnormal(N1)) {
        AssignFP16(&A, N1);
        AssignFP16(&B, N2);
    }
    else {
        AssignFP16(&A, N2);
        AssignFP16(&B, N1);
    }
    //A - субнормальное В - нормальное
    int flag = 0;
    int temp2 = A.in.man * B.in.man; //сохранять, тут последниий бит разложить на / и %
    int temp = ((A.in.man * B.in.man) / (1 << manLength)) + A.in.man;
    int temp_exp = B.in.exp - shiftExp + 1;
    if (temp_exp<= 0 )flag += temp & 1;
    while (temp_exp <= 0) {
        temp_exp++;
        temp /= 2;
    }
    Res->in.exp = temp_exp;

    if (Res->in.exp == 0) {
        if (temp >= (1 << manLength)) {
            temp -= (1 << manLength);
            Res->in.exp++;
        }
        Res->in.man = temp;
    }
    else {
        if (Res->in.exp == 1) {
            if (temp < (1 << manLength)) {
                temp += (1 << manLength);
                Res->in.exp--;
            }
            if (temp >= (1 << (manLength + 1))) {
                flag += temp & 1;
                temp = (temp / 2);
                Res->in.exp++;
            }
        }
        else {
            if (temp < (1 << manLength)) {
                temp *= 2;
                Res->in.exp--;
            }

            if (temp >= (1 << (manLength + 1))) {
                flag += temp & 1;
                temp = (temp / 2);
                Res->in.exp++;
            }
        }
        Res->in.man = temp - (1 << (manLength));
    }
    
    //округление по отрезанным битам
    if (((temp2 % (1 << manLength)) >= (1 << (manLength - 1))) || flag) (Res->in.man += 1);
    Res->in.sign = N1.in.sign ^ N2.in.sign;
    if (IsNull(N1) || IsNull(N2)) Res->in.man = (Res->in.exp = 0);
}


void MulFP16_S_S(FP16 N1, FP16 N2, FP16* Res) {
    Res->in.sign = N1.in.sign ^ N2.in.sign;
    Res->in.man = Res->in.exp = 0;
}


void MulFP16(FP16 N1, FP16 N2, FP16* Res) {
    if ((IsSubnormal(N1) + IsSubnormal(N2)) == 0) {
        int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
        int temp2 = (N1.in.man * N2.in.man); //сохранять, тут последниий бит разложить на / и %
        int temp = ((N1.in.man * N2.in.man) / (1 << manLength)) + N1.in.man + N2.in.man + (1 << manLength);

         int temp_exp = ((N1.in.exp - shiftExp) + (N2.in.exp - shiftExp) + shiftExp);
        if (temp_exp <=0){
            temp = temp >> (-1*temp_exp + 1);
            temp_exp=0;
            //получается сабнормал
            //переход к нормалу
            if (temp>=(1<<manLength)){
                temp-=(1<<manLength);
                temp_exp++;
            }
        }
         Res->in.exp = temp_exp;

        int flag = 0;
        if (temp >= (1 << (manLength + 1))) {
            flag += (temp & 1);
            temp /= 2; //и здесь один бит теряется
            Res->in.exp++;
        }

        //округление по отрезанным битам
        if (((temp2 % (1 << manLength)) >= (1 << (manLength - 1))) || flag) (temp += 1);

        Res->in.man = (temp);
        Res->in.sign = N1.in.sign + N2.in.sign;
        if (IsNull(N1) || IsNull(N2)) Res->in.man = (Res->in.exp = 0);
    }
    else {
        if ((IsSubnormal(N1) + IsSubnormal(N2)) == 1) MulFP16_N_S(N1, N2, Res);
        else MulFP16_S_S(N1, N2, Res);
    }
}

/*
void GetX0(FP16 D){
    FP16 x0;
    x0.in.sign=0;
    x0.in.exp=20;
    x0.in.man=512;
    FP16 A;
    FP16 t_48;
    FP16 t_1_17;
    FP16 t_32;
    FP16 temp;
    ConvertftoFP16(48, &t_48);
    ConvertftoFP16((1/17), &t_1_17);
    ConvertftoFP16(32, &t_32);
    MulFP16(t_48, t_1_17, &temp);
    AddFP16(x0, temp, &x0);
    MulFP16(t_32, t_1_17, &temp);
    MulFP16(temp, D, &temp);
    AddFP16(x0, temp, &x0);
    PrintFP16_ed(x0);
    printf("\n");
    PrintFP16_f(x0);
    printf("\n");
}
*/
