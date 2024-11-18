// FP16.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>

//распределение бит для FP16 1:5:10
const int manLength = 10;
const int expLength = 5;
const int signLength = 1;
const int shiftExp = (1 << (expLength - 1)) - 1;

class FP16 {

public:
    unsigned int man : manLength; //мантисса       %1024
    unsigned int exp : expLength; //порядок         %32
    unsigned int sign : signLength; //знак           %2


    //
    FP16(int _sign = 0, int _exp = 0, int _man = 0) :sign(_sign), exp(_exp), man(_man) {}
    
    //
    float ConvertFP16tof() const {
        float f = 0;
        unsigned int t = 0;

        t |= (((unsigned int)sign) << 31);
        if (IsSubnormal() == 0) { //если не сабнормал
            if (((int)(exp)-shiftExp) > 0) {
                t |= 1 << 30;
                t |= (((exp & ((1 << (expLength - 1)) - 1)) << 23));
            }
            else {
                if (exp != 0) {
                    for (int i = 0; i < 3; ++i) {
                        t += (1 << (27 + i));
                    }
                    t += ((exp & ((1 << (expLength - 1)) - 1)) << 23);
                }
            }
            t |= (man) << 13;
            f = *((float*)(&t));
        }
        else { //если сабнормал
            for (int i = 0; i < manLength; ++i) {
                if ((man & (1 << (manLength - 1 - i))) != 0) f += (float)(1 / (float)(1<<(15+i)));
            }
            if (sign == 1) f *= -1;
        }
        return f;
    }

    //
    float GetFloat() const { return ConvertFP16tof(); }

    //
    void ConvertftoFP16(const float f) {
        char bytes[sizeof(float)];
        char bites[sizeof(float) * 8];

        for (int i = 0; i < sizeof(float) / sizeof(char); i++)
        {
            bytes[i] = *((char*)(&f) + i * sizeof(char));
        }
        for (int i = 0; i < (sizeof(float) * 8); ++i) {
            bites[31 - i] = (bytes[i / 8] >> (i % 8)) & 1;
        }

        sign = bites[0];

        int manf = 0;
        for (int i = 0; i < 10; ++i) {
            manf += (bites[9 + i] * (1 << (manLength - 1 - i)));
        }

        int expf = 0;
        for (int i = 0; i < 8; ++i) {
            expf += (bites[1 + i]) * (1 << (7 - i));
        }
        expf -= ((1 << (8 - 1)) - 1); //смещение float

        if (expf > 16) expf = 16;

        if (expf == -15) {
            manf = manf >> 1;
            manf += (1 << (manLength - 1));
        }

        if (expf < -15) {
            manf = (manf >> (-shiftExp - expf));
            expf = -15;
        }

        man = manf;
        exp = (expf + shiftExp);
    }

    //
    FP16(const float f) {
        ConvertftoFP16(f);
    }

    //
    void PrintFP16() const {
        int temp;
        printf("%c", (sign + '0'));
        for (int i = 0; i < expLength; ++i) {
            if ((exp & (1 << (expLength - i - 1))) != 0) temp = 1; else temp = 0;
            printf("%c", (temp + '0'));
        }
        for (int i = 0; i < manLength; ++i) {
            if ((man & (1 << (manLength - i - 1))) != 0) temp = 1; else temp = 0;
            printf("%c", (temp + '0'));
        }
    }

    //
    void PrintFP16_ed() const{
        int temp;
        printf("%c_", (sign + '0'));
        for (int i = 0; i < expLength; ++i) {
            if ((exp & (1 << (expLength - i - 1))) != 0) temp = 1; else temp = 0;
            printf("%c", (temp + '0'));
        }
        printf("_");
        for (int i = 0; i < manLength; ++i) {
            if ((man & (1 << (manLength - i - 1))) != 0) temp = 1; else temp = 0;
            printf("%c", (temp + '0'));
        }
    }

    //
    void PrintFP16_f() const{
        float f = this->ConvertFP16tof();
        printf("%f", f);
    }

    //
    bool IsNull() const noexcept{
        if ((exp == 0) && (man == 0)) return true;
        return false;
    }

    //
    bool IsSubnormal() const noexcept{
        if ((exp != 0) || (man == 0)) return false;
        return true;
    }

    //
    FP16& operator=(FP16 N) {
        sign = N.sign;
        exp = N.exp;
        man = N.man;
        return *this;
    }

    //
    FP16 operator+(FP16 N2) {
        if ((IsSubnormal() + N2.IsSubnormal()) == 0) {
            return AddFP16_N_N(*this, N2);
        }
        else {
            if ((IsSubnormal() + N2.IsSubnormal()) == 1) return AddFP16_N_S(*this, N2);
            else return AddFP16_S_S(*this, N2);
        }
    }

    //
    FP16 operator-(FP16 N2) const {
        FP16 temp(*this);
        ++N2.sign;
        return (temp + N2);
    }

    //
    FP16 operator-() const {
        FP16 temp;
        temp = *this;
        ++temp.sign;
        return (temp);
    }
    
    //
    FP16 operator*(FP16 N2) {
        if ((IsSubnormal() + N2.IsSubnormal()) == 0) {
            return MulFP16_N_N(*this, N2);
        }
        else {
            if ((IsSubnormal() + N2.IsSubnormal()) == 1) return MulFP16_N_S(*this, N2);
            else return MulFP16_S_S(*this, N2);
        }
    }

    FP16 operator/(FP16 N2) {
        FP16 Res;
        FP16 N, D;
        FP16 a, b, z, t,s;
        N = *this;
        D = N2;
        

        //D на отрезке от 0.5 до 1
        while ((D.exp - shiftExp) < -1) {
            D.exp++;
            N.exp++;
        }
        while ((D.exp - shiftExp) > -1) {
            D.exp--;
            N.exp--;
        }
        
        FP16 X = FP16(48.f / 17.f) - (FP16(32.f / 17.f) * D);
        //FP16 X = FP16(1.f/17.f)*(FP16(48.f) - (FP16(32.f) * D));
        
        //X = (FP16(2.0f) * X) - (D * X * X);
        X = X * ((FP16(2.0f)) - (D * X));

        a = (D - FP16(0.5f)) * X;
        b = FP16(0.5f) * X;
        s = a + b;
        z = s - a;
        t = b - z;

        X = (FP16(2.0f) * X) - (s * X) - (t * X);


         //добавить округление при 0.5 в зависимости от последнего бита
         //если пред четное, то не округляется вверх,
        Res = N * X;
        return Res;
    }


    float GetLastBit() {
        float temp = pow(2, abs((int)exp - shiftExp - manLength));
        if (((int)exp - shiftExp - manLength) < 0) temp = 1 / temp;
        return temp;
    }



    //
protected:
    FP16 AddFP16_N_N(FP16 N1, FP16 N2);
    FP16 AddFP16_N_S(FP16 N1, FP16 N2);
    FP16 AddFP16_S_S(FP16 N1, FP16 N2);
    FP16 MulFP16_N_N(FP16 N1, FP16 N2);
    FP16 MulFP16_N_S(FP16 N1, FP16 N2);
    FP16 MulFP16_S_S(FP16 N1, FP16 N2);

    //

};

























































FP16 FP16::AddFP16_N_N(FP16 N1, FP16 N2) {
    FP16 Res;
    if (N1.IsNull() + N2.IsNull() == 0) {
        int diff = abs(int(N1.exp - N2.exp)); //разница порядков
        int MaxExp;
        int sign_Min, sign_Max;
        int flag = 0; //для остатка
        if (N1.exp >= N2.exp) {
            MaxExp = N1.exp;
            sign_Min = N2.sign;
            sign_Max = N1.sign;
        }
        else {
            MaxExp = N2.exp;
            sign_Min = N1.sign;
            sign_Max = N2.sign;
        }

        int temp;
        int shift;
        if ((manLength - diff) <= 0) shift = 0;
        else shift = (1 << (manLength - diff));
        if (manLength == diff) shift = 1;

        temp = (pow(-1, N1.sign) * (N1.man >> (MaxExp - N1.exp))) + (pow(-1, N2.sign) * (N2.man >> (MaxExp - N2.exp))) + (pow(-1, sign_Min) * shift) + (pow(-1, sign_Max) * (1 << manLength));


        if (N1.exp > N2.exp) {
            flag = 1 & (N2.man >> (MaxExp - N2.exp - 1));
        }
        if (N1.exp < N2.exp) {
            flag = 1 & (N1.man >> (MaxExp - N1.exp - 1));
        }
        if (temp < 0) Res.sign = 1;
        else Res.sign = 0;
        temp = abs(temp);

        if (temp == 0) MaxExp = 0;

        while ((temp < (1 << manLength)) && (MaxExp > 0)) {
            temp *= 2;
            MaxExp--;
        }

        if ((temp >= (1 << (manLength + 1)) && (MaxExp < (1 << (expLength))))) {
            flag += temp & 1;
            temp = (temp / 2);
            MaxExp++;
        }

        Res.exp = MaxExp;
        Res.man = (temp - (1 << manLength));
        //обработка последнего бита
        if (flag && (N1.exp > N2.exp)) {
            //проверяем меньшую мантиссу, которую приводили к порядку большей
            if ((N2.man % (1 << (MaxExp - N2.exp))) >= (1 << (MaxExp - N2.exp - 1)))  
                Res.man = Res.man + pow(-1,N2.sign);
        }
        else {
            if (flag && (N1.exp <= N2.exp)) {
                if ((N1.man % (1 << (MaxExp - N1.exp))) >= (1 << (MaxExp - N1.exp - 1)))  
                    Res.man = Res.man + pow(-1, N1.sign);
            }
        }
    }
    else {
        Res.exp = N1.exp + N2.exp;
        Res.man = N1.man + N2.man;
        Res.sign = N1.sign;
        if (N1.IsNull()) Res.sign = N2.sign;
        if (N2.IsNull()) Res.sign = N1.sign;
    }
    return Res;
}

//
FP16 FP16::AddFP16_N_S(FP16 N1, FP16 N2) {
    FP16 A, B, Res;
    int shiftExp = (1 << (expLength - 1)) - 1;
    if (N1.IsSubnormal()) {
        A = N1;
        B = N2;
    }
    else {
        A = N2;
        B = N1;
    }
    //A - субнормальное В - нормальное
    Res.sign = B.sign;
    int diff = (B.exp - shiftExp) + (shiftExp - 1);
    int temp;

    if (diff != 0) temp = pow(-1, B.sign) * ((1 << manLength) + B.man) + (pow(-1, A.sign) * (A.man >> diff));
    else temp = pow(-1, B.sign) * ((1 << manLength) + B.man) + (pow(-1, A.sign) * A.man);

    //int flag = 0;

    Res.exp = B.exp;
    if ((temp >= (1 << (manLength + 1)) && (Res.exp < ((1 << (expLength)) - 1)))) {
        temp = (temp / 2);
        Res.exp++;
    }
    if ((temp < (1 << manLength)) && (Res.exp > 0)) {
        temp *= 2;
        Res.exp--;
    }
    Res.man = (temp - (1 << manLength));
    return Res;
}

//
FP16 FP16::AddFP16_S_S(FP16 N1, FP16 N2) {
    FP16 Res;

    if (N1.man > N2.man) Res.sign = N1.sign;
    else Res.sign = N2.sign;

    Res.exp = 0;
    unsigned int temp = (pow(-1, N1.sign) * N1.man) + (pow(-1, N2.sign) * N2.man);
    if (temp >= (1 << manLength)) {
        temp -= 1024;
        Res.exp++;
    }
    Res.man = temp;
    return Res;
}


//
FP16 FP16::MulFP16_N_N(FP16 N1, FP16 N2) {
    FP16 Res;
    unsigned int shiftExp = (1 << (expLength - 1)) - 1; //смещение   -15 ... 16
    unsigned int temp2 = (N1.man * N2.man); //сохранять, тут последниий бит разложить на / и %
    unsigned int temp = ((N1.man * N2.man) / (1 << manLength)) + N1.man + N2.man + (1 << manLength);

    int temp_exp = ((N1.exp - shiftExp) + (N2.exp - shiftExp) + shiftExp);
    if (temp_exp <= 0) {
        temp = temp >> (-1 * temp_exp + 1);
        temp_exp = 0;
        //получается сабнормал
        // 
        //переход к нормалу
        if (temp >= (1 << manLength)) {
            temp -= (1 << manLength);
            temp_exp++;
        }
    }
    

    int flag = 0;
    if (temp >= (1 << (manLength + 1))) {
        flag += (temp & 1);
        temp /= 2; //и здесь один бит теряется
        temp_exp++;
    }

    //округление по отрезанным битам
    //теперь к ближайшему четному
    if (((temp2 & (1 << (manLength))) || flag) && ((temp_exp & 1) == 1))
        (temp += 1);
    if (temp == (1 << (manLength+1)) || ((temp == (1 << (manLength))) && (temp_exp == 0))) {
        temp_exp++;
        temp = 0;
    }
    Res.exp = temp_exp;
    Res.man = (temp);
    if (temp_exp>=(1<<expLength)){
        Res.exp = (1<<manLength)-1;
        Res.man = 0;
    }
    Res.sign = N1.sign + N2.sign;
    if (N1.IsNull() || N2.IsNull()) Res.man = (Res.exp = 0);
    return Res;
}

//
FP16 FP16::MulFP16_N_S(FP16 N1, FP16 N2) {
    FP16 A, B, Res;
    // Определяем, какое число нормальное, а какое субнормальное
    if (N1.IsSubnormal()) {
        A = N1;
        B = N2; 
    }
    else {
        A = N2; 
        B = N1; 
    }
    // A - субнормальное, Б - нормальное
    Res.sign = A.sign ^ B.sign; // Определяем знак результата
    int temp_exp = B.exp - shiftExp + 1;
    unsigned int mul_man = A.man * B.man;
    unsigned int temp = (mul_man >> manLength) + A.man;

    while (temp_exp < 0) {
        temp /= 2;
        temp_exp++;
    }
    
    //проверка на ноль
    if (temp == 0) {
        Res.exp = Res.man = 0; return Res;
    }

    if (temp_exp == 0) {
        if (temp >= (1 << manLength)) {
            temp_exp++;
        }
    }

    if (temp >= (1 << (manLength + 1))) {
        temp /= 2;
        temp_exp++;
    }

    while ((temp_exp > 0) && (temp < (1 << manLength))) {
        if (temp_exp == 1) temp_exp--;
        else {
            temp_exp--;
            temp *= 2;
        }
    }


    Res.exp = temp_exp;
    if (temp_exp == 0) {
        Res.man = temp;
    }
    else Res.man = temp - (1 << manLength);

    if (N1.IsNull() || N2.IsNull()) Res.exp = Res.man = 0;
    return Res;
}

//
FP16 FP16::MulFP16_S_S(FP16 N1, FP16 N2) {
    FP16 Res;
    Res.sign = N1.sign ^ N2.sign;
    Res.man = Res.exp = 0;
    return Res;
}











void Test_Mult();
void Test_Add();
void Test_Sub();


/*1 - A   0.124695 - B
0.875977 - my   0.875305 - need
0_01110_1100000010 - my
0_01110_1100000000 - need

1 - A   0.234985 - B
0.765625 - my   0.765015 - need
0_01110_1000100000 - my
0_01110_1000011110 - need

1 - A   0.235962 - B
0.764648 - my   0.764038 - need
0_01110_1000011110 - my
0_01110_1000011100 - need

1 - A   0.236938 - B
0.763672 - my   0.763062 - need
0_01110_1000011100 - my
0_01110_1000011010 - need

1 - A   0.237915 - B
0.762695 - my   0.762085 - need
0_01110_1000011010 - my
0_01110_1000011000 - need*/

using namespace std;
int main() {
    //Test_Mult();

    float f1 = 1;
    float f2 = 0.235962;
    
    FP16 A(f1);
    FP16 B(f2);
    f1 = A.GetFloat();
    f2 = B.GetFloat();
    A.PrintFP16_ed();
    cout << endl;
    printf("A = %.15f", A.GetFloat());
    cout << endl;
    B.PrintFP16_ed();
    printf("\nB = %.15f", B.GetFloat());
    cout << endl;

    

    cout << endl;
    cout << "Operator +";
    cout << endl << (A + B).GetFloat() << " - my" << endl;
    cout << f1 + f2 << " - need" << endl;
    cout << ((A + B).GetFloat() - (f1 + f2))<< " - calculation error" << endl;
    cout << (A + B).GetLastBit()<<"- the value of the last bit" << endl;
    (A + B).PrintFP16_ed();
    cout << " - my" << endl;
    (FP16(f1 + f2)).PrintFP16_ed();
    cout << " - need" << endl;
    cout << endl;


    cout << endl;
    cout << "Operator -";
    cout << endl << (A - B).GetFloat() << " - my" << endl;
    cout << f1 - f2 << " - need" << endl;
    cout << ((A - B).GetFloat() - (f1 - f2)) << " - calculation error" << endl;
    cout << (A - B).GetLastBit() << "- the value of the last bit" << endl;
    (A - B).PrintFP16_ed();
    cout << " - my" << endl;
    (FP16(f1 - f2)).PrintFP16_ed();
    cout << " - need" << endl;
    cout << endl;


    cout << endl;
    cout << "Operator *";
    cout << endl << (A * B).GetFloat() << " - my" << endl;
    cout << f1 * f2 << " - need" << endl;
    cout << ((A * B).GetFloat() - (f1 * f2)) << " - calculation error" << endl;
    cout << (A * B).GetLastBit() << "- the value of the last bit" << endl;
    (A * B).PrintFP16_ed();
    cout << " - my" << endl;
    (FP16(f1 * f2)).PrintFP16_ed();
    cout << " - need" << endl;
    cout << endl;

    
    cout << endl;
    cout << "Operator /";
    cout << endl << (A / B).GetFloat() << " - my" << endl;
    cout << f1 / f2 << " - need" << endl;
    cout << ((A / B).GetFloat() - (f1 / f2)) << " - calculation error" << endl;
    cout << (A / B).GetLastBit() << "- the value of the last bit" << endl;
    (A / B).PrintFP16_ed();
    cout << " - my" << endl;
    (FP16(f1 / f2)).PrintFP16_ed();
    cout << " - need" << endl;
    cout << endl;


    
    return 0;
}


void Test_Mult() {
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 15; exp1 < (1 << expLength); ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength); ++exp2) {
                        for (size_t man2 = 900; man2 < (1 << manLength); ++man2) {
                            FP16 A(sign1, exp1, man1);
                            FP16 B(sign2, exp2, man2);
                            if (((A * B).GetFloat() - (A.GetFloat() * B.GetFloat())) > ((A * B).GetLastBit())) {
                                cout << A.GetFloat() << " - A\t" << B.GetFloat() <<" - B" << endl;
                                cout << (A * B).GetFloat() << " - my\t" << (A.GetFloat() * B.GetFloat()) <<" - need" << endl;
                                (A * B).PrintFP16_ed();
                                cout << " - my" << endl;
                                FP16((A.GetFloat() * B.GetFloat())).PrintFP16_ed();
                                cout << " - need" << endl;
                                cout << endl;
                            }
                        }
                    }
                }
            }
            cout << "Exp = " << exp1 << endl;
        }
    }
}



void Test_Add() {
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 0; exp1 < (1 << expLength); ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength); ++exp2) {
                        for (size_t man2 = 900; man2 < (1 << manLength); ++man2) {
                            FP16 A(sign1, exp1, man1);
                            FP16 B(sign2, exp2, man2);
                            if (((A + B).GetFloat() - (A.GetFloat() + B.GetFloat())) > ((A + B).GetLastBit())) {
                                cout << A.GetFloat() << " - A\t" << B.GetFloat() << " - B" << endl;
                                cout << (A + B).GetFloat() << " - my\t" << (A.GetFloat() + B.GetFloat()) << " - need" << endl;
                                (A + B).PrintFP16_ed();
                                cout << " - my" << endl;
                                FP16((A.GetFloat() + B.GetFloat())).PrintFP16_ed();
                                cout << " - need" << endl;
                                cout << endl;
                            }
                        }
                    }
                }
            }
            cout << "Exp = " << exp1 << endl;
        }
    }
}



void Test_Sub() {
    for (size_t sign1 = 0; sign1 < 1; ++sign1) {
        for (size_t exp1 = 16; exp1 < (1 << expLength); ++exp1) {
            for (size_t man1 = 0; man1 < (1 << manLength); ++man1) {
                for (size_t sign2 = 0; sign2 < 1; ++sign2) {
                    for (size_t exp2 = 0; exp2 < (1 << expLength); ++exp2) {
                        for (size_t man2 = 900; man2 < (1 << manLength); ++man2) {
                            FP16 A(sign1, exp1, man1);
                            FP16 B(sign2, exp2, man2);
                            if (((A - B).GetFloat() - (A.GetFloat() - B.GetFloat())) > ((A - B).GetLastBit())) {
                                cout << A.GetFloat() << " - A\t" << B.GetFloat() << " - B" << endl;
                                cout << (A - B).GetFloat() << " - my\t" << (A.GetFloat() - B.GetFloat()) << " - need" << endl;
                                (A - B).PrintFP16_ed();
                                cout << " - my" << endl;
                                FP16((A.GetFloat() - B.GetFloat())).PrintFP16_ed();
                                cout << " - need" << endl;
                                cout << endl;
                            }
                        }
                    }
                }
            }
            cout << "Exp = " << exp1 << endl;
        }
    }
}
