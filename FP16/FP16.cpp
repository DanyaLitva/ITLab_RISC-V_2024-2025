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
    FP16::FP16(int _sign = 0, int _exp = 0, int _man = 0) :sign(_sign), exp(_exp), man(_man) {}
    
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
    FP16 operator-(FP16 N2) {
        FP16 temp;
        temp = *this;
        ++temp.sign;
        return (temp + N2);
    }

    //
    FP16 operator-() {
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

        temp = pow(-1, N1.sign) * (N1.man >> (MaxExp - N1.exp)) + pow(-1, N2.sign) * (N2.man >> (MaxExp - N2.exp)) + pow(-1, sign_Min) * shift + pow(-1, sign_Max) * (1 << manLength);

        if (temp < 0) Res.sign = 1;
        else Res.sign = 0;
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

        Res.exp = MaxExp;
        Res.man = (temp - (1 << manLength));

        //обработка последнего бита
        if (flag || N1.exp > N2.exp) {
            //проверяем меньшую мантиссу, которую приводили к порядку большей
            if ((N2.man % (1 << (MaxExp - N2.exp))) >= (1 << (MaxExp - N2.exp - 1)))  Res.man++;
        }
        else {
            if (flag || N1.exp < N2.exp) {
                if ((N1.man % (1 << (MaxExp - N1.exp))) >= (1 << (MaxExp - N1.exp - 1)))  Res.man++;
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
        //AssignFP16(&A, N1);
        B = N2;
        //AssignFP16(&B, N2);
    }
    else {
        A = N2;
        //AssignFP16(&A, N2);
        B = N1;
        //AssignFP16(&B, N1);
    }
    //A - субнормальное В - нормальное
    Res.sign = B.sign;
    int diff = (B.exp - shiftExp) + (shiftExp - 1);
    int temp;
    if (diff != 0) temp = pow(-1, B.sign) * ((1 << manLength) + B.man) + (pow(-1, A.sign) * (A.man >> diff));
    else temp = pow(-1, B.sign) * ((1 << manLength) + B.man) + (pow(-1, A.sign) * A.man);

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
    Res.exp = temp_exp;

    int flag = 0;
    if (temp >= (1 << (manLength + 1))) {
        flag += (temp & 1);
        temp /= 2; //и здесь один бит теряется
        Res.exp++;
    }

    //округление по отрезанным битам
    if ((temp2&(1<<(manLength))) || flag) (temp += 1);
    if (temp == (1 << manLength)) {
        Res.exp++;
        temp = 0;
    }

    Res.man = (temp);
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

















using namespace std;
int main() {
    float f1 = 0.0005f;
    float f2 = 0.00031f;
    FP16 A(f1);
    FP16 B(f2);
    A.PrintFP16_ed();
    cout << endl;
    printf("%.15f", A.GetFloat());
    cout << endl;
    B.PrintFP16_ed();
    cout << endl;
    printf("%.15f", B.GetFloat());
    cout << endl << endl << (A + B).GetFloat() << " - my" << endl;
    cout << f1 + f2 << " - need" << endl;
    cout << ((A + B).GetFloat() - (f1 + f2)) << endl;
    (A + B).PrintFP16_ed();
    cout << endl;
    cout << endl << (A * B).GetFloat() << " - my" << endl;
    cout << f1 * f2 << " - need" << endl;
    cout << ((A * B).GetFloat() - f1 * f2) << endl;
    (A * B).PrintFP16_ed();
    cout << " - my" << endl;
    (FP16(f1*f2)).PrintFP16_ed();
    cout << " - need" << endl;
    cout << endl;
   
    return 0;
}