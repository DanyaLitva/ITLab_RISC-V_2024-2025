#include <iostream>

using float16_t = _Float16;

//распределение бит для FP16 1:5:10
const int manLength = 10; //биты для мантиссы 0-1023
const int expLength = 5; //биты для экспоненты 0-31
const int signLength = 1; //биты для знаков 0-1
const int shiftExp = (1 << (expLength - 1)) - 1;  //смещение мантиссы на 15

class FP16 {
private:
    unsigned int man : manLength; //мантисса       %1024
    unsigned int exp : expLength; //порядок         %32
    unsigned int sign : signLength; //знак           %2

    //нормал: 2^(exp-15) * 1.man, exp!=0
    //сабнормал: 2^-14 * 0.man, exp=0
    
public:
    FP16(int _sign = 0, int _exp = 0, int _man = 0) :sign(_sign), exp(_exp), man(_man) {}
    
    bool IsSubnormal() const{
        if (exp==0) return true;
        else return false;
    }
    
    bool IsInf() const{
        if (exp==((1<<expLength)-1)) return true;
        else return false;
    }
    
    void printFP16_bites(){
        std::cout<<sign<<"_";
        for(char i = 0; i<expLength;++i){
            if ((exp & (1 << (expLength - i - 1))) == 0) std::cout<<0;
            else std::cout<<1;
        }
        std::cout<<"_";
        for(char i = 0; i<manLength;++i){
            if ((man & (1 << (manLength - i - 1))) == 0) std::cout<<0;
            else std::cout<<1;
        }
    }

    
    
};

void print_float16_bites(float16_t f){
    uint8_t bytes[2];
    bool bites[16];
    
    for (int i = 0; i < sizeof(f); i++) bytes[i] = *((uint8_t*)(&f) + i);
    
    for (int i = 0; i < 16; ++i) bites[15 - i] = (bytes[i / 8] >> (i % 8)) & 1;
    
    for (int i = 0; i < 16; ++i) {
        if (i==1 || i==6) std::cout<<"_";
        std::cout << bites[i];
    }
}


unsigned short float16_bites_to_short(float16_t f){
    unsigned short temp = 0;
    uint8_t bytes[2];
    bool bites[16];
    
    for (int i = 0; i < sizeof(f); i++) bytes[i] = *((uint8_t*)(&f) + i);
    
    for (int i = 0; i < 16; ++i) bites[15 - i] = (bytes[i / 8] >> (i % 8)) & 1;

    for (int i = 0; i < 16; ++i) temp |= (bites[15-i]<<i);
    
    return temp;
}






using namespace std;
int main(){
    float A = 1.1754944e-38;
    float B = 0;
    float16_t a = A;
    float16_t b = B;

    

    cout<<endl;
    return 0;
}



