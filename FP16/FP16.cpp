// FP16.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>

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
};


int main(){
    
    return 0;
}
