#ifndef qtCWENO_CellData_H //include guard
#define qtCWENO_CellData_H

#include <array>

class CellData
{
    //TODO: разобраться, почему без инициализации (Qn = {}) получаются неверные расчеты (возможно, изменилось при переходе на std::array)
    std::array<std::array<double, EQ_NUM>, RK_ORDER_MAX + 1> Qn = {}; //хранимые консервативные величины на подшагах Рунге-Кутты
    double y = 1.0; //для осесимметричных координат, в декартовых задается равным 1.0

public:
    //вычисления по консервативным величинам
    double rho(rkStep rk = 0) const;
    double u(rkStep rk = 0) const;
    double v(rkStep rk = 0) const;
    double p(rkStep rk = 0) const;
    double Mach(rkStep rk = 0) const;
    double lambdaGLFX(rkStep rk = 0) const; //|u|+a
    double lambdaGLFY(rkStep rk = 0) const; //|v|+a
    double lambdaX(Equation eq, rkStep rk = 0) const; //depends on equation
    double lambdaY(Equation eq, rkStep rk = 0) const;
    double F(Equation eq, rkStep rk = 0) const;
    double G(Equation eq, rkStep rk = 0) const;

    void clear(); //обнуление данных
    void set(const CellData& d); //запись данных
    void setY(double _y); //задание y
    void add(const CellData& d); //добавление к текущим консервативным величинам
    void divide(double divisor); //деление консервативных величин
    //void flipVelocity(Neighbour n); //изменение знака компоненты скорости
    std::string dumpQn() const; //дамп всех Qn

    friend std::ostream& operator<<(std::ostream& os, CellData& d); //output overload
};


#endif
