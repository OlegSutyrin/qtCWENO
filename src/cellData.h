#ifndef qtCWENO_cellData_H //include guard
#define qtCWENO_cellData_H

class cellData
{
    //TODO: разобраться, почему без инициализации (Qn = {}) получаются неверные расчеты
    double Qn[RK_ORDER_MAX + 1][EQ_NUM] = {}; //хранимые консервативные величины на подшагах Рунге-Кутты
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
    void setY(double _y); //задание y
    //void flipVelocity(Neighbour n); //изменение знака компоненты скорости
    std::string dumpQn() const; //дамп всех Qn

    friend std::ostream& operator<<(std::ostream& os, cellData& d); //output overload
};


#endif
