#ifndef qtCWENO_CellData_H //include guard
#define qtCWENO_CellData_H

#include <array>

//struct вектор консервативных величин
struct ConservativeVector
{
    std::array<double, EQ_NUM> Q = {};
    //вычисления по консервативным величинам
    double rho(double y) const;
    double u(double y) const;
    double v(double y) const;
    double p(double y) const;
    double Mach(double y) const;
    double lambdaGLFX(double y) const; //|u|+a
    double lambdaGLFY(double y) const; //|v|+a
    double lambdaX(Equation eq, double y) const; //depends on equation
    double lambdaY(Equation eq, double y) const;
    double F(Equation eq, double y) const;
    double G(Equation eq, double y) const;

    void clear(); //обнуление вектора
    void set(Equation eq, double value); //задание компоненты

    const double operator()(Equation eq) const; //() operator: компонента по уравнению
    ConservativeVector& operator+=(const ConservativeVector& rhs); //+= overload
    ConservativeVector& operator/=(const double rhs); ///= overload
    friend std::ostream& operator<<(std::ostream& os, const ConservativeVector& rhs); //output overload
};

//данные ячейки на подшагах Рунге-Кутты + координата y для осесимметричных расчетов
class CellData
{
    //TODO: разобраться, почему без инициализации (Qn = {}) получаются неверные расчеты (возможно, изменилось при переходе на std::array)
    std::array<ConservativeVector, RK_ORDER_MAX + 1> Qn = {}; //хранимые консервативные величины на подшагах Рунге-Кутты
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
    void setQ0(const ConservativeVector Q); //запись вектора консервативный переменных только для rk=0
    void setY(double _y); //задание y
    void add(const CellData& d); //добавление к текущим консервативным величинам
    void divide(double divisor); //деление консервативных величин
    //void flipVelocity(Neighbour n); //изменение знака компоненты скорости
    std::string dumpQn() const; //дамп всех Qn

    friend std::ostream& operator<<(std::ostream& os, const CellData& d); //output overload
};


#endif
