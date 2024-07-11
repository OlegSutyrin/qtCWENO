#include "main.h"
#include "cellData.h"



//вычисления по консервативным величинам
double cellData::rho(rkStep rk) const { return Qn[rk][static_cast<int>(Equation::density)] / y; }; //деление на y для осесимметричных координат (в декартовых y в cellData задается равным 1)
double cellData::u(rkStep rk) const { return Qn[rk][static_cast<int>(Equation::momentum_x)] / y / Qn[rk][static_cast<int>(Equation::density)]; };
double cellData::v(rkStep rk) const { return Qn[rk][static_cast<int>(Equation::momentum_y)] / y / Qn[rk][static_cast<int>(Equation::density)]; };
double cellData::p(rkStep rk) const
{
    return (config.gamma - 1.0) * (Qn[rk][static_cast<int>(Equation::energy)] / y - 0.5 * (Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::momentum_x)]
        + Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::momentum_y)]) / Qn[rk][static_cast<int>(Equation::density)]);
}
double cellData::Mach(rkStep rk) const { return sqrt((u(rk) * u(rk) + v(rk) * v(rk)) / config.gamma / p(rk) * rho(rk)); };
double cellData::lambdaGLFX(rkStep rk) const { return fabs(u(rk)) + sqrt(config.gamma * p(rk) / rho(rk)); }; //|u|+a
double cellData::lambdaGLFY(rkStep rk) const { return fabs(v(rk)) + sqrt(config.gamma * p(rk) / rho(rk)); }; //|v|+a
double cellData::lambdaX(Equation eq, rkStep rk) const //depends on equation
{
    switch (eq)
    {
    case Equation::density:
        return u(rk) - sqrt(config.gamma * p(rk) / rho(rk));
        break;
    case Equation::momentum_y:
        return u(rk) + sqrt(config.gamma * p(rk) / rho(rk));
        break;
    default:
        return u(rk);
        break;
    }
}
double cellData::lambdaY(Equation eq, rkStep rk) const
{
    switch (eq)
    {
    case Equation::density:
        return v(rk) - sqrt(config.gamma * p(rk) / rho(rk));
        break;
    case Equation::momentum_y:
        return v(rk) + sqrt(config.gamma * p(rk) / rho(rk));
        break;
    default:
        return v(rk);
        break;
    }
}
double cellData::F(Equation eq, rkStep rk) const
{
    switch (eq)
    {
    case Equation::density:
        return Qn[rk][static_cast<int>(Equation::momentum_x)] * y; //rho u (*r)
        break;
    case Equation::momentum_x:
        return (0.5 * (3.0 - config.gamma) * Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::momentum_x)] / Qn[rk][static_cast<int>(Equation::density)]
            + 0.5 * (1.0 - config.gamma) * Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::momentum_y)] / Qn[rk][static_cast<int>(Equation::density)]
            + (config.gamma - 1.0) * Qn[rk][static_cast<int>(Equation::energy)]) * y; //rho*u*u + p (*r)
    case Equation::momentum_y:
        return Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::momentum_y)] / Qn[rk][static_cast<int>(Equation::density)] * y; //rho*u*v (*r)
        break;
    case Equation::energy:
        return (config.gamma * Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::energy)] / Qn[rk][static_cast<int>(Equation::density)]
            - 0.5 * (config.gamma - 1.0) * Qn[rk][static_cast<int>(Equation::momentum_x)] * (Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::momentum_x)]
                + Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::momentum_y)]) / Qn[rk][static_cast<int>(Equation::density)] / Qn[rk][static_cast<int>(Equation::density)]) * y; //(rhoe+p)*u (*r)
        break;
    default: //to suppress warning
        return 0;
    }
}
double cellData::G(Equation eq, rkStep rk) const
{
    switch (eq)
    {
    case Equation::density:
        return Qn[rk][static_cast<int>(Equation::momentum_y)] * y; //rho v (*r)
        break;
    case Equation::momentum_x:
        return Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::momentum_x)] / Qn[rk][static_cast<int>(Equation::density)] * y; //rho*v*u (*r)
    case Equation::momentum_y:
        return (0.5 * (1.0 - config.gamma) * Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::momentum_x)] / Qn[rk][static_cast<int>(Equation::density)]
            + 0.5 * (3.0 - config.gamma) * Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::momentum_y)] / Qn[rk][static_cast<int>(Equation::density)]
            + (config.gamma - 1.0) * Qn[rk][static_cast<int>(Equation::energy)]) * y; //rho*v*v + p (*r)
        break;
    case Equation::energy:
        return (config.gamma * Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::energy)] / Qn[rk][static_cast<int>(Equation::density)]
            - 0.5 * (config.gamma - 1.0) * Qn[rk][static_cast<int>(Equation::momentum_y)] * (Qn[rk][static_cast<int>(Equation::momentum_x)] * Qn[rk][static_cast<int>(Equation::momentum_x)]
                + Qn[rk][static_cast<int>(Equation::momentum_y)] * Qn[rk][static_cast<int>(Equation::momentum_y)]) / Qn[rk][static_cast<int>(Equation::density)] / Qn[rk][static_cast<int>(Equation::density)]) * y; //(rhoe+p)*v (*r)
        break;
    default:
        return 0;
    }
}

void cellData::clear() //обнуление данных
{
    for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
        for (auto& eq : Equations)
        {
            Qn[rk][static_cast<int>(eq)] = 0;
        }
    //y = 1.0; не нужно?
}

void cellData::setY(double _y) { y = _y; } //задание y

/*void cellData::flipVelocity(Neighbour n)  //изменение знака компоненты скорости
{
    if (n == NEIGHBOUR_LEFT || n == NEIGHBOUR_RIGHT) //горизонтальный сосед
    {
        for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
            Qn[rk][Equation::momentum_x] *= -1;
    }
    else //вертикальный сосед
    {
        for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
            Qn[rk][Equation::momentum_y] *= -1;
    }
}*/

std::string cellData::dumpQn() const //дамп всех Qn
{
    std::stringstream buffer;
    buffer << "( ";
    for (rkStep rk = 0; rk <= config.rk_order; rk++)
    {
        for (auto eq : Equations)
            buffer << Qn[rk][static_cast<int>(eq)] << " ";
        buffer << "| ";
    }
    buffer << ")" << endl;
    return buffer.str();

}

std::ostream& operator<<(std::ostream& os, cellData& d) //output operator overload
{
    os << "(" << d.rho() << ", " << d.u() << ", " << d.v() << ", " << d.p() << ")";
    return os;
}