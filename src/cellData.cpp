#include "main.h"
#include "CellData.h"

//------------------------------------------------------------------- ConservativeVector -------------------------------------------------------------------
double ConservativeVector::rho(double y) const { return Q[static_cast<int>(Equation::density)] / y; }; //������� �� y ��� ��������������� ��������� (� ���������� y � ConservativeVector �������� ������ 1)
double ConservativeVector::u(double y) const { return Q[static_cast<int>(Equation::momentum_x)] / y / Q[static_cast<int>(Equation::density)]; };
double ConservativeVector::v(double y) const { return Q[static_cast<int>(Equation::momentum_y)] / y / Q[static_cast<int>(Equation::density)]; };
double ConservativeVector::p(double y) const
{
    return (config.gamma - 1.0) * (Q[static_cast<int>(Equation::energy)] / y - 0.5 * (Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::momentum_x)]
        + Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::momentum_y)]) / Q[static_cast<int>(Equation::density)]);
}
double ConservativeVector::Mach(double y) const { return sqrt((u(y) * u(y) + v(y) * v(y)) / config.gamma / p(y) * rho(y)); };
double ConservativeVector::lambdaGLFX(double y) const { return fabs(u(y)) + sqrt(config.gamma * p(y) / rho(y)); }; //|u|+a
double ConservativeVector::lambdaGLFY(double y) const { return fabs(v(y)) + sqrt(config.gamma * p(y) / rho(y)); }; //|v|+a
double ConservativeVector::lambdaX(Equation eq, double y) const //depends on equation
{
    switch (eq)
    {
    case Equation::density:
        return u(y) - sqrt(config.gamma * p(y) / rho(y));
        break;
    case Equation::momentum_y:
        return u(y) + sqrt(config.gamma * p(y) / rho(y));
        break;
    default:
        return u(y);
        break;
    }
}
double ConservativeVector::lambdaY(Equation eq, double y) const
{
    switch (eq)
    {
    case Equation::density:
        return v(y) - sqrt(config.gamma * p(y) / rho(y));
        break;
    case Equation::momentum_y:
        return v(y) + sqrt(config.gamma * p(y) / rho(y));
        break;
    default:
        return v(y);
        break;
    }
}
double ConservativeVector::F(Equation eq, double y) const
{
    switch (eq)
    {
    case Equation::density:
        return Q[static_cast<int>(Equation::momentum_x)] * y; //rho u (*r)
        break;
    case Equation::momentum_x:
        return (0.5 * (3.0 - config.gamma) * Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::momentum_x)] / Q[static_cast<int>(Equation::density)]
            + 0.5 * (1.0 - config.gamma) * Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::momentum_y)] / Q[static_cast<int>(Equation::density)]
            + (config.gamma - 1.0) * Q[static_cast<int>(Equation::energy)]) * y; //rho*u*u + p (*r)
    case Equation::momentum_y:
        return Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::momentum_y)] / Q[static_cast<int>(Equation::density)] * y; //rho*u*v (*r)
        break;
    case Equation::energy:
        return (config.gamma * Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::energy)] / Q[static_cast<int>(Equation::density)]
            - 0.5 * (config.gamma - 1.0) * Q[static_cast<int>(Equation::momentum_x)] * (Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::momentum_x)]
                + Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::momentum_y)]) / Q[static_cast<int>(Equation::density)] / Q[static_cast<int>(Equation::density)]) * y; //(rhoe+p)*u (*r)
        break;
    default: //to suppress warning
        return 0;
    }
}
double ConservativeVector::G(Equation eq, double y) const
{
    switch (eq)
    {
    case Equation::density:
        return Q[static_cast<int>(Equation::momentum_y)] * y; //rho v (*r)
        break;
    case Equation::momentum_x:
        return Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::momentum_x)] / Q[static_cast<int>(Equation::density)] * y; //rho*v*u (*r)
    case Equation::momentum_y:
        return (0.5 * (1.0 - config.gamma) * Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::momentum_x)] / Q[static_cast<int>(Equation::density)]
            + 0.5 * (3.0 - config.gamma) * Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::momentum_y)] / Q[static_cast<int>(Equation::density)]
            + (config.gamma - 1.0) * Q[static_cast<int>(Equation::energy)]) * y; //rho*v*v + p (*r)
        break;
    case Equation::energy:
        return (config.gamma * Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::energy)] / Q[static_cast<int>(Equation::density)]
            - 0.5 * (config.gamma - 1.0) * Q[static_cast<int>(Equation::momentum_y)] * (Q[static_cast<int>(Equation::momentum_x)] * Q[static_cast<int>(Equation::momentum_x)]
                + Q[static_cast<int>(Equation::momentum_y)] * Q[static_cast<int>(Equation::momentum_y)]) / Q[static_cast<int>(Equation::density)] / Q[static_cast<int>(Equation::density)]) * y; //(rhoe+p)*v (*r)
        break;
    default:
        return 0;
    }
}

void ConservativeVector::clear() { Q.fill(0); } //��������� �������
void ConservativeVector::set(Equation eq, double value) { Q[static_cast<int>(eq)] = value; } //������� ����������
void ConservativeVector::flipVelocity(Orientation ori) //��������� ����� ���������� ��������
{
    if (ori == Orientation::horizontal) //�������� �������������� ����������
        Q[static_cast<int>(Equation::momentum_x)] *= -1;
    else
        Q[static_cast<int>(Equation::momentum_y)] *= -1;
}

const double ConservativeVector::operator()(Equation eq) const { return Q[static_cast<int>(eq)]; } //���������� �� ���������
ConservativeVector& ConservativeVector::operator=(const ConservativeVector& rhs) //= overload
{
    if (this == &rhs)
        return *this;
    Q = rhs.Q; //std::array assignment
    return *this;
}
ConservativeVector& ConservativeVector::operator+=(const ConservativeVector& rhs) //+= overload
{
    for(auto eq : Equations)
        Q[static_cast<int>(eq)] += rhs.Q[static_cast<int>(eq)];
    return *this;
}
ConservativeVector& ConservativeVector::operator/=(const double rhs) ///= overload
{
    for (auto eq : Equations)
        Q[static_cast<int>(eq)] /= rhs;
    return *this;
}
std::ostream& operator<<(std::ostream& os, const ConservativeVector& rhs) //output operator overload
{
    os << "(";
    for (auto eq : Equations)
        os << rhs.Q[static_cast<int>(eq)] << " ";
    os << ")";
    return os;
}

//------------------------------------------------------------------- CellData -------------------------------------------------------------------
ConservativeVector& CellData::Qref(rkStep rk) { return Qn[static_cast<int>(rk)]; } //������ �� ������ ����������
double CellData::rho(rkStep rk) const { return Qn[rk].rho(y); }; //������� �� y ��� ��������������� ��������� (� ���������� y � CellData �������� ������ 1)
double CellData::u(rkStep rk) const { return Qn[rk].u(y); };
double CellData::v(rkStep rk) const { return Qn[rk].v(y); };
double CellData::p(rkStep rk) const { return Qn[rk].p(y); };
double CellData::Mach(rkStep rk) const { return Qn[rk].Mach(y); };
double CellData::lambdaGLFX(rkStep rk) const { return Qn[rk].lambdaGLFX(y); };
double CellData::lambdaGLFY(rkStep rk) const { return Qn[rk].lambdaGLFY(y); };
double CellData::lambdaX(Equation eq, rkStep rk) const { return Qn[rk].lambdaX(eq, y); }
double CellData::lambdaY(Equation eq, rkStep rk) const { return Qn[rk].lambdaY(eq, y); }
double CellData::F(Equation eq, rkStep rk) const { return Qn[rk].F(eq, y); }
double CellData::G(Equation eq, rkStep rk) const { return Qn[rk].G(eq, y); }

void CellData::clear() //��������� ������
{
    for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
        Qn[rk].clear();
    //y = 1.0; �� �����?
}

void CellData::set(const CellData& d) //������ ������
{
    for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
        Qn[rk] = d.Qn[rk]; //std::array assignment
    y = d.y;
}
void CellData::setQ(const ConservativeVector Q, rkStep rk) { Qn[rk] = Q; } //������ ������� �������������� ���������� � ���� [rk]
void CellData::setY(double _y) { y = _y; } //������� y

void CellData::add(const CellData& d) //���������� � ������� �������������� ���������
{
    for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
        Qn[rk] += d.Qn[rk];
}
void CellData::divide(double divisor) //������� �������������� �������
{
    for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
        Qn[rk] /= divisor;
}

/*void CellData::flipVelocity(Neighbour n)  //��������� ����� ���������� ��������
{
    if (n == NEIGHBOUR_LEFT || n == NEIGHBOUR_RIGHT) //�������������� �����
    {
        for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
            Qn[rk][Equation::momentum_x] *= -1;
    }
    else //������������ �����
    {
        for (rkStep rk = 0; rk < RK_ORDER_MAX; rk++)
            Qn[rk][Equation::momentum_y] *= -1;
    }
}*/

void CellData::putQn(rkStep rk) { Qn[0] = Qn[static_cast<int>(rk)]; } //�������� Qn[rk_order] -> Qn[0]

std::string CellData::dumpQn() const //���� ���� Qn
{
    std::stringstream buffer;
    buffer << "( ";
    for (rkStep rk = 0; rk <= config.rk_order; rk++)
    {
        buffer << Qn[rk] << "| ";
    }
    buffer << ")" << endl;
    return buffer.str();

}

std::ostream& operator<<(std::ostream& os, const CellData& d) //output operator overload
{
    os << "(" << d.rho() << ", " << d.u() << ", " << d.v() << ", " << d.p() << ")";
    return os;
}