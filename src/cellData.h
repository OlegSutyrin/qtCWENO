#ifndef qtCWENO_CellData_H //include guard
#define qtCWENO_CellData_H

#include <array>

//struct ������ �������������� �������
struct ConservativeVector
{
    std::array<double, EQ_NUM> Q = {};
    //���������� �� �������������� ���������
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

    void clear(); //��������� �������
    void set(Equation eq, double value); //������� ����������

    const double operator()(Equation eq) const; //() operator: ���������� �� ���������
    ConservativeVector& operator+=(const ConservativeVector& rhs); //+= overload
    ConservativeVector& operator/=(const double rhs); ///= overload
    friend std::ostream& operator<<(std::ostream& os, const ConservativeVector& rhs); //output overload
};

//������ ������ �� �������� �����-����� + ���������� y ��� ��������������� ��������
class CellData
{
    //TODO: �����������, ������ ��� ������������� (Qn = {}) ���������� �������� ������� (��������, ���������� ��� �������� �� std::array)
    std::array<ConservativeVector, RK_ORDER_MAX + 1> Qn = {}; //�������� �������������� �������� �� �������� �����-�����
    double y = 1.0; //��� ��������������� ���������, � ���������� �������� ������ 1.0

public:
    //���������� �� �������������� ���������
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

    void clear(); //��������� ������
    void set(const CellData& d); //������ ������
    void setQ0(const ConservativeVector Q); //������ ������� �������������� ���������� ������ ��� rk=0
    void setY(double _y); //������� y
    void add(const CellData& d); //���������� � ������� �������������� ���������
    void divide(double divisor); //������� �������������� �������
    //void flipVelocity(Neighbour n); //��������� ����� ���������� ��������
    std::string dumpQn() const; //���� ���� Qn

    friend std::ostream& operator<<(std::ostream& os, const CellData& d); //output overload
};


#endif
