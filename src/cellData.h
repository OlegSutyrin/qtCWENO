#ifndef qtCWENO_cellData_H //include guard
#define qtCWENO_cellData_H

#include <array>

class cellData
{
    //TODO: �����������, ������ ��� ������������� (Qn = {}) ���������� �������� �������
    std::array<std::array<double, EQ_NUM>, RK_ORDER_MAX + 1> Qn = {}; //�������� �������������� �������� �� �������� �����-�����
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
    void set(cellData d); //������ ������
    void setY(double _y); //������� y
    void add(cellData d); //���������� � ������� �������������� ���������
    void divide(double divisor); //������� �������������� �������
    //void flipVelocity(Neighbour n); //��������� ����� ���������� ��������
    std::string dumpQn() const; //���� ���� Qn

    friend std::ostream& operator<<(std::ostream& os, cellData& d); //output overload
};


#endif
