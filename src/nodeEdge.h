#ifndef qtCWENO_nodeEdge_H //include guard
#define qtCWENO_nodeEdge_H

#include "main.h"
#include "NodeTag.h"

//����� - ������� (��� � ��������) ������, ��� �������� �������
class NodeEdge {
    bool is_deleted = false;
    Orientation orientation_ = Orientation::horizontal;
    NodeTag n1_{}, n2_{}; //���� �����-�������, ����������� "�������, ������" ��� "�����, ������"
    std::array<double, EQ_NUM> FQ_ = {}; //����� ����� �������
    std::array<Point, QUADRATURE_POINTS_NUM> qps = {}; //����� ��� ����������, ����������� "������, �������" ��� "�����, ������"

public:
    NodeEdge() {}; //default contructor
    NodeEdge(NodeTag n1, NodeTag n2, Orientation orientation) : n1_(n1), n2_(n2), orientation_(orientation) {}; //����������� �� ���� �������
    
    //accessors
    const NodeTag n1() const;
    const NodeTag n2() const;
    double FQ(Equation eq) const; //���������� ������
    Point middle(); //����� �����
    double length() const; //���������� ����� ����� �� ������ ���������� TODO: ������� ������� � ����� � ����� �����

    //mutators
    void markDeleted();
    void setN1(NodeTag ntag);
    void setN2(NodeTag ntag);
    
    //inspectors
    bool isDeleted();

    //other
    void computeQuadraturePoints(); //������ ����� ����������
    void computeFluxLF(rkStep rk); //������ ������ (Lax-Friedrich flux)
    void computeFluxHLLC(rkStep rk); //������ ������ (HLLC Rieman solver)
    void computeFluxRiemannExact(rkStep rk); //������ ������ (exact Rieman solver)

    bool operator==(const NodeEdge& rhs); //equality overload

    //output
    std::string dumpNeigboursVectors(); //���� ������� � ���� ��������
    std::string dumpQuadraturePoints(); //���� ����� ����������
    std::string dumpFlux(); //���� �������� �������

    friend std::ostream& operator<<(std::ostream& os, const NodeEdge& edge); //output overload
};

#endif
