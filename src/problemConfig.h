#ifndef qtCWENO_ProblemConfig_H //include guard
#define qtCWENO_ProblemConfig_H

#include "main.h"

//compile-time ���������
const double DOUBLE_EPS12 = 1e-12; //������� ��� �������� ��������� ����� ���� double
const int QUADRATURE_POINTS_NUM = 2; //����� ����� ��������� ���������� �� �����
const int POLY_COEFF_NUM = 5; //����� ������������� ����������� CWENO: px, py, pxx, pxy, pyy
const int TECPLOT_FIELDS_NUMBER = 2 + EQ_NUM + 3; //x,y + rho,u,v,p + Mach,level,gradrho
const int RK_ORDER_MAX = 3; //������������ ������� ������ �����-�����

//��������� ��� �������� ������ �� �������� � �����
const bool INCLUDE_GHOSTS = true;
const bool SKIP_GHOSTS = false;
const bool INCLUDE_BRANCHES = true;
const bool SKIP_BRANCHES = false;


#include "CellBox.h"
class ProblemConfig //runtime ���������
{
public:
    size_t Nx = 0, Ny = 0; //����� �������� �� ����
    int max_depth = 0; //���� ������� ��������
    CellBox global_box; //���������� ������� �������

    FlowType problem = FlowType::bubble; //��� �������
    double shock_position_x = 0.0; //��������� ��������� ������ �� ��� x
    double layer_right = 0.0; //������� ����
    double layer_bottom = -999.0;
    double layer_top = 0.0;
    double bubble_axle_x = 0.2; //������� ������
    double bubble_axle_y = 0.2;
    double wedge_angle_bottom = 0.0; //����� �����
    double wedge_angle_top = 90.0;
    double gamma = 1.4; //���������� ��������
    double Mach = 1.0; //����� ���� ��������� ������
    double Atwood = 0.0; //����� ������ ���� ������ ����/������
    CoordType coord_type = CoordType::cartesian; //��� ���������
    BCType boundary_conditions[DIRECTIONS_NUM] = { BCType::ddn_zero, BCType::ddn_zero, BCType::ddn_zero, BCType::ddn_zero }; //��������� �������
    double end_time = 1.0; //����� ��������� �������

    rkStep rk_order = 1; //������� ������ �����-�����
    double Courant = 0.5; //����� �������
    FluxType flux = FluxType::LF; //����� ����������� �������
    size_t meshRefinePeriod = 1; //��� ����� ��������� �����

    double meshRefinelevel0 = 0.001; //��� mesh refine: ��������� ������� � ��������� ��� ����������� �������
    double meshRefinefactor = 3;
    double meshInitialRefinePadding = 0.01; //������ ������ ��������������� ����������� ������ ����� ���������� �����������
    bool meshRefineAll = false; //true = ���������� �� �� ����� � ����� �� ���������

    double exportPeriod = 0.1; //��� ����� �������� � ����
    bool exportScatter = false; //�������� �� scatter (ASCII)
    bool exportNeighbours = false; //�������� �� ������� (ASCII)
    bool exportNeighbours12 = false; //�������� �� ������� � ������ ������������ � ���������� ����� (ASCII)
    bool exportEdges = false; //�������� �� ����� (ASCII)
    bool exportNodeEdges = false; //�������� �� ����� ����� (ASCII)
    bool exportEdgeFluxes = false; //�������� �� ������ (ASCII)
    bool exportPlt = true; //�������� �� �������� .plt-�����

    ProblemConfig() {}; //��������� �����������
    ProblemConfig(std::string filename); //����������� �� json-�����
    friend std::ostream& operator<<(std::ostream& os, const ProblemConfig& c); //output overload
};

#endif