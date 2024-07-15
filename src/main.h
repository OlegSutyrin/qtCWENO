#ifndef qtCWENO_main_H //include guard
#define qtCWENO_main_H

#include <iostream>
#include <vector>
#include <fstream>

#pragma warning(error:4706) //������ ���������� ������ if() ("warning C4706: assignment within conditional expression")
using std::cout;
using std::endl;

using quadTreeId = size_t; //TODO: �������� �� "enum class quadTreeId : size_t {}", �� ����� ������ �������� ��������� static_cast'��
using treeNodeId = size_t; //(� �������� ����� ����� overload []?)
using cellDataId = size_t;
using nodeEdgeId = size_t;
const size_t null = size_t(-1); //���������� �������� ��� ��������� �����
using rkStep = size_t;

enum class Directions { //������������ �����������
    up, right, down, left, count //����� �����
};
enum class Orientation { //���������� �����, ����� ����� � �.�.
    horizontal, vertical
};
enum class Quadrant { //��������� (���� ����)
    top_left, top_right, bottom_right, bottom_left, count
};
const Quadrant Quadrants[] = { Quadrant::top_left, Quadrant::top_right, Quadrant::bottom_right, Quadrant::bottom_left }; //��� ������ �� ����������
enum class Neighbour { //������ (������������ ��� ������ �� ���������� ������)
    top, right, bottom, left, count
};
const Neighbour Neighbours[] = { Neighbour::top, Neighbour::right, Neighbour::bottom, Neighbour::left }; //��� ������ �� �������
Neighbour opposite(Neighbour n); //��������������� �����

enum class Neighbour12 { //������ � ������ ���������� � ��������� ����� (������������ ��� ����� � CWENO)
    top1, top2, top_right, right1, right2, bottom_right, bottom1, bottom2, bottom_left, left1, left2, top_left, count
};
const Neighbour12 Neighbours12[] = { Neighbour12::top1, Neighbour12::top2, Neighbour12::top_right, Neighbour12::right1, Neighbour12::right2, Neighbour12::bottom_right, Neighbour12::bottom1, Neighbour12::bottom2, Neighbour12::bottom_left, Neighbour12::left1, Neighbour12::left2, Neighbour12::top_left }; //��� ������
Neighbour12 opposite(Neighbour12 n12); //��������������� �����
//Neighbour12 next(Neighbour12 n12);


enum class Equation : unsigned int { //���������
    density, momentum_x, momentum_y, energy, count
};
const Equation Equations[] = { Equation::density, Equation::momentum_x, Equation::momentum_y, Equation::energy }; //��� ������ �� ����������

enum class CoordType { //���� ���������
    cartesian, //���������
    axisymmetric //��������������� ��������������
};
enum class BCType { //���� ��������� �������
    ddn_zero, //"������������" 1�� �������
    wall //������ (=������� ���������)
};
enum class FluxType { //������ ���������� ������� �� ������
    LF, //Lax-Friedrich
    Riemann_exact, //exact Riemann solver
    HLLC //HLLC Riemann solver
};

#include "output.h"
class Globals
{
public:
    size_t timestep_number = 0;
    size_t export_number = 0;
    double time = 0.0;
    double dt = 0.001;
    double last_exported_time = 0.0;
    std::ofstream file_output;
    std::vector<double> refineLvls;
};

//���������� �������
#include "ProblemConfig.h"
extern ProblemConfig config; //������ ������
extern Globals globals; //���������� ����������

#include "QuadTreeForest.h"
extern QuadTreeForest forest; //��� ��������




#endif