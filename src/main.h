#ifndef qtCWENO_main_H //include guard
#define qtCWENO_main_H

#include <iostream>
#include <vector>
#include <fstream>

#pragma warning(error:4706) //������ ���������� ������ if() ("warning C4706: assignment within conditional expression")
using std::cout;
using std::endl;

using treeId = size_t;
using nodeId = size_t;
using dataId = size_t;
using edgeId = size_t;
const size_t null = size_t(-1); //���������� �������� ��� ��������� �����
using rkStep = std::size_t;

enum class Directions { //������������ �����������
    up,
    right,
    down,
    left,
    count //����� �����
};
enum class Orientation { //���������� �����, ����� ����� � �.�.
    horizontal,
    vertical
};
enum class Quadrant { //��������� (���� ����)
    top_left,
    top_right,
    bottom_right,
    bottom_left,
    count
};
const Quadrant Quadrants[] = { Quadrant::top_left, Quadrant::top_right, Quadrant::bottom_right, Quadrant::bottom_left }; //��� ������ �� ����������

enum class Equation : unsigned int { //���������
    density,
    momentum_x,
    momentum_y,
    energy,
    count
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
struct Globals
{
    size_t timestep_number = 0;
    size_t export_number = 0;
    double time = 0.0;
    double dt = 0.001;
    double last_exported_time = 0.0;
    std::ofstream file_output;
    std::vector<double> refineLvls;
};

//���������� �������
#include "problemConfig.h"
extern problemConfig config; //������ ������
extern Globals globals; //���������� ����������
//extern quadTreeForest forest; //��� ��������




#endif