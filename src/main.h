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
using ushorty = std::uint16_t; //16-������ ��������������� ����� (8-������ �� ��������� ��������� ��� ������, ��� ��������)
using rkStep = std::uint16_t;

enum class Directions { //������������ �����������
    up,
    right,
    down,
    left
};
enum class Orientation { //���������� �����, ����� ����� � �.�.
    horizontal,
    vertical
};
enum class Quadrant { //��������� (���� ����)
    top_left,
    top_right,
    bottom_right,
    bottom_left
};
const Quadrant Quadrants[] = { Quadrant::top_left, Quadrant::top_right, Quadrant::bottom_right, Quadrant::bottom_left }; //��� ������ �� ����������

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

//���������� �������
#include "problemConfig.h"
extern problemConfig config; //������ ������
//extern Globals globals; //���������� ����������
//extern quadTreeForest forest; //��� ��������




#endif