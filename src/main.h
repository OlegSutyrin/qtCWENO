#ifndef qtCWENO_main_H //include guard
#define qtCWENO_main_H

#include <iostream>
#include <array>
#include <vector>
#include <fstream>

#pragma warning(error:4706) //запрет присовения внутри if() ("warning C4706: assignment within conditional expression")
using std::cout;
using std::endl;

using quadTreeId = size_t; //TODO: заменить на "enum class quadTreeId : size_t {}", но найти способ избежать множества static_cast'ов
using treeNodeId = size_t; //(в массивах можно через overload []?)
using cellDataId = size_t;
using nodeEdgeId = size_t;
const size_t null = size_t(-1); //наибольшее значение для служебных целей
using rkStep = size_t;

enum class Directions { //кардинальные направления
    up, right, down, left, count //общее число
};
const int DIRECTIONS_NUM = static_cast<int>(Directions::count); //число кардинальных направлений системы координат

enum class Orientation { //ориентации линий, ребер ячеек и т.д.
    horizontal, vertical
};

enum class Quadrant { //квадранты (дети ноды)
    top_left, top_right, bottom_right, bottom_left, count
};
const int QUADRANTS_NUM = static_cast<int>(Quadrant::count); //число вершин ячейки
const std::array<Quadrant, QUADRANTS_NUM> Quadrants = { Quadrant::top_left, Quadrant::top_right, Quadrant::bottom_right, Quadrant::bottom_left }; //для циклов по квадрантам

enum class Neighbour { //соседи (используются в основном для структуры дерева)
    top, right, bottom, left, count
};
const int NEIGHBOURS_NUM = static_cast<int>(Neighbour::count); //число соседей
const std::array <Neighbour, NEIGHBOURS_NUM> Neighbours = { Neighbour::top, Neighbour::right, Neighbour::bottom, Neighbour::left }; //для циклов по соседям

enum class Neighbour12 { //соседи с учетом диагоналей и разбиения ребер (используются для CWENO)
    top1, top2, top_right, right1, right2, bottom_right, bottom1, bottom2, bottom_left, left1, left2, top_left, count
};
const int MAX_NEIGHBOURS12_NUM = static_cast<int>(Neighbour12::count); //число всех возможных соседей с учетом диагональных
const std::array <Neighbour12, MAX_NEIGHBOURS12_NUM> Neighbours12 = { Neighbour12::top1, Neighbour12::top2, Neighbour12::top_right, Neighbour12::right1, Neighbour12::right2, Neighbour12::bottom_right, Neighbour12::bottom1, Neighbour12::bottom2, Neighbour12::bottom_left, Neighbour12::left1, Neighbour12::left2, Neighbour12::top_left }; //для циклов
const std::array <Neighbour12, MAX_NEIGHBOURS12_NUM - 4> Neighbours12Cardinal = { Neighbour12::top1, Neighbour12::top2, Neighbour12::right1, Neighbour12::right2, Neighbour12::bottom1, Neighbour12::bottom2, Neighbour12::left1, Neighbour12::left2 }; //для циклов по соседям по стороне
const std::array <Neighbour12, 4> Neighbours12Diagonal = { Neighbour12::top_right, Neighbour12::bottom_right, Neighbour12::bottom_left, Neighbour12::top_left }; //для циклов по диагональным соседям

enum class Edge { //ребра ячейки, до двух на каждую сторону
    top1, top2, right1, right2, bottom1, bottom2, left1, left2, count
};
const int EDGES_NUM = static_cast<int>(Edge::count); //число возможных ребер в ячейке
const std::array <Edge, EDGES_NUM> Edges = { Edge::top1, Edge::top2, Edge::right1, Edge::right2, Edge::bottom1, Edge::bottom2, Edge::left1, Edge::left2 }; //для циклов

enum class Equation : unsigned int { //уравнения
    density, momentum_x, momentum_y, energy, count
};
const int EQ_NUM = static_cast<int>(Equation::count); //число уравнений
const std::array <Equation, EQ_NUM> Equations = { Equation::density, Equation::momentum_x, Equation::momentum_y, Equation::energy }; //для циклов по уравнениям

enum class CoordType { //типы координат
    cartesian, //декартовы
    axisymmetric //осесимметричные цилиндрические
};
enum class BCType { //типы граничных условий
    ddn_zero, //"неотражающее" 1го порядка
    wall //стенка (=плоская симметрия)
};
enum class FluxType { //методы вычисления потоков на ребрах
    LF, //Lax-Friedrich
    Riemann_exact, //exact Riemann solver
    HLLC //HLLC Riemann solver
};

//преобразования enum class'ов
Neighbour opposite(Neighbour n); //противоположный сосед
Neighbour12 opposite(Neighbour12 n12); //противоположный (по ребру или вершине) сосед
Quadrant toQuadrant(Neighbour12 n12); //квадрант по соседу12
Neighbour12 toNeighbour12(Neighbour n); //сосед12 по соседу4
Edge toEdge(Neighbour n); //ребро по соседу4
Edge next(Edge e); //следующее ребро
Edge opposite(Edge e); //противоположное (по стороне) ребро

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

//глобальные объекты
#include "ProblemConfig.h"
extern ProblemConfig config; //конфиг задачи
extern Globals globals; //глобальные переменные

#include "QuadTreeForest.h"
extern QuadTreeForest forest; //лес деревьев




#endif