#ifndef qtCWENO_main_H //include guard
#define qtCWENO_main_H

#include <iostream>
#include <vector>
#include <fstream>

#pragma warning(error:4706) //запрет присовения внутри if() ("warning C4706: assignment within conditional expression")
using std::cout;
using std::endl;

using treeId = size_t;
using nodeId = size_t;
using dataId = size_t;
using edgeId = size_t;
const size_t null = size_t(-1); //наибольшее значение для служебных целей
using ushorty = std::uint16_t; //16-битное неотрицательное целое (8-битное по умолчанию выводится как символ, что неудобно)
using rkStep = std::uint16_t;

enum class Directions { //кардинальные направления
    up,
    right,
    down,
    left
};
enum class Orientation { //ориентации линий, ребер ячеек и т.д.
    horizontal,
    vertical
};
enum class Quadrant { //квадранты (дети ноды)
    top_left,
    top_right,
    bottom_right,
    bottom_left
};
const Quadrant Quadrants[] = { Quadrant::top_left, Quadrant::top_right, Quadrant::bottom_right, Quadrant::bottom_left }; //для циклов по квадрантам

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

//глобальные объекты
#include "problemConfig.h"
extern problemConfig config; //конфиг задачи
//extern Globals globals; //глобальные переменные
//extern quadTreeForest forest; //лес деревьев




#endif