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
using rkStep = std::size_t;

enum class Directions { //кардинальные направления
    up,
    right,
    down,
    left,
    count //общее число
};
enum class Orientation { //ориентации линий, ребер ячеек и т.д.
    horizontal,
    vertical
};
enum class Quadrant { //квадранты (дети ноды)
    top_left,
    top_right,
    bottom_right,
    bottom_left,
    count
};
const Quadrant Quadrants[] = { Quadrant::top_left, Quadrant::top_right, Quadrant::bottom_right, Quadrant::bottom_left }; //для циклов по квадрантам

enum class Equation : unsigned int { //уравнения
    density,
    momentum_x,
    momentum_y,
    energy,
    count
};
const Equation Equations[] = { Equation::density, Equation::momentum_x, Equation::momentum_y, Equation::energy }; //для циклов по уравнениям

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

//глобальные объекты
#include "problemConfig.h"
extern problemConfig config; //конфиг задачи
extern Globals globals; //глобальные переменные
//extern quadTreeForest forest; //лес деревьев




#endif