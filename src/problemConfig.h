#ifndef qtCWENO_ProblemConfig_H //include guard
#define qtCWENO_ProblemConfig_H

#include "main.h"

//compile-time константы
const double DOUBLE_EPS12 = 1e-12; //эпсилон для проверки равенства чисел типа double
const int QUADRATURE_POINTS_NUM = 2; //число точек Гауссовой квадратуры на ребре
const int POLY_COEFF_NUM = 5; //число коэффициентов параболоида CWENO: px, py, pxx, pxy, pyy
const int TECPLOT_FIELDS_NUMBER = 2 + EQ_NUM + 3; //x,y + rho,u,v,p + Mach,level,gradrho
const int RK_ORDER_MAX = 3; //максимальный порядок метода Рунге-Кутты

//константы для шаблонов циклов по деревьям и нодам
const bool INCLUDE_GHOSTS = true;
const bool SKIP_GHOSTS = false;
const bool INCLUDE_BRANCHES = true;
const bool SKIP_BRANCHES = false;


#include "CellBox.h"
class ProblemConfig //runtime константы
{
public:
    size_t Nx = 0, Ny = 0; //число деревьев по осям
    int max_depth = 0; //макс глубина деревьев
    CellBox global_box; //глобальная область расчета

    FlowType problem = FlowType::bubble; //тип течения
    double shock_position_x = 0.0; //положение исходного скачка по оси x
    double layer_right = 0.0; //границы слоя
    double layer_bottom = -999.0;
    double layer_top = 0.0;
    double bubble_axle_x = 0.2; //полуоси пузыря
    double bubble_axle_y = 0.2;
    double wedge_angle_bottom = 0.0; //грани клина
    double wedge_angle_top = 90.0;
    double gamma = 1.4; //показатель адиабаты
    double Mach = 1.0; //число Маха исходного скачка
    double Atwood = 0.0; //число Атвуда газа внутри слоя/пузыря
    CoordType coord_type = CoordType::cartesian; //тип координат
    BCType boundary_conditions[DIRECTIONS_NUM] = { BCType::ddn_zero, BCType::ddn_zero, BCType::ddn_zero, BCType::ddn_zero }; //граничные условия
    double end_time = 1.0; //время окончания расчета

    rkStep rk_order = 1; //порядок метода Рунге-Кутты
    double Courant = 0.5; //число Куранта
    FluxType flux = FluxType::LF; //метод расщепления потоков
    size_t meshRefinePeriod = 1; //как часто обновлять сетку

    double meshRefinelevel0 = 0.001; //для mesh refine: стартовый уровень и множитель для последующих уровней
    double meshRefinefactor = 3;
    double meshInitialRefinePadding = 0.01; //ширина полосы дополнительного измельчения вокруг линий начального измельчения
    bool meshRefineAll = false; //true = измельчить всё до упора и далее не склеивать

    double exportPeriod = 0.1; //как часто выводить в файл
    bool exportScatter = false; //выводить ли scatter (ASCII)
    bool exportNeighbours = false; //выводить ли соседей (ASCII)
    bool exportNeighbours12 = false; //выводить ли соседей с учетом диагональных и разделения ребер (ASCII)
    bool exportEdges = false; //выводить ли ребра (ASCII)
    bool exportNodeEdges = false; //выводить ли ребра ячеек (ASCII)
    bool exportEdgeFluxes = false; //выводить ли потоки (ASCII)
    bool exportPlt = true; //выводить ли бинарные .plt-файлы

    ProblemConfig() {}; //дефолтный конструктор
    ProblemConfig(std::string filename); //конструктор по json-файлу
    friend std::ostream& operator<<(std::ostream& os, const ProblemConfig& c); //output overload
};

#endif