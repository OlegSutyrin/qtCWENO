#ifndef qtCWENO_QuadTreeForest_H //include guard
#define qtCWENO_QuadTreeForest_H

#include <vector>
#include <fstream>

#include "QuadTree.h"
#include "NodeEdge.h"

class QuadTreeForest {
    std::vector<QuadTree> trees; //деревья
    std::vector<std::vector<NodeTag>> toRefine; //поуровневый список нод, подлежащих балансировке
    std::vector<NodeEdge> edges; //стороны ячеек
    std::vector<nodeEdgeId> vacant_edge_ids; //список свободных мест в edges

public:
    //accessors
    QuadTree& treeRef(quadTreeId id); //ссылка на дерево по id
    const QuadTree& treeRefConst(quadTreeId id) const; //const версия
    QuadTree& getTreeByCoords(Point p); //поиск дерева по координатам точки
    NodeEdge& edgeRef(nodeEdgeId eid); //ссылка на ребро по id

    //mutators
    void initialize(); //выделение памяти под деревья и прочее
    void addTree(QuadTree tree); //добавление дерева в список
    void addNodeToRefine(const NodeTag& t); //добавление ноды в список на дробление
    nodeEdgeId addEdge(NodeEdge edge); //внесение ребра в список
    nodeEdgeId addEdgeUnique(NodeEdge edge); //внесение ребра в список с проверкой уникальности
    void updateEdge(nodeEdgeId eid, NodeTag n1, NodeTag n2); //обновление данных ребра
    void removeEdge(nodeEdgeId eid); //удаление ребра из списка

    //inspectors
    size_t activeNodesNumber() const; //подсчет активных (неудаленных) нод
    size_t leavesNumber() const; //подсчет листьев

    //other
    dataExtrema getExtrema(); //сбор экстремумов всех величин для вывода в Tecplot
    nodeEdgeId getVacantEdgeId(); //получение номера вакантной ячейки или содание новой в векторе edges
    void meshApplyRefineList(); //дробление ячеек из списка toRefine и балансировка
    void meshRefineInitial(); //начальное дробление сетки
    void meshCoarsenInitial(); //начальное склеивание сетки (для теста)
    void meshUpdate(); //обновление сетки
    void computeQuadraturePoints(); //расчет точек квадратуры во всех ребрах
    void updateEigenObjects(); //создание или обновление Eigen объектов для всех ячеек
    void initialCondition(); //начальные условия
    double CFLTimestepSize(); //величина шага по времени
    void putQn(rkStep rk); //переброс Qn[rk_order] -> Qn[0]
    void computePolynomialCoeffs(rkStep rk); //вычисление коэффициентов CWENO полинома во всех ячейках
    void computeFluxesCWENO(rkStep rk); //расчет потоков на всех ребрах
    void advanceTime(); //полный шаг по времени: вычисление Qn[rk_order]
    void boundaryConditions(rkStep rk); //граничные условия
    void boundaryConditionsAll(); //граничные условия для всех rk

    //output
    void exportScatter(std::string filename); //вывод в файл (Tecplot ASCII scatter)
    void exportNeighbours(std::string filename); //вывод соседей в файл (Tecplot ASCII vectors)
    void exportNeighbours12(std::string filename); //вывод соседей12 в файл (Tecplot ASCII vectors)
    void exportEdges(std::string filename); //вывод ребер в файл (Tecplot ASCII scatter and vectors)
    void exportNodeEdges(std::string filename); //вывод ребер ячеек в файл (Tecplot ASCII vectors)

    friend class QuadTree;
    friend class TreeNode;

    //шаблон цикла по деревьям
    /*template<typename Func>
    inline void forAllTrees(Func f, bool include_ghosts = SKIP_GHOSTS)
    {
        for (auto& tree : trees) //по ссылке
        {
            if (!include_ghosts && tree.isGhost()) //пропуск ghost-деревьев
                continue;
            f(tree);
        }
    }

    //шаблон цикла по нодам
    template<typename Func>
    inline void forAllNodes(Func f, bool include_ghosts = SKIP_GHOSTS, bool include_branches = SKIP_BRANCHES)
    {
        for (auto& tree : trees) //по ссылке
        {
            if (!include_ghosts && tree.isGhost()) //пропуск ghost-деревьев
                continue;
            tree.forAllNodes(f, include_branches);
        }
    }*/
};

//коэффициенты Рунге-Кутты
const double RKcoeff[RK_ORDER_MAX][RK_ORDER_MAX][RK_ORDER_MAX + 1] = { //[order-1][step][Qn], [..][..][RK_ORDER_MAX] - производные
    {
        {1.0, 0.0, 0.0, 1.0}, //1 порядок: explicit Euler
        {},
        {}
    },
    {
        {1.0, 0.0, 0.0, 1.0}, //2 Heun's method
        {0.5, 0.5, 0.0, 0.5}, //(Butcher tableau: a21 = 1 | b1 = 1/2, b2 = 1/2)
        {}
    },
    {
        {1.0, 0.0, 0.0, 1.0}, //3 Strong stability preserving (TVD) Runge-Kutta (SSPRK3)
        {0.75, 0.25, 0.0, 0.25},
        {1.0 / 3.0, 0.0, 2.0 / 3.0, 2.0 / 3.0} //(Butcher tableau: a21 = 1 | a31 = 1/4, a32 = 1/4 | b1 = 1/6, b2 = 1/6, b3 = 2/3)
    }
};


#endif