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
    void initialCondition(); //начальные условия

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

#endif