#ifndef qtCWENO_QuadTreeForest_H //include guard
#define qtCWENO_QuadTreeForest_H

#include <vector>
#include <fstream>

#include "QuadTree.h"
#include "NodeEdge.h"

class QuadTreeForest {
    std::vector<QuadTree> trees; //деревья
    std::vector<std::vector<NodeTag>> toRefine; //поуровневый список нод, подлежащих балансировке
    //std::vector<nodeEdge> edges; //стороны ячеек
    //std::vector<nodeEdgeId> vacant_edge_ids; //список свободных мест в edges

public:
    //accessors
    QuadTree& treeRef(quadTreeId id); //ссылка на дерево по id
    const QuadTree& treeRefConst(quadTreeId id) const; //const версия
    QuadTree& getTreeByCoords(Point p); //поиск дерева по координатам точки

    //mutators
    void initialize(); //выделение памяти под деревья и прочее
    void addTree(QuadTree tree); //добавление дерева в список
    void addNodeToRefine(const NodeTag& t); //добавление ноды в список на дробление
    //edgeId addEdge(nodeEdge edge); //внесение ребра в список
    //edgeId addEdgeUnique(nodeEdge edge); //внесение ребра в список с проверкой уникальности
    //int updateEdge(edgeId eid, NodeTag n1, NodeTag n2); //обновление данных ребра
    //int removeEdge(edgeId eid); //удаление ребра из списка

    //inspectors
    size_t activeNodesNumber(); //подсчет активных (неудаленных) нод
    size_t leavesNumber(); //подсчет листьев

    //other
    dataExtrema getExtrema(); //сбор экстремумов всех величин для вывода в Tecplot
    //edgeId getVacantEdgeId(); //получение номера вакантной ячейки или содание новой в векторе edges
    void meshRefineInitial(); //начальное дробление сетки
    void meshApplyRefineList(); //дробление ячеек из списка toRefine

    //output
    void exportScatter(std::string filename); //вывод в файл (Tecplot ASCII scatter)
    void exportNeighbours(std::string filename); //вывод в файл (Tecplot ASCII vectors)

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