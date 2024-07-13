#ifndef qtCWENO_quadTreeForest_H //include guard
#define qtCWENO_quadTreeForest_H

#include <vector>
#include <fstream>

#include "quadTree.h"
#include "nodeEdge.h"

class quadTreeForest {
    std::vector<quadTree> trees; //деревья
    std::vector<std::vector<nodeTag>> toRefine; //поуровневый список нод, подлежащих балансировке

public: //векторы сделаны публичными для упрощения циклов по уровням и ячейкам (TODO: разобраться, как сделать - и нужно ли - их приватными)
    //std::vector<nodeEdge> edges; //стороны ячеек
    //std::vector<nodeEdgeId> vacant_edge_ids; //список свободных мест в edges

    void initialize(); //выделение памяти под деревья и прочее
    void addTree(quadTree tree); //доабвление дерева в список
    quadTree& getTreeByCoords(point p); //поиск дерева по координатам точки
    //edgeId getVacantEdgeId(); //получение номера вакантной ячейки или содание новой в векторе edges
    //edgeId addEdge(nodeEdge edge); //внесение ребра в список
    //edgeId addEdgeUnique(nodeEdge edge); //внесение ребра в список с проверкой уникальности
    //int updateEdge(edgeId eid, nodeTag n1, nodeTag n2); //обновление данных ребра
    //int removeEdge(edgeId eid); //удаление ребра из списка
    dataExtrema getExtrema(); //сбор экстремумов всех величин для вывода в Tecplot
    size_t activeNodesNumber(); //подсчет активных (неудаленных) нод
    size_t leavesNumber(); //подсчет листьев

    void exportForestScatter(std::string filename); //вывод в файл (Tecplot ASCII scatter)

    friend class quadTree;
    friend class treeNode;


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