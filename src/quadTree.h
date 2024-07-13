#ifndef qtCWENO_quadTree_H //include guard
#define qtCWENO_quadTree_H

#include "main.h"
#include "treeNode.h"

class quadTree
{
    const int FIRST_LEVEL = 1; //вместо define (static по умолчанию?)
    const treeNodeId FIRST_ID = 0;
    quadTreeId id_;
    bool is_ghost = false; //является ли дерево ghost'ом
    int depth = 0;

    std::vector<std::vector<treeNode>> nodes; //поуровневая структура дерева
    std::vector<std::vector<treeNodeId>> vacant_node_ids; //поуровневый список свободных мест в nodes
    std::vector<size_t> active_nodes_num; //число активных (неудаленных) нод по уровням
    std::vector<size_t> leaf_nodes_num; //число листовых нод по уровням
    std::vector<cellData> data; //одномерный вектор физических данных для листьев
    std::vector<cellDataId> vacant_data_ids; //список свободных мест в data

public:
    quadTree(quadTreeId _id); //конструктор дерева с одной нодой
    treeNodeId id() const;
    const treeNode& root() const; //ссылка на (неизменяюему) корневую ноду
    bool isGhost() const; //является ли дерево ghost'ом
    bool isGhostCorner() const; //является ли дерево угловым ghost'ом
    static bool isTreeGhost(quadTreeId id); //является ли дерево ghost'ом (по id)
    void initNewLevel(); //инициализация нового уровня дерева (cross-check with coarsenTreeNode)
    cellBox generateBox(cellBox global_box) const; //вычисление bounding box для дерева по параметрам начальной сетки и индексу дерева
    quadTreeId calcNeighbourTreeId(Neighbour Neighbour) const; //вычисление id соседнего дерева по id текущего
    //quadTreeId getNeighbour12TreeId(Neighbour12 Neighbour); //вычисление id соседнего дерева по id текущего
    treeNodeId getVacantNodeId(int depth); //получение номера вакантной ячейки или создание новой в векторе nodes[depth]
    cellDataId getVacantDataId(); //получение номера вакантной ячейки или создание новой в векторе data
    treeNode& getNodeByCoords(point p) const; //ссылка на ноду по координатам
    const std::string dump() const; //дамп дерева в строку

    friend class quadTreeForest;
    friend class treeNode;

    //шаблон цикла по нодам
    /*template<typename Func>
    inline void forAllNodes(Func f, bool include_branches = SKIP_BRANCHES)
    {
        for (auto& nodes_level : nodes) //по ссылке
        {
            for (auto& cnode : nodes_level) //по ссылке
            {
                if (cnode.isDeleted() || !include_branches && !cnode.isLeaf())
                    continue;
                f(cnode);
            }
        }
    }*/
};

#endif
