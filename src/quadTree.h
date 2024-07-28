#ifndef qtCWENO_QuadTree_H //include guard
#define qtCWENO_QuadTree_H

#include "main.h"
#include "TreeNode.h"

class QuadTree
{
    const int FIRST_LEVEL = 1; //вместо define (static по умолчанию?)
    const treeNodeId FIRST_ID = 0;
    quadTreeId id_;
    bool is_ghost = false; //является ли дерево ghost'ом
    int depth_ = 0;

    std::vector<std::vector<TreeNode>> nodes; //поуровневая структура дерева
    std::vector<std::vector<treeNodeId>> vacant_node_ids; //поуровневый список свободных мест в nodes
    std::vector<size_t> active_nodes_num; //число активных (неудаленных) нод по уровням
    std::vector<size_t> leaf_nodes_num; //число листовых нод по уровням
    std::vector<CellData> data_; //одномерный вектор физических данных для листьев
    std::vector<cellDataId> vacant_data_ids; //список свободных мест в data

public:
    explicit QuadTree(quadTreeId _id); //конструктор дерева с одной нодой

    //accessors
    treeNodeId id() const;
    int depth() const;
    TreeNode& rootRef(); //ссылка на корневую ноду
    const TreeNode& rootRefConst() const; //ссылка на (const) корневую ноду
    TreeNode& nodeRef(int depth, treeNodeId id); //ссылка на ноду
    CellData& dataRef(cellDataId id); //ссылка на CellData по id
    CellData data(cellDataId id) const; //копия CellData по id
    TreeNode& getNodeByCoords(Point p) const; //ссылка на ноду по координатам

    //mutators
    void initNewLevel(); //инициализация нового уровня дерева (cross-check with coarsenTreeNode)
    void deleteLevelIfEmpty(int depth); //удаление уровня, если на нем не осталось ячеек
    void vacateNodeGroup(int depth, treeNodeId id); //пометка группы нод свободной
    void vacateData(cellDataId did); //пометка ячейки данных свободной
    void incrementCounterNodes(int dpth, int amount); //изменение счетчика активных нод
    void incrementCounterLeaves(int dpth, int amount); //изменение счетчика листьев

    //inspectors
    bool isGhost() const; //является ли дерево ghost'ом
    bool isGhostCorner() const; //является ли дерево угловым ghost'ом
    //static bool isTreeGhost(quadTreeId id); //является ли дерево ghost'ом (по id)

    //other
    CellBox generateBox(const CellBox& global_box) const; //вычисление bounding box для дерева по параметрам начальной сетки и индексу дерева
    quadTreeId calcNeighbourTreeId(Neighbour n) const; //вычисление id соседнего дерева по id текущего
    quadTreeId calcNeighbour12TreeId(Neighbour12 n12) const; //вычисление id соседнего дерева по id текущего
    treeNodeId getVacantNodeId(int depth); //получение номера вакантной ячейки или создание новой в векторе nodes[depth]
    cellDataId getVacantDataId(); //получение номера вакантной ячейки или создание новой в векторе data


    //output
    const std::string dump() const; //дамп дерева в строку

    friend class QuadTreeForest;
    friend class TreeNode;

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
