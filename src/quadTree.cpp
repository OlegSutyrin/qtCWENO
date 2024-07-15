#include "main.h"
#include "CellBox.h"
#include "CellData.h"
#include "TreeNode.h"
#include "QuadTree.h"

#include <sstream>


QuadTree::QuadTree(quadTreeId _id) //конструктор дерева с одной нодой
{
    id_ = _id;
    //нулевой уровень не используется для простоты нумерации
    nodes.emplace_back();
    active_nodes_num.push_back(1); //на нулевом уровне хранится общее число активных нод
    leaf_nodes_num.push_back(1); //общее число листьев
    vacant_node_ids.emplace_back();
    //первый слой - одна нода
    initNewLevel(); //выделение места под первый слой
    nodes[FIRST_LEVEL].emplace_back(); //выделение места под одну ноду (и инициализация значениями по умолчанию)
    TreeNode& root = nodes[FIRST_LEVEL][FIRST_ID]; //по ссылке
    root.setTag(NodeTag(id_, FIRST_LEVEL, FIRST_ID));
    CellBox root_box = generateBox(config.global_box);
    root.setBox(root_box);
    //root.is_leaf = true; (=true по умолчанию)
    //root.is_deleted = false; (=false по умолчанию)
    //root.parent = null; (=null по умолчанию)
    active_nodes_num[FIRST_LEVEL] = 1; //счетчик активных нод
    leaf_nodes_num[FIRST_LEVEL] = 1; //счетчик листовых нод

    //NodeTag ctag = { root.tree_id, root.depth, root.id }; //тег корневой ноды
    //соседи и ребра
    for (auto n : Neighbours)
    {
        auto nid = calcNeighbourTreeId(n);
        if (nid == null)
            is_ghost = true;
        else
        {
            //запись о соседе (корень соседнего дерева)
            root.setNeighbour(n, NodeTag(nid, FIRST_LEVEL, FIRST_ID));
            //создание ребра между ячейками
            /*nodeEdge edge;
            switch (n)
            {
            case NEIGHBOUR_TOP:
                edge = { root.neighbours[n], ctag, ORIENTATION_HORIZONTAL };
                break;
            case NEIGHBOUR_RIGHT:
                edge = { ctag, root.neighbours[n], ORIENTATION_VERTICAL };
                break;
            case NEIGHBOUR_BOTTOM:
                edge = { ctag, root.neighbours[n], ORIENTATION_HORIZONTAL };
                break;
            case NEIGHBOUR_LEFT:
                edge = { root.neighbours[n], ctag, ORIENTATION_VERTICAL };
                break;
            }*/
            //добавление ребра в список с проверкой уникальности
            //root.edges[2 * n] = forest.addEdgeUnique(edge); //умножение на 2 из-за нумерации ребер
            //root.edges[2 * n + 1] = null; //"вторая половина" неразделенной стороны
            //cout << "created edge " << forest.edges[root.edges[2 * n]] << endl;
        }
    }

    //соседи с учетом диагональных
    for (auto n : Neighbours12)
    {
        auto nid = calcNeighbour12TreeId(n);
        root.setNeighbour12(n, NodeTag(nid, FIRST_LEVEL, FIRST_ID));
    }

    //выделение памяти под физические данные
    root.setDataId(getVacantDataId());
    if (config.coord_type == CoordType::axisymmetric) //внесение r в CellData для осесимметричных координат
    {
        data_[root.dataId()].setY(root_box.center().y);
    }
}

treeNodeId QuadTree::id() const { return id_; }
bool QuadTree::isGhost() const { return is_ghost; }

const TreeNode& QuadTree::root() const { return nodes[FIRST_LEVEL][FIRST_ID]; } //ссылка на (неизменяюему) корневую ноду
const int ERROR_NO_NODE = -1;
TreeNode& QuadTree::nodeRef(int d, treeNodeId id)
{
    try {
        if (d <= depth && id < nodes[d].size()) //есть такая нода (м.б. удаленная)
        {
            return nodes[d][id];
        }
        else
        {
            throw ERROR_NO_NODE;
        }
    }
    catch (int err_code)
    {
        cout << "tree.nodeRef() error: " << err_code << "tag: (" << id_ << ", " << d << ", " << id << ")" << endl;
        return nodes[FIRST_LEVEL][FIRST_ID]; //to suppress warning
    }
}

const int ERROR_BAD_DATA_ID = -1;
CellData& QuadTree::dataRef(cellDataId id) //ссылка на CellData по id
{
    try {
        if (id < data_.size()) //есть такие данные (м.б. удаленные)
        {
            return data_[id];
        }
        else
        {
            throw ERROR_BAD_DATA_ID;
        }
    }
    catch (int err_code)
    {
        cout << "tree.dataRef() error: " << err_code << "data id: " << id << endl;
        return data_[0]; //to suppress warning
    }
}
CellData QuadTree::data(cellDataId id) //копия CellData по id
{
    try {
        if (id < data_.size()) //есть такие данные (м.б. удаленные)
        {
            return data_[id];
        }
        else
        {
            throw ERROR_BAD_DATA_ID;
        }
    }
    catch (int err_code)
    {
        cout << "tree.data() error: " << err_code << "data id: " << id << endl;
        return data_[0]; //to suppress warning
    }
}

void QuadTree::initNewLevel()  //инициализация нового уровня дерева (cross-check with coarsenTreeNode)
{
    std::vector<TreeNode> nodes_level;
    nodes.push_back(nodes_level); //новый слой нод
    active_nodes_num.push_back(0); //новая запись о числе активных нод на слое
    leaf_nodes_num.push_back(0); //новая запись о числе листьев на слое
    std::vector<treeNodeId> vacant_node_ids_level;
    vacant_node_ids.push_back(vacant_node_ids_level); //новый слой записей о вакантных нодах
    depth++; //глубина дерева +1
    return;
}

CellBox QuadTree::generateBox(const CellBox& global_box) const //вычисление bounding box для дерева по параметрам начальной сетки и индексу дерева
{
    double dx = (global_box.bottom_right().x - global_box.bottom_left().x) / config.Nx;
    double dy = (global_box.top_left().y - global_box.bottom_left().y) / config.Ny;
    size_t discrete_coord_x = (size_t)(id_ / config.Ny) % config.Nx;
    size_t discrete_coord_y = (size_t)id_ % config.Ny;
    return CellBox({ global_box.bottom_left().x + discrete_coord_x * dx,global_box.bottom_left().y + discrete_coord_y * dy }, 
        { global_box.bottom_left().x + (discrete_coord_x + 1) * dx,global_box.bottom_left().y + (discrete_coord_y + 1) * dy });
}

quadTreeId QuadTree::calcNeighbourTreeId(Neighbour n) const //вычисление id соседнего дерева по id текущего
{
    quadTreeId ret = null;
    size_t discrete_coord_x = (size_t)(id_ / config.Ny) % config.Nx;
    size_t discrete_coord_y = id_ % config.Ny;
    switch (n)
    {
    case Neighbour::top:
        if (discrete_coord_y < config.Ny - 1)
            ret = id_ + 1;
        break;
    case Neighbour::right:
        if (discrete_coord_x < config.Nx - 1)
            ret = id_ + config.Ny;
        break;
    case Neighbour::bottom:
        if (discrete_coord_y > 0)
            ret = id_ - 1;
        break;
    case Neighbour::left:
        if (discrete_coord_x > 0)
            ret = id_ - config.Ny;
        break;
    }
    return ret;
}
quadTreeId QuadTree::calcNeighbour12TreeId(Neighbour12 n12) const
{
    quadTreeId ret = null;
    size_t discrete_coord_x = (size_t)(id_ / config.Ny) % config.Nx;
    size_t discrete_coord_y = id_ % config.Ny;
    switch (n12)
    {
    case Neighbour12::top1: //единственный верхний записывается как top1 = <tag>, top2 = null
        if (discrete_coord_y < config.Ny - 1)
            ret = id_ + 1;
        break;
    case Neighbour12::top_right:
        if (discrete_coord_y < config.Ny - 1 && discrete_coord_x < config.Nx - 1)
            ret = id_ + 1 + config.Ny;
        break;
    case Neighbour12::right1: //единственный правый
        if (discrete_coord_x < config.Nx - 1)
            ret = id_ + config.Ny;
        break;
    case Neighbour12::bottom_right:
        if (discrete_coord_y > 0 && discrete_coord_x < config.Nx - 1)
            ret = id_ - 1 + config.Ny;
        break;
    case Neighbour12::bottom1: //единственный нижний
        if (discrete_coord_y > 0)
            ret = id_ - 1;
        break;
    case Neighbour12::bottom_left:
        if (discrete_coord_y > 0 && discrete_coord_x > 0)
            ret = id_ - 1 - config.Ny;
        break;
    case Neighbour12::left1: //единственный левый
        if (discrete_coord_x > 0)
            ret = id_ - config.Ny;
        break;
    case Neighbour12::top_left:
        if (discrete_coord_y < config.Ny - 1 && discrete_coord_x > 0)
            ret = id_ + 1 - config.Ny;
        break;
    default:
        ret = null;
    }
    return ret;
}

treeNodeId QuadTree::getVacantNodeId(int depth) //получение номера вакантной ячейки или создание новой в векторе nodes[depth]
{
    treeNodeId ret;
    //если нет вакантных мест - писать в конец массива
    if (vacant_node_ids[depth].empty())
    {
        ret = nodes[depth].size(); //на 1 больше номера последнего элемента
        for (auto q : Quadrants) //можно сделать одним вызовом?
            nodes[depth].emplace_back(); //выделение памяти, должно вызываться после .size() выше
    }
    else
    {
        ret = vacant_node_ids[depth].back(); //номер последней (по порядку маркирования) вакантной ячейки в массиве nodes[depth]
        vacant_node_ids[depth].pop_back(); //убрать этот номер из списка вакантных ячеек
    }
    return ret;
}

cellDataId QuadTree::getVacantDataId() //получение номера вакантной ячейки или создание новой в векторе data
{
    cellDataId ret;
    //если нет вакантных мест - писать в конец массива
    if (vacant_data_ids.empty())
    {
        ret = data_.size(); //на 1 больше номера последнего элемента
        data_.emplace_back(); //выделение памяти (должно вызываться после .size())
    }
    else
    {
        ret = vacant_data_ids.back(); //номер последней (по времени маркирования) вакантной ячейки в массиве data
        vacant_data_ids.pop_back(); //убрать из списка вакантных ячеек
    }
    return ret;
}

TreeNode& QuadTree::getNodeByCoords(Point p) const //ссылка на ноду по координатам
{
    auto& tree = forest.getTreeByCoords(p);
    return tree.nodes[FIRST_LEVEL][FIRST_ID].getChildByCoords(p);
}


const std::string QuadTree::dump() const //дамп дерева в строку (без нод)
{
    std::stringstream ret;
    ret << "Tree " << id_ << ", " << "depth " << depth << " ( ";
    for (auto& n : active_nodes_num)
    {
        ret << n << " ";
    }
    ret << "):" << endl;
    return ret.str();
}
