#include "main.h"
#include "cellBox.h"
#include "cellData.h"
#include "treeNode.h"
#include "quadTree.h"

#include <sstream>


quadTree::quadTree(quadTreeId _id)  //конструктор дерева с одной нодой
{
    id = _id;
    //нулевой уровень не используется для простоты нумерации
    nodes.emplace_back();
    active_nodes_num.push_back(1); //на нулевом уровне хранится общее число активных нод
    leaf_nodes_num.push_back(1); //общее число листьев
    vacant_node_ids.emplace_back();
    //первый слой - одна нода
    initNewLevel(); //выделение места под первый слой
    nodes[FIRST_LEVEL].emplace_back(); //выделение места под одну ноду (и инициализация значениями по умолчанию)
    treeNode& root = nodes[FIRST_LEVEL][FIRST_ID]; //по ссылке
    root.setTag({ id, FIRST_LEVEL, FIRST_ID });
    cellBox root_box = generateBox(config.global_box);
    root.setBox(root_box);
    //root.is_leaf = true; (=true по умолчанию)
    //root.is_deleted = false; (=false по умолчанию)
    //root.parent = null; (=null по умолчанию)
    active_nodes_num[FIRST_LEVEL] = 1; //счетчик активных нод
    leaf_nodes_num[FIRST_LEVEL] = 1; //счетчик листовых нод

    //nodeTag ctag = { root.tree_id, root.depth, root.id }; //тег корневой ноды
    //соседи и ребра
    for (auto n : Neighbours)
    {
        auto nid = calcNeighbourTreeId(n);
        if (nid == null)
            is_ghost = true;
        else
        {
            //запись о соседе (корень соседнего дерева)
            root.setNeighbour(n, { nid, FIRST_LEVEL, FIRST_ID });
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
    /*for (auto n : Neighbours12)
    {
        auto nid = getNeighbour12TreeId(n);
        root.neighbours12[n] = { nid, FIRST_LEVEL, FIRST_ID };
    }*/

    //выделение памяти под физические данные
    root.setDataId(getVacantDataId());
    if (config.coord_type == CoordType::axisymmetric) //внесение r в cellData для осесимметричных координат
    {
        data[root.dataId()].setY(root_box.center().y);
    }
}

void quadTree::initNewLevel()  //инициализация нового уровня дерева (cross-check with coarsenTreeNode)
{
    std::vector<treeNode> nodes_level;
    nodes.push_back(nodes_level); //новый слой нод
    active_nodes_num.push_back(0); //новая запись о числе активных нод на слое
    leaf_nodes_num.push_back(0); //новая запись о числе листьев на слое
    std::vector<treeNodeId> vacant_node_ids_level;
    vacant_node_ids.push_back(vacant_node_ids_level); //новый слой записей о вакантных нодых
    depth++; //глубина дерева +1
    return;
}

cellBox quadTree::generateBox(cellBox global_box) const //вычисление bounding box для дерева по параметрам начальной сетки и индексу дерева
{
    double dx = (global_box.bottom_right().x - global_box.bottom_left().x) / config.Nx;
    double dy = (global_box.top_left().y - global_box.bottom_left().y) / config.Ny;
    size_t discrete_coord_x = (size_t)(id / config.Ny) % config.Nx;
    size_t discrete_coord_y = (size_t)id % config.Ny;
    return cellBox({ global_box.bottom_left().x + discrete_coord_x * dx,global_box.bottom_left().y + discrete_coord_y * dy }, 
        { global_box.bottom_left().x + (discrete_coord_x + 1) * dx,global_box.bottom_left().y + (discrete_coord_y + 1) * dy });
}

quadTreeId quadTree::calcNeighbourTreeId(Neighbour Neighbour) const //вычисление id соседнего дерева по id текущего
{
    quadTreeId ret = null;
    size_t discrete_coord_x = (size_t)(id / config.Ny) % config.Nx;
    size_t discrete_coord_y = id % config.Ny;
    switch (Neighbour)
    {
    case Neighbour::top:
        if (discrete_coord_y < config.Ny - 1)
            ret = id + 1;
        break;
    case Neighbour::right:
        if (discrete_coord_x < config.Nx - 1)
            ret = id + config.Ny;
        break;
    case Neighbour::bottom:
        if (discrete_coord_y > 0)
            ret = id - 1;
        break;
    case Neighbour::left:
        if (discrete_coord_x > 0)
            ret = id - config.Ny;
        break;
    }
    return ret;
}

cellDataId quadTree::getVacantDataId() //получение номера вакантной ячейки или создание новой в векторе data
{
    cellDataId ret;
    //если нет вакантных мест - писать в конец массива
    if (vacant_data_ids.empty())
    {
        ret = data.size(); //на 1 больше номера последнего элемента
        data.emplace_back(); //выделение памяти (должно вызываться после .size())
    }
    else
    {
        ret = vacant_data_ids.back(); //номер последней (по времени маркирования) вакантной ячейки в массиве data
        vacant_data_ids.pop_back(); //убрать из списка вакантных ячеек
    }
    return ret;
}