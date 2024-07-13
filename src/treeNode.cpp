#include "main.h"
#include "cellBox.h"
#include "cellData.h"
#include "nodeEdge.h"
#include "nodeTag.h"
#include "quadTree.h"
#include "treeNode.h"

#include "Eigen\Dense"
#include <iomanip>

bool treeNode::isDeleted() const { return is_deleted; }
bool treeNode::isLeaf() const { return is_leaf; }
nodeTag treeNode::tag() const { return tag_; } //получение тэга 
void treeNode::setTag(nodeTag t) { tag_ = t; } //задание тэга
cellDataId treeNode::dataId() const { return dataId_; } //получение dataId
void treeNode::setDataId(cellDataId id) { dataId_ = id; } //задание dataId
cellBox treeNode::box() const {	return box_; } //получение box'a
void treeNode::setBox(cellBox b) { box_ = b; } //задание box'а
nodeTag treeNode::getNeighbour(Neighbour n) const {	return neighbours[static_cast<int>(n)]; }
void treeNode::setNeighbour(Neighbour n, nodeTag t) { neighbours[static_cast<int>(n)] = t; }

/*
const int ERROR_NODE_NO_DATA = -1;
cellData& treeNode::dataRef() const //ссылка на данные
{
    try {
        if (dataId_ != null)
            return forest.dataRef(tag_.tree(), dataId_);
        else
            throw ERROR_NODE_NO_DATA;
    }
    catch (int err_code)
    {
        cout << "getDataRef error: " << err_code << ", tag = " << tag_ << endl;
        return forest.dataRef(0, 0); //to suppress warning
    }
}

double treeNode::magGradRho() const //примерный градиент плотности
{
    if (!is_leaf || is_deleted)
        return 0.0;

    cellData& data = getDataRef();
    if (data.rho() < DOUBLE_EPS10)
        return 0.0;
    double drhosum = 0.0;
    int n12num = 0;
    for (auto n12 : Neighbours12)
    {
        if (hasNeighbour12(n12) && !quadTree::isTreeGhost(neighbours12[n12].tree_id))
        {
            auto& nnode = getNode(neighbours12[n12]);
            auto& ndata = nnode.getDataRef(); //по ссылке, т.к. в neighbour12 лежат только листья
            drhosum += fabs(data.rho() - ndata.rho()) / distance(box.center, nnode.box.center);
            n12num++;
        }
    }
    return drhosum / n12num / data.rho();
    return 0.0;
}*/

const int ERROR_NODE_NO_DATA = -1;
cellData& treeNode::dataRef() const //ссылка на данные
{
    try {
        if (dataId_ != null)
            return forest.trees[tag().tree()].data[dataId_];
        else
            throw ERROR_NODE_NO_DATA;
    }
    catch (int err_code)
    {
        cout << "getDataRef error: " << err_code << ", tag = " << tag() << endl;
        return forest.trees[0].data[0]; //to suppress warning
    }
}

const int ERROR_NODE_NO_CHILDREN = -1;
treeNode& treeNode::getChild(Quadrant q) //ссылка на ребенка по квадранту
{
    try {
        if (childrenId_ != null)
        {
            return forest.trees[tag().tree()].nodes[tag().depth() + 1][childrenId_ + static_cast<int>(q)];
        }
        else
        {
            throw ERROR_NODE_NO_CHILDREN;
        }
    }
    catch (int err_code)
    {
        cout << "getChild error: " << err_code << " " << tag().tree() << " " << tag().depth() << " " << childrenId_ << endl;
        return *this; //to suppress warning
    }
}

treeNode& treeNode::getChildByCoords(point p) //ссылка на ребенка по координатам
{
    if (is_leaf)
        return *this;
    for (auto q : Quadrants)
    {
        auto& child = getChild(q);
        if (child.box().isPointInside(p))
            return child.getChildByCoords(p);
    }
    return forest.trees[0].nodes[1][0]; //to suppress warning
}

double treeNode::magGradRho() const //примерный градиент плотности
{
    if (!is_leaf || is_deleted)
        return 0.0;

    cellData& data = dataRef();
    if (data.rho() < DOUBLE_EPS12)
        return 0.0;
    double drhosum = 0.0;
    int n12num = 0;
    /*for (auto n12 : Neighbours12)
    {
        if (hasNeighbour12(n12) && !quadTree::isTreeGhost(neighbours12[n12].tree_id))
        {
            auto& nnode = getNode(neighbours12[n12]);
            auto& ndata = nnode.getDataRef(); //по ссылке, т.к. в neighbour12 лежат только листья
            drhosum += fabs(data.rho() - ndata.rho()) / distance(box.center, nnode.box.center);
            n12num++;
        }
    }*/


    
    
    return 0.0;
    


    
    //return drhosum / n12num / data.rho();
}

std::string treeNode::dump() const //дамп ноды в строку
{
    std::stringstream buffer;
    buffer << "(" << std::setw(3) << (tag_.tree() == null ? " - " : std::to_string(tag_.tree())) << ", " << tag_.depth() << ", " << std::setw(3) << (tag_.id() == null ? " - " : std::to_string(tag_.id())) << ", [" << std::setw(3) << (parentId_ == null ? " - " : std::to_string(parentId_)) << "]): ";
    /*for (auto n : Neighbours)
    {
        if (neighbours[static_cast<int>(n)].id() != null)
            buffer << Dirs[n] << ": " << std::setw(3) << neighbours[n].tree_id << ", " << neighbours[n].depth << ", " << std::setw(3) << neighbours[n].id << " | ";
        else
            buffer << "               | ";
    }*/
    if (is_deleted)
        buffer << " deleted";
    buffer << endl;

    /*buffer << "                     [";
    for (auto n : Neighbours)
    {
        if (neighbours[n].tree_id != null && neighbours[n].id != null)
        {
            auto& nnode = getNode(neighbours[n]);
            Neighbour dir = (Neighbour)((int)(n + 2) % 4); //противоположное направление
            buffer << Dirs[dir] << ": " << std::setw(3) << nnode.neighbours[dir].tree_id << ", " << nnode.neighbours[dir].depth << ", " << std::setw(3) << nnode.neighbours[dir].id << " | ";
        }
        else
            buffer << "               | ";
    }
    buffer << "]" << endl;*/

    return buffer.str();
}