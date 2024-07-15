#include "main.h"
#include "CellBox.h"
#include "CellData.h"
#include "NodeEdge.h"
#include "NodeTag.h"
#include "QuadTree.h"
#include "TreeNode.h"

#include "Eigen\Dense"
#include <iomanip>

//accessors -----------------------
NodeTag TreeNode::tag() const { return tag_; } //получение тэга 
cellDataId TreeNode::dataId() const { return dataId_; } //получение dataId
CellBox TreeNode::box() const {	return box_; } //получение box'a
NodeTag TreeNode::getNeighbour(Neighbour n) const {	return neighbours[static_cast<int>(n)]; }
NodeTag TreeNode::getNeighbour12(Neighbour n12) const { return neighbours12[static_cast<int>(n12)]; }

//const int ERROR_NODE_NO_DATA = -1;
CellData& TreeNode::dataRef() const //ссылка на данные
{
    try {
        if (dataId_ != null)
            return forest.treeRef(tag().tree()).dataRef(dataId_);
        else
            throw std::invalid_argument("null data reference in this node!");
    }
    catch (const std::invalid_argument& e)
    {
        cout << "getDataRef error: " << e.what() << ", tag = " << tag() << endl;
        return forest.treeRef(0).dataRef(0); //to suppress warning
    }
}

CellData TreeNode::data() const //копия данных
{
    try {
        if (is_leaf)
        {
            if (dataId_ != null)
            {
                return forest.treeRef(tag().tree()).data(dataId_);
            }
            else
            {
                throw std::invalid_argument("null data reference in this node!");
            }
        }
        else //нужно собрать данные с детей
        {
            CellData d; //инициализируется нулями
            for (auto q : Quadrants) //рекурсивный консервативный сбор данных с детей
            {
                d.add(childRef(q).data());
            }
            d.divide(static_cast<double>(QUADRANTS_NUM));
            return d;
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "getData error: " << e.what() << ", tag = " << tag() << endl;
        return forest.treeRef(0).data(0); //to suppress warning
    }
}

const int ERROR_NODE_NO_CHILDREN = -1;
TreeNode& TreeNode::childRef(Quadrant q) //ссылка на ребенка по квадранту
{
    try {
        if (childrenId_ != null)
        {
            NodeTag t = NodeTag(tag().tree(), tag().depth() + 1, (treeNodeId)(childrenId_ + static_cast<int>(q)));
            return nodeRef(t);
        }
        else
        {
            throw ERROR_NODE_NO_CHILDREN;
        }
    }
    catch (int err_code)
    {
        cout << "node.childRef error: " << err_code << " " << tag().tree() << " " << tag().depth() + 1 << " " << childrenId_ + static_cast<int>(q) << endl;
        return *this; //to suppress warning
    }
}

const TreeNode& TreeNode::childRef(Quadrant q) const //(const) ссылка на ребенка для функции TreeNode::data() const
{
    try {
        if (childrenId_ != null)
        {
            NodeTag t = NodeTag(tag().tree(), tag().depth() + 1, (treeNodeId)(childrenId_ + static_cast<int>(q)));
            return nodeRef(t);
        }
        else
        {
            throw ERROR_NODE_NO_CHILDREN;
        }
    }
    catch (int err_code)
    {
        cout << "node.childRef error: " << err_code << " " << tag().tree() << " " << tag().depth() + 1 << " " << childrenId_ + static_cast<int>(q) << endl;
        return *this; //to suppress warning
    }
}

TreeNode& TreeNode::getChildByCoords(Point p) //ссылка на ребенка по координатам
{
    if (is_leaf)
        return *this;
    for (auto q : Quadrants)
    {
        auto& rchild = childRef(q);
        if (rchild.box().isPointInside(p))
            return rchild.getChildByCoords(p);
    }
    return nodeRef({}); //to suppress warning
}

//mutators -----------------------
void TreeNode::setTag(const NodeTag& t) { tag_ = t; } //задание тэга
void TreeNode::setDataId(cellDataId id) { dataId_ = id; } //задание dataId
void TreeNode::setBox(const CellBox& b) { box_ = b; } //задание box'а
void TreeNode::setNeighbour(Neighbour n, NodeTag t) { neighbours[static_cast<int>(n)] = t; }
void TreeNode::setNeighbour12(Neighbour12 n12, NodeTag t) { neighbours12[static_cast<int>(n12)] = t; }

void TreeNode::setData(const CellData& data) //запись данных
{
    try {
        if (dataId_ == null)
        {
            throw std::invalid_argument("null data reference in this node!");
        }
        else
        {
            dataRef().set(data);
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "setData error: " << e.what() << ", tag = " << tag() << endl;
    }
    return;
}

//inspectors -----------------------
bool TreeNode::isDeleted() const { return is_deleted; }
bool TreeNode::isLeaf() const { return is_leaf; }
bool TreeNode::hasChildren() const { return childrenId_ != null; } //есть ли дети
bool TreeNode::hasGrandChildren() const { return has_grandchildren; } //есть ли внуки

//uncategorized -----------------------
const int ERROR_NODE_IS_DELETED = -1;
TreeNode& TreeNode::nodeRef(const NodeTag& tag) //ссылка на ноду по тэгу
{
    try {
        TreeNode& rnode = forest.treeRef(tag.tree()).nodeRef(tag.depth(), tag.id());
        if(rnode.isDeleted())
        {
            throw ERROR_NODE_IS_DELETED;
        }
        return rnode;
    }
    catch (int err_code)
    {
        cout << "TreeNode.nodeRef error: " << err_code << ", tag: " << tag << endl;
    }
    return nodeRef({}); //to suppress warning
}

double TreeNode::magGradRho() const //примерный градиент плотности
{
    if (!is_leaf || is_deleted)
        return 0.0;

    CellData& rd = dataRef();
    if (rd.rho() < DOUBLE_EPS12)
        return 0.0;
    double drhosum = 0.0;
    int n12num = 0;
    for (auto& n12 : neighbours12)
    {
        if (!n12.isNull() && !forest.treeRef(n12.tree()).isGhost())
        {
            auto& rnnode = nodeRef(n12);
            auto& rndata = rnnode.dataRef(); //можно по ссылке, т.к. в neighbour12 лежат только листья
            drhosum += fabs(rd.rho() - rndata.rho()) / distance(box().center(), rnnode.box().center());
            n12num++;
        }
    }
    return drhosum / n12num / rd.rho();
}

//output -----------------------
std::string TreeNode::dump() const //дамп ноды в строку
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