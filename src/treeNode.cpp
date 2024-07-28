#include "main.h"
#include "CellBox.h"
#include "CellData.h"
#include "NodeEdge.h"
#include "NodeTag.h"
#include "QuadTree.h"
#include "TreeNode.h"

#include "Eigen\Dense"
#include <iomanip>

//------------------------------------------------------------------- ChildrenTags -------------------------------------------------------------------
const NodeTag& ChildrenTags::operator()(Quadrant q) const { return tags[static_cast<int>(q)]; } //() operator: тэг по квадранту

//------------------------------------------------------------------- TreeNode -------------------------------------------------------------------
const NodeTag TreeNode::tag() const { return tag_; } //получение тэга 
cellDataId TreeNode::dataId() const { return dataId_; } //получение dataId
CellBox TreeNode::box() const {	return box_; } //получение box'a
//QuadTree& TreeNode::treeRef() { return forest.treeRef(tag().tree()); } // ссылка на дерево, содержущее ноду
//const QuadTree& TreeNode::treeRefConst() const { return forest.treeRefConst(tag().tree()); } //(const) ссылка на дерево, содержущее ноду
const NodeTag TreeNode::neighbour(Neighbour n) const { return neighbours[static_cast<int>(n)]; }
const NodeTag TreeNode::neighbour12(Neighbour12 n12) const { return neighbours12[static_cast<int>(n12)]; }

CellData& TreeNode::dataRef() //ссылка на данные
{
    try {
        if (dataId_ != null)
        {
            return forest.treeRef(tag().tree()).dataRef(dataId_);
            //return treeRef().dataRef(dataId_);
        }
        else
        {
            throw std::invalid_argument("null data reference in this node");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "TreeNode.dataRef() error: " << e.what() << ", tag = " << tag() << endl;
        return forest.treeRef(0).dataRef(0); //to suppress warning
    }
}
const CellData& TreeNode::dataRefConst() const
{
    try {
        if (dataId_ != null)
        {
            return forest.treeRef(tag().tree()).dataRef(dataId_);
            //return treeRef().dataRef(dataId_);
        }
        else
        {
            throw std::invalid_argument("null data reference in this node");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "TreeNode.dataRef() const error: " << e.what() << ", tag = " << tag() << endl;
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
                //return treeRefConst().data(dataId_);
            }
            else
            {
                throw std::invalid_argument("null data id in this node");
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
        cout << "TreeNode.data() error: " << e.what() << ", tag = " << tag() << endl;
        return forest.treeRef(0).data(0); //to suppress warning
    }
}

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
            throw std::invalid_argument("no children");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "node.childRef() error: " << e.what() << " " << tag().tree() << " " << tag().depth() + 1 << " " << childrenId_ + static_cast<int>(q) << endl;
        return *this; //to suppress warning
    }
}
const TreeNode& TreeNode::childRef(Quadrant q) const //const версия
{
    try {
        if (childrenId_ != null)
        {
            NodeTag t = NodeTag(tag().tree(), tag().depth() + 1, (treeNodeId)(childrenId_ + static_cast<int>(q)));
            return nodeRef(t);
        }
        else
        {
            throw std::invalid_argument("no children");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "node.childRef() const error: " << e.what() << " " << tag().tree() << " " << tag().depth() + 1 << " " << childrenId_ + static_cast<int>(q) << endl;
        return *this; //to suppress warning
    }
}
const ChildrenTags TreeNode::childrenTags() const //тэги всех детей
{
    return ChildrenTags(childRef(Quadrant::top_left).tag(), childRef(Quadrant::top_right).tag(), childRef(Quadrant::bottom_right).tag(), childRef(Quadrant::bottom_left).tag());
}

const nodeEdgeId TreeNode::edge(Edge e) const { return edges[static_cast<int>(e)]; } //id ребра
NodeEdge& TreeNode::edgeRef(Edge e) { return forest.edgeRef(edge(e)); } //ссылка на ребро
double TreeNode::polyCoeff(Equation eq, int p) const { return polyCoeffs[static_cast<int>(eq)][p]; } //коэффициент полинома CWENO

TreeNode& TreeNode::getChildOrSelfByCoords(Point p) //ссылка на ребенка по координатам
{
    if (is_leaf)
        return *this;
    for (auto q : Quadrants)
    {
        auto& rchild = childRef(q);
        if (rchild.box().isPointInside(p))
            return rchild.getChildOrSelfByCoords(p);
    }
    return nodeRef({}); //to suppress warning
}

//mutators -----------------------
void TreeNode::setTag(const NodeTag& t) { tag_ = t; } //задание тэга
void TreeNode::markDeleted() { is_deleted = true; } //пометка удаленной
void TreeNode::setDataId(cellDataId id) { dataId_ = id; } //задание dataId
void TreeNode::setBox(const CellBox& b) { box_ = b; } //задание box'а
void TreeNode::setNeighbour(Neighbour n, NodeTag ntag) //задание соседа ячейке и детям с нужной стороны
{ 
    neighbours[static_cast<int>(n)] = ntag; //запись себе
    if (hasChildren()) //детям
    {
        switch (n) 
        {
        case Neighbour::top: //пришел тэг соседа сверху
            childRef(Quadrant::top_left).setNeighbour(Neighbour::top, ntag);
            childRef(Quadrant::top_right).setNeighbour(Neighbour::top, ntag);
            break;
        case Neighbour::right: //справа
            childRef(Quadrant::top_right).setNeighbour(Neighbour::right, ntag);
            childRef(Quadrant::bottom_right).setNeighbour(Neighbour::right, ntag);
            break;
        case Neighbour::bottom: //снизу
            childRef(Quadrant::bottom_right).setNeighbour(Neighbour::bottom, ntag);
            childRef(Quadrant::bottom_left).setNeighbour(Neighbour::bottom, ntag);
            break;
        case Neighbour::left: //слева
            childRef(Quadrant::top_left).setNeighbour(Neighbour::left, ntag);
            childRef(Quadrant::bottom_left).setNeighbour(Neighbour::left, ntag);
            break;
        }

    }
}
void TreeNode::setChildrenNeighbours(Neighbour n, ChildrenTags tags) //внесение данных о соседях для детей
{
    if (!hasChildren()) //этот метод не должен вызываться для листьев
    {
        cout << "TreeNode.setChildrenNeighbours() error: no children, tag: " << tag() << endl;
        return;
    }
    switch (n)
    {
    case Neighbour::top: //данные пришли от соседа сверху
        childRef(Quadrant::top_left).setNeighbour(Neighbour::top, tags(Quadrant::bottom_left));
        childRef(Quadrant::top_right).setNeighbour(Neighbour::top, tags(Quadrant::bottom_right));
        break;
    case Neighbour::right: //справа
        childRef(Quadrant::top_right).setNeighbour(Neighbour::right, tags(Quadrant::top_left));
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::right, tags(Quadrant::bottom_left));
        break;
    case Neighbour::bottom: //снизу
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::bottom, tags(Quadrant::top_right));
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::bottom, tags(Quadrant::top_left));
        break;
    case Neighbour::left: //слева
        childRef(Quadrant::top_left).setNeighbour(Neighbour::left, tags(Quadrant::top_right));
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::left, tags(Quadrant::bottom_right));
        break;
    }
}

void TreeNode::clearNeighbours12() //удаление всех записей о соседях12 после дробления
{
    for (auto n12 : Neighbours12)
    {
        neighbours12[static_cast<int>(n12)] = NodeTag();
    }
}
void TreeNode::setNeighbour12(Neighbour12 n12, NodeTag ntag) { neighbours12[static_cast<int>(n12)] = ntag; }

void TreeNode::setChildrenOrSelfNeighbour12(Neighbour n, NodeTag ntag) //внесение данных о соседе12 для себя или детей после его склейки
{
    if (!hasChildren()) //нет детей
    {
        switch (n)
        {
        case Neighbour::top: //данные пришли от соседа сверху
            setNeighbour12(Neighbour12::top1, ntag);
            setNeighbour12(Neighbour12::top2, {}); //второго соседа с той стороны теперь нет
            break;
        case Neighbour::right: //справа
            setNeighbour12(Neighbour12::right1, ntag);
            setNeighbour12(Neighbour12::right2, {});
            break;
        case Neighbour::bottom: //снизу
            setNeighbour12(Neighbour12::bottom1, ntag);
            setNeighbour12(Neighbour12::bottom2, {});
            break;
        case Neighbour::left: //слева
            setNeighbour12(Neighbour12::left1, ntag);
            setNeighbour12(Neighbour12::left2, {});
            break;
        }
    }
    else //есть дети
    {
        switch (n)
        {
        case Neighbour::top: //данные пришли от соседа сверху
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top1, ntag); //единственный сосед по стороне
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top_right, {}); //по диагонали теперь нет соседа
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top1, ntag);
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top_left, {});
            break;
        case Neighbour::right: //справа
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::right1, ntag);
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::bottom_right, {});
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::right1, ntag);
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::top_right, {});
            break;
        case Neighbour::bottom: //снизу
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom1, ntag);
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom_left, {});
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom1, ntag);
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom_right, {});
            break;
        case Neighbour::left: //слева
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::left1, ntag);
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::bottom_left, {});
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::left1, ntag);
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::top_left, {});
            break;
        }
    }
}

void TreeNode::setNeighbours12(Neighbour n, ChildrenTags tags)
{
    switch (n)
    {
    case Neighbour::top: //данные пришли от соседа сверху
        setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_left));
        setNeighbour12(Neighbour12::top2, tags(Quadrant::bottom_right));
        break;
    case Neighbour::right: //справа
        setNeighbour12(Neighbour12::right1, tags(Quadrant::top_left));
        setNeighbour12(Neighbour12::right2, tags(Quadrant::bottom_left));
        break;
    case Neighbour::bottom: //снизу
        setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_right));
        setNeighbour12(Neighbour12::bottom2, tags(Quadrant::top_left));
        break;
    case Neighbour::left: //слева
        setNeighbour12(Neighbour12::left1, tags(Quadrant::bottom_right));
        setNeighbour12(Neighbour12::left2, tags(Quadrant::top_right));
        break;
    }
}

void TreeNode::setChildrenOrSelfNeighbours12(Neighbour n, ChildrenTags tags) //внесение данных о соседях12 для себя или детей
{
    if (!hasChildren()) //нет детей
    {
        switch (n)
        {
        case Neighbour::top: //данные пришли от соседа сверху
            setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_left));
            setNeighbour12(Neighbour12::top2, tags(Quadrant::bottom_right));
            break;
        case Neighbour::right: //справа
            setNeighbour12(Neighbour12::right1, tags(Quadrant::top_left));
            setNeighbour12(Neighbour12::right2, tags(Quadrant::bottom_left));
            break;
        case Neighbour::bottom: //снизу
            setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_right));
            setNeighbour12(Neighbour12::bottom2, tags(Quadrant::top_left));
            break;
        case Neighbour::left: //слева
            setNeighbour12(Neighbour12::left1, tags(Quadrant::bottom_right));
            setNeighbour12(Neighbour12::left2, tags(Quadrant::top_right));
            break;
        }
    }
    else //есть дети
    {
        switch (n)
        {
        case Neighbour::top: //данные пришли от соседа сверху
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_left)); //единственный сосед по стороне
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top_right, tags(Quadrant::bottom_right)); //по диагонали
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_right));
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top_left, tags(Quadrant::bottom_left));
            break;
        case Neighbour::right: //справа
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::right1, tags(Quadrant::top_left));
            childRef(Quadrant::top_right).setNeighbour12(Neighbour12::bottom_right, tags(Quadrant::bottom_left));
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::right1, tags(Quadrant::bottom_left));
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::top_right, tags(Quadrant::top_left));
            break;
        case Neighbour::bottom: //снизу
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_right));
            childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom_left, tags(Quadrant::top_left));
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_left));
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom_right, tags(Quadrant::top_right));
            break;
        case Neighbour::left: //слева
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::left1, tags(Quadrant::top_right));
            childRef(Quadrant::top_left).setNeighbour12(Neighbour12::bottom_left, tags(Quadrant::bottom_right));
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::left1, tags(Quadrant::bottom_right));
            childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::top_left, tags(Quadrant::top_right));
            break;
        }
    }
}
void TreeNode::setChildrenCommonNeighbour12(Neighbour n, NodeTag ntag) //внесение данных об общем соседе12 для детей
{
    if (!hasChildren()) //этот метод не должен вызываться для листьев
    {
        cout << "TreeNode.setChildrenCommonNeighbour12() error: no children, tag: " << tag() << endl;
        return;
    }
    switch (n)
    {
    case Neighbour::top: //данные пришли от соседа сверху
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top1, ntag); //единственный сосед по стороне
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top1, ntag);
        break;
    case Neighbour::right: //справа
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::right1, ntag);
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::right1, ntag);
        break;
    case Neighbour::bottom: //снизу
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom1, ntag);
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom1, ntag);
        break;
    case Neighbour::left: //слева
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::left1, ntag);
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::left1, ntag);
        break;
    }
}

void TreeNode::setData(const CellData& data) //запись данных
{
    try {
        if (dataId_ == null)
        {
            throw std::invalid_argument("null data reference in this node");
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
void TreeNode::setGrandParency(bool status) { has_grandchildren = status; }  //отметка о наличии внуков

void TreeNode::updateGrandParency()
{
    bool no_grandchildren = true;
    for (auto q : Quadrants)
    {
        if (childRef(q).hasChildren())
        {
            no_grandchildren = false;
            break;
        }
    }
    if (no_grandchildren)
        has_grandchildren = false;
    return;
}

/*void TreeNode::setChildrenCommonNeighbour(Neighbour n, NodeTag ntag) //внесение данных об общем соседе для детей
{
    if (!hasChildren()) //этот метод не должен вызываться для листьев
    {
        cout << "TreeNode.setChildrenCommonNeighbour() error: no children, tag: " << tag() << endl;
        return;
    }
    switch (n)
    {
    case Neighbour::top: //данные пришли от соседа сверху
        childRef(Quadrant::top_left).setNeighbour(Neighbour::top, ntag);
        childRef(Quadrant::top_right).setNeighbour(Neighbour::top, ntag);
        break;
    case Neighbour::right: //от соседа справа
        childRef(Quadrant::top_right).setNeighbour(Neighbour::right, ntag);
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::right, ntag);
        break;
    case Neighbour::bottom: //снизу
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::bottom, ntag);
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::bottom, ntag);
        break;
    case Neighbour::left: //слева
        childRef(Quadrant::top_left).setNeighbour(Neighbour::left, ntag);
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::left, ntag);
        break;
    }
}*/

void TreeNode::setEdge(Edge etype, nodeEdgeId eid) { edges[static_cast<int>(etype)] = eid; } //задание ребра
void TreeNode::splitEdge(Neighbour n) //разделить ребро и обновить записи у себя и (бездетного) соседа
{
    Edge etype = toEdge(n);
    nodeEdgeId eid = edge(etype);
    auto& rnnode = nodeRef(neighbour(n));
    ChildrenRefs crefs(childRef(Quadrant::top_left), childRef(Quadrant::top_right), childRef(Quadrant::bottom_right), childRef(Quadrant::bottom_left));
    ChildrenTags ctags(crefs(Quadrant::top_left).tag(), crefs(Quadrant::top_right).tag(), crefs(Quadrant::bottom_right).tag(), crefs(Quadrant::bottom_left).tag());
    switch (n)
    {
    case Neighbour::top: //общее ребро с соседом сверху
        forest.updateEdge(eid, neighbour(n), ctags(Quadrant::top_left)); //запись нового маленького ребра поверх большого старого
        crefs(Quadrant::top_left).setEdge(etype, eid); //запись в ребенке
        rnnode.setEdge(opposite(etype), eid); //запись в соседней ячейке
        crefs(Quadrant::top_right).setEdge(etype, forest.addEdge(NodeEdge(neighbour(n), ctags(Quadrant::top_right), Orientation::horizontal))); //создание второго ребра и запись в ребенке
        rnnode.setEdge(opposite(next(etype)), crefs(Quadrant::top_right).edge(etype)); //запись о втором ребре в соседней ячейке
        break;
    case Neighbour::right: //справа
        forest.updateEdge(eid, ctags(Quadrant::top_right), neighbour(n));
        crefs(Quadrant::top_right).setEdge(etype, eid);
        rnnode.setEdge(opposite(etype), eid);
        crefs(Quadrant::bottom_right).setEdge(etype, forest.addEdge(NodeEdge(ctags(Quadrant::bottom_right), neighbour(n), Orientation::vertical)));
        rnnode.setEdge(opposite(next(etype)), crefs(Quadrant::bottom_right).edge(etype));
        break;
    case Neighbour::bottom: //снизу
        forest.updateEdge(eid, ctags(Quadrant::bottom_right), neighbour(n));
        crefs(Quadrant::bottom_right).setEdge(etype, eid);
        rnnode.setEdge(opposite(etype), eid);
        crefs(Quadrant::bottom_left).setEdge(etype, forest.addEdge(NodeEdge(ctags(Quadrant::bottom_left), neighbour(n), Orientation::horizontal)));
        rnnode.setEdge(opposite(next(etype)), crefs(Quadrant::bottom_left).edge(etype));
        break;
    case Neighbour::left: //слева
        forest.updateEdge(eid, neighbour(n), ctags(Quadrant::bottom_left));
        crefs(Quadrant::bottom_left).setEdge(etype, eid);
        rnnode.setEdge(opposite(etype), eid);
        crefs(Quadrant::top_left).setEdge(etype, forest.addEdge(NodeEdge(neighbour(n), ctags(Quadrant::top_left), Orientation::vertical)));
        rnnode.setEdge(opposite(next(etype)), crefs(Quadrant::top_left).edge(etype));
        break;
    }
}
void TreeNode::joinEdge(Neighbour n) //склеить ребро и обновить записи у себя и (бездетного) соседа
{
    Edge etype = toEdge(n);
    nodeEdgeId eid = null; //определяется ниже
    switch (n)
    {
    case Neighbour::top: //общее ребро с соседом сверху
        forest.removeEdge(childRef(Quadrant::top_left).edge(etype)); //удаление первого ребра
        eid = childRef(Quadrant::top_right).edge(etype); //берется второе ребро, т.к. у бездетного соседа оно записано в нужном слоте
        forest.updateEdge(eid, {}, tag()); //перезапись данных в склеенном ребре
        break;
    case Neighbour::right: //справа
        forest.removeEdge(childRef(Quadrant::top_right).edge(etype));
        eid = childRef(Quadrant::bottom_right).edge(etype);
        forest.updateEdge(eid, tag(), {});
        break;
    case Neighbour::bottom: //снизу
        forest.removeEdge(childRef(Quadrant::bottom_right).edge(etype));
        eid = childRef(Quadrant::bottom_left).edge(etype);
        forest.updateEdge(eid, tag(), {});
        break;
    case Neighbour::left: //слева
        forest.removeEdge(childRef(Quadrant::bottom_left).edge(etype));
        eid = childRef(Quadrant::top_left).edge(etype);
        forest.updateEdge(eid, {}, tag());
        break;
    }
    setEdge(etype, eid); //перенос записи о ребре в первый слот
    setEdge(next(etype), null); //удаление записи о ребре во втором слоте
    nodeRef(neighbour(n)).setEdge(opposite(etype), null); //удаление записи о втором ребре у соседа
}

void TreeNode::updateChildrenEdges(Neighbour n) //обновить ребра у своих и соседских детей
{
    Edge etype = toEdge(n);
    nodeEdgeId eid = edge(etype);
    ChildrenRefs crefs(childRef(Quadrant::top_left), childRef(Quadrant::top_right), childRef(Quadrant::bottom_right), childRef(Quadrant::bottom_left));
    ChildrenTags ctags(crefs(Quadrant::top_left).tag(), crefs(Quadrant::top_right).tag(), crefs(Quadrant::bottom_right).tag(), crefs(Quadrant::bottom_left).tag());
    switch (n)
    {
    case Neighbour::top: //общее ребро с соседом сверху
        forest.updateEdge(eid, {}, ctags(Quadrant::top_left)); //частичное обновление ребра
        crefs(Quadrant::top_left).setEdge(etype, eid); //запись в ребенке
        eid = edge(next(etype));
        forest.updateEdge(eid, {}, ctags(Quadrant::top_right)); //второе ребро
        crefs(Quadrant::top_right).setEdge(etype, eid);
        break;
    case Neighbour::right: //справа
        forest.updateEdge(eid, ctags(Quadrant::top_right), {});
        crefs(Quadrant::top_right).setEdge(etype, eid);
        eid = edge(next(etype));
        forest.updateEdge(eid, ctags(Quadrant::bottom_right), {});
        crefs(Quadrant::bottom_right).setEdge(etype, eid);
        break;
    case Neighbour::bottom: //снизу
        forest.updateEdge(eid, ctags(Quadrant::bottom_right), {});
        crefs(Quadrant::bottom_right).setEdge(etype, eid);
        eid = edge(next(etype));
        forest.updateEdge(eid, ctags(Quadrant::bottom_left), {});
        crefs(Quadrant::bottom_left).setEdge(etype, eid);
        break;
    case Neighbour::left: //слева
        forest.updateEdge(eid, {}, ctags(Quadrant::bottom_left));
        crefs(Quadrant::bottom_left).setEdge(etype, eid);
        eid = edge(next(etype));
        forest.updateEdge(eid, {}, ctags(Quadrant::top_left));
        crefs(Quadrant::top_left).setEdge(etype, eid);
        break;
    }
}

void TreeNode::gatherEdgesFromChildren(Neighbour n) //записать детские ребра себе и обновить ребра у соседских детей
{
    Edge etype = toEdge(n);
    nodeEdgeId eid = null; //определяется ниже
    switch (n)
    {
    case Neighbour::top: //общее ребро с соседом сверху
        eid = childRef(Quadrant::top_left).edge(etype); //ребро ребенка
        setEdge(etype, eid); //запись себе
        forest.updateEdge(eid, {}, tag()); //обновление одного тега в ребре на тег данной ячейки
        eid = childRef(Quadrant::top_right).edge(etype); //второе ребро
        setEdge(next(etype), eid);
        forest.updateEdge(eid, {}, tag());
        break;
    case Neighbour::right: //справа
        eid = childRef(Quadrant::top_right).edge(etype);
        setEdge(etype, eid);
        forest.updateEdge(eid, tag(), {});
        eid = childRef(Quadrant::bottom_right).edge(etype);
        setEdge(next(etype), eid);
        forest.updateEdge(eid, tag(), {});
        break;
    case Neighbour::bottom: //снизу
        eid = childRef(Quadrant::bottom_right).edge(etype);
        setEdge(etype, eid);
        forest.updateEdge(eid, tag(), {});
        eid = childRef(Quadrant::bottom_left).edge(etype);
        setEdge(next(etype), eid);
        forest.updateEdge(eid, tag(), {});
        break;
    case Neighbour::left: //слева
        eid = childRef(Quadrant::bottom_left).edge(etype);
        setEdge(etype, eid);
        forest.updateEdge(eid, {}, tag());
        eid = childRef(Quadrant::top_left).edge(etype);
        setEdge(next(etype), eid);
        forest.updateEdge(eid, {}, tag());
        break;
    }
}

int TreeNode::markToRefine() //пометка ячейки к дроблению
{
    //удаленные ячейки нельзя разделять
    if (is_deleted)
        return INT_ERROR_CODE_CANT_REFINE_DELETED;

    //ячейка уже высшей глубины
    if (tag().depth() >= config.max_depth)
        return INT_ERROR_CODE_CANT_REFINE_DEEPEST_LEVEL;

    //ячейка уже разделена
    if (!is_leaf)
        return INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED;

    forest.addNodeToRefine(tag());

    //балансировка
    if(!forest.treeRef(tag().tree()).isGhost()) //для ghost-деревьев балансировка не нужна
    {
        for (auto& rntag : neighbours)
        {
            if (!rntag.isNull())
            {
                //(?)проверка на ghost-овость соседа тут не нужна, т.к. его размер выравнивается ниже
                if (rntag.depth() < tag().depth()) //сосед на 2+ уровня больше будущих детей данной ноды
                    nodeRef(rntag).markToRefine();
                //выравнивание размера соседних ghost'ов
                if (forest.treeRef(rntag.tree()).isGhost() && rntag.depth() <= tag().depth()) //ghost-сосед крупнее будущих детей данной ноды
                    nodeRef(rntag).markToRefine();
            }
        }
    }
    return 0;
}

int TreeNode::refine() //дробление ячейки
{
    //удаленные ячейки нельзя разделять
    if (is_deleted)
        return INT_ERROR_CODE_CANT_REFINE_DELETED;

    //ячейка уже высшей глубины
    if (tag().depth() >= config.max_depth)
        return INT_ERROR_CODE_CANT_REFINE_DEEPEST_LEVEL;

    //ячейка уже разделена
    if (!is_leaf)
        return INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED;

    auto& rtree = forest.treeRef(tag().tree());
    CellData cdata = this->data(); //по значению, чтобы корректно копировалось

    //выделение памяти под новый уровень (если нужно)
    if (rtree.depth() < tag().depth() + 1)
        rtree.initNewLevel();

    //запись нод-детей
    childrenId_ = rtree.getVacantNodeId(tag().depth() + 1); //получение id свободного блока мест и его запись в эту ноду
    for (auto q : Quadrants)
    {
        auto& rchild = rtree.nodeRef(tag().depth() + 1, childrenId_ + static_cast<int>(q)); //TODO: сделать генерацию детей сразу через инициализацию с нужными данными
        auto cbox = box().quarterBox(q);
        //запись через конструктор, все неуказанные поля инициализируются значениями по умолчанию
        rchild = TreeNode(
            NodeTag(tag().tree(), tag().depth() + 1, childrenId_ + static_cast<int>(q)), 
            tag().id(), 
            rtree.getVacantDataId(), 
            cbox);

        //консервативное разделение данных на детей
        if (config.coord_type == CoordType::axisymmetric) //учет r в cellData'х детей для осесиммтричных координат
            cdata.setY(cbox.center().y);
        rchild.setData(cdata);
    }

    is_leaf = false; //нода более не лист
    rtree.vacateData(dataId_); //очистка данных
    dataId_ = null;
    rtree.incrementCounterNodes(tag().depth() + 1, QUADRANTS_NUM); //активных нод на следующем уровне стало на 4 больше
    rtree.incrementCounterLeaves(tag().depth() + 1, QUADRANTS_NUM); //листьев на следующем уровне стало на 4 больше
    rtree.incrementCounterLeaves(tag().depth(), -1); //листьев на текущем уровне стало на 1 меньше

    //обновление данных о внуках родителя
    if (parentId_ != null)
        rtree.nodeRef(tag().depth() - 1, parentId_).setGrandParency(true);

    //тэги и ссылки на детей
    ChildrenRefs crefs(childRef(Quadrant::top_left), childRef(Quadrant::top_right), childRef(Quadrant::bottom_right), childRef(Quadrant::bottom_left));
    ChildrenTags ctags(crefs(Quadrant::top_left).tag(), crefs(Quadrant::top_right).tag(), crefs(Quadrant::bottom_right).tag(), crefs(Quadrant::bottom_left).tag());

    //обновление данных о соседях
    crefs(Quadrant::top_left).setNeighbour(Neighbour::right, ctags(Quadrant::top_right)); //siblings
    crefs(Quadrant::top_left).setNeighbour(Neighbour::bottom, ctags(Quadrant::bottom_left));
    crefs(Quadrant::top_right).setNeighbour(Neighbour::left, ctags(Quadrant::top_left));
    crefs(Quadrant::top_right).setNeighbour(Neighbour::bottom, ctags(Quadrant::bottom_right));
    crefs(Quadrant::bottom_right).setNeighbour(Neighbour::top, ctags(Quadrant::top_right));
    crefs(Quadrant::bottom_right).setNeighbour(Neighbour::left, ctags(Quadrant::bottom_left));
    crefs(Quadrant::bottom_left).setNeighbour(Neighbour::top, ctags(Quadrant::top_left));
    crefs(Quadrant::bottom_left).setNeighbour(Neighbour::right, ctags(Quadrant::bottom_right));
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n)) //есть сосед
        {
            auto& rnnode = nodeRef(neighbour(n));
            if (rnnode.hasChildren()) //у соседа есть дети
            {
                rnnode.setChildrenNeighbours(opposite(n), ctags); //обновление для детей соседа
                ChildrenTags ntags = rnnode.childrenTags();
                setChildrenNeighbours(n, ntags); //для детей этой ноды
            }
            else //можно не посылать данные соседу-листу
            { 
                setNeighbour(n, neighbour(n)); //для детей этой ноды
            }
        }
    }

    //обновление данных о соседях12
    crefs(Quadrant::top_left).setNeighbour12(Neighbour12::right1, ctags(Quadrant::top_right)); //siblings
    crefs(Quadrant::top_left).setNeighbour12(Neighbour12::bottom1, ctags(Quadrant::bottom_left));
    crefs(Quadrant::top_left).setNeighbour12(Neighbour12::bottom_right, ctags(Quadrant::bottom_right));
    crefs(Quadrant::top_right).setNeighbour12(Neighbour12::left1, ctags(Quadrant::top_left));
    crefs(Quadrant::top_right).setNeighbour12(Neighbour12::bottom1, ctags(Quadrant::bottom_right));
    crefs(Quadrant::top_right).setNeighbour12(Neighbour12::bottom_left, ctags(Quadrant::bottom_left));
    crefs(Quadrant::bottom_right).setNeighbour12(Neighbour12::top1, ctags(Quadrant::top_right));
    crefs(Quadrant::bottom_right).setNeighbour12(Neighbour12::left1, ctags(Quadrant::bottom_left));
    crefs(Quadrant::bottom_right).setNeighbour12(Neighbour12::top_left, ctags(Quadrant::top_left));
    crefs(Quadrant::bottom_left).setNeighbour12(Neighbour12::top1, ctags(Quadrant::top_left));
    crefs(Quadrant::bottom_left).setNeighbour12(Neighbour12::right1, ctags(Quadrant::bottom_right));
    crefs(Quadrant::bottom_left).setNeighbour12(Neighbour12::top_right, ctags(Quadrant::top_right));
    for (auto n : Neighbours) //соседи по сторонам
    {
        if (hasNeighbour(n)) //есть сосед
        {
            auto& rnnode = nodeRef(neighbour(n));
            rnnode.setChildrenOrSelfNeighbours12(opposite(n), ctags); //обновление для соседа или его детей
            if (rnnode.hasChildren()) //у соседа есть дети
            {
                ChildrenTags ntags = rnnode.childrenTags(); //тэги детей соседа
                setChildrenOrSelfNeighbours12(n, ntags); //для детей этой ноды
            }
            else
            {
                setChildrenCommonNeighbour12(n, neighbour(n)); //для детей этой ноды
            }
        }
    }    
    for (auto n12 : Neighbours12Diagonal) //соседи по вершинам
    {
        if (hasNeighbour12(n12)) //есть сосед
        {
            auto& rnnode = nodeRef(neighbour12(n12));
            rnnode.setNeighbour12(opposite(n12), ctags(toQuadrant(n12))); //для соседа
            crefs(toQuadrant(n12)).setNeighbour12(n12, neighbour12(n12)); //для ребенка
        }
    }
    clearNeighbours12(); //удаление записей в этой ноде (необязательно?)

    //обработка ребер: создание внутренних ребер
    nodeEdgeId eid = forest.addEdge(NodeEdge(ctags(Quadrant::top_left), ctags(Quadrant::top_right), Orientation::vertical)); //создание и внесение ребра в список
    crefs(Quadrant::top_left).setEdge(toEdge(Neighbour::right), eid); //запись в ячейки
    crefs(Quadrant::top_right).setEdge(toEdge(Neighbour::left), eid);
    eid = forest.addEdge(NodeEdge(ctags(Quadrant::top_right), ctags(Quadrant::bottom_right), Orientation::horizontal));
    crefs(Quadrant::top_right).setEdge(toEdge(Neighbour::bottom), eid);
    crefs(Quadrant::bottom_right).setEdge(toEdge(Neighbour::top), eid);
    eid = forest.addEdge(NodeEdge(ctags(Quadrant::bottom_left), ctags(Quadrant::bottom_right), Orientation::vertical));
    crefs(Quadrant::bottom_left).setEdge(toEdge(Neighbour::right), eid);
    crefs(Quadrant::bottom_right).setEdge(toEdge(Neighbour::left), eid);
    eid = forest.addEdge(NodeEdge(ctags(Quadrant::top_left), ctags(Quadrant::bottom_left), Orientation::horizontal));
    crefs(Quadrant::top_left).setEdge(toEdge(Neighbour::bottom), eid);
    crefs(Quadrant::bottom_left).setEdge(toEdge(Neighbour::top), eid);
    
    //внешние ребра
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n))
        {
            auto& rnnode = nodeRef(neighbour(n));
            if (rnnode.hasChildren()) //у соседа есть дети, ребро уже разделено, нужно только обновить данные в ребрах и внести их своим детям
            {
                updateChildrenEdges(n);
            }
            else //детей нет, нужно разделить ребро
            {
                splitEdge(n);
            }
        }
    }
    //очистка записей о ребрах в данной ячейке
    for (size_t i = 0; i < EDGES_NUM; i++)
        edges[i] = null;
    return 0;
}

void TreeNode::gatherDataFromChildren() //сбор данных из детей для объединения
{
    try {
        if (dataId_ == null)
            throw std::invalid_argument("node is deleted");
        else
        {
            auto& rd = dataRef();
            rd.clear(); //обнуление
            //осреднение консервативных величин из детей
            for (auto q : Quadrants)
            {
                CellData cd = childRef(q).data();
                rd.add(cd);
            }
            rd.divide(static_cast<double>(QUADRANTS_NUM));
            if (config.coord_type == CoordType::axisymmetric) //внесение r в cellData для осесимметричных координат
                rd.setY(box().center().y);
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "TreeNode.gatherDataFromChildren() error: " << e.what() << ", tag: " << tag() << endl;
    }
    return;
}

const int INT_ERROR_CODE_CANT_COARSEN_DELETED = -1;
const int INT_ERROR_CODE_CANT_COARSEN_LEAF = -2;
const int INT_ERROR_CODE_CANT_COARSEN_GRANDPARENT = -3;
const int INT_ERROR_CODE_CANT_COARSEN_UNBALANCED = -4;
int TreeNode::coarsen() //склейка ячейки
{
    //удаленные ячейки, листы или дедушек нельзя склеивать
    if (is_deleted)
        return INT_ERROR_CODE_CANT_COARSEN_DELETED;
    if (is_leaf)
        return INT_ERROR_CODE_CANT_COARSEN_LEAF;
    if (has_grandchildren)
        return INT_ERROR_CODE_CANT_COARSEN_GRANDPARENT;

    //проверка на баланc
    if (!forest.treeRef(tag().tree()).isGhost()) //в ghost-деревьях баланс не нужен
    {
        for (auto n : Neighbours)
        {
            if (hasNeighbour(n) && nodeRef(neighbour(n)).hasGrandChildren(opposite(n))) //у соседа есть внуки с ближней стороны
            {
                return INT_ERROR_CODE_CANT_COARSEN_UNBALANCED;
            }
        }
    }

    auto& rtree = forest.treeRef(tag().tree());
    //получение номера ячейки данных для текущей ноды и сбор данных в неё
    dataId_ = rtree.getVacantDataId();
    gatherDataFromChildren();

    //обновление данных в соседних ячейках
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n))
        {
            nodeRef(neighbour(n)).setNeighbour(opposite(n), tag()); //для соседа и его детей
        }
    }

    //ссылки на детей
    ChildrenRefs crefs(childRef(Quadrant::top_left), childRef(Quadrant::top_right), childRef(Quadrant::bottom_right), childRef(Quadrant::bottom_left));

    //обновление данных о соседях12
    for (auto n : Neighbours) //соседи по сторонам
    {
        if (hasNeighbour(n)) //есть сосед
        {
            auto& rnnode = nodeRef(neighbour(n));
            rnnode.setChildrenOrSelfNeighbour12(opposite(n), tag()); //для соседа или его детей
            if (rnnode.hasChildren()) //у соседа есть дети
            {
                ChildrenTags ntags = rnnode.childrenTags(); //тэги детей соседа
                setNeighbours12(n, ntags); //два соседа12 для этой ноды
            }
            else
            {
                setNeighbour12(toNeighbour12(n), rnnode.tag()); //один сосед12 для этой ноды
            }
        }
    }
    for (auto n12 : Neighbours12Diagonal) //соседи по вершинам
    {
        if (crefs(toQuadrant(n12)).hasNeighbour12(n12)) //у ребенка есть сосед по диагонали
        {
            auto& rnnode = nodeRef(crefs(toQuadrant(n12)).neighbour12(n12));
            rnnode.setNeighbour12(opposite(n12), tag()); //для соседа
            setNeighbour12(n12, rnnode.tag()); //для этой ноды
        }
    }

    //обработка ребер
    forest.removeEdge(crefs(Quadrant::top_left).edge(Edge::right1)); //удаление внутренних ребер
    forest.removeEdge(crefs(Quadrant::top_left).edge(Edge::bottom1));
    forest.removeEdge(crefs(Quadrant::bottom_right).edge(Edge::top1));
    forest.removeEdge(crefs(Quadrant::bottom_right).edge(Edge::left1));
    for (auto n : Neighbours) //боковые ребра
    {
        if (hasNeighbour(n))
        {
            auto& rnnode = nodeRef(neighbour(n));
            if (rnnode.hasChildren()) //у соседа есть дети - только обновляем данные о ребрах и в ребрах 
            {
                gatherEdgesFromChildren(n);
            }
            else //нет детей - склеиваем и обновляем записи о ребре и в ребре
            {
                joinEdge(n);
            }
        }
    }

    //удаление детей
    rtree.vacateNodeGroup(tag().depth()+1, childrenId_); //пометка группы ячеек для нод свободными
    for (auto q : Quadrants)
    {
        auto& rchild = childRef(q);
        rchild.markDeleted();
        rtree.vacateData(rchild.dataId()); //пометка ячеек в массиве tree.data свободными
    }
    rtree.incrementCounterNodes(tag().depth() + 1, -QUADRANTS_NUM); //активных нод на следующем уровне стало на 4 меньше
    rtree.incrementCounterLeaves(tag().depth() + 1, -QUADRANTS_NUM); //листьев на следующем уровне стало на 4 меньше
    rtree.incrementCounterLeaves(tag().depth(), 1); //листьев на текущем уровне стало на 1 больше
        
    //ячейка стала листом
    is_leaf = true;
    childrenId_ = null;

    //обновление данных о внуках родителя
    if (parentId_ != null)
        nodeRef(NodeTag(tag().tree(), tag().depth() - 1, parentId_)).updateGrandParency();

    //удаление слоя дерева при необходимости
    rtree.deleteLevelIfEmpty(tag().depth());

    //выравнивание размера соседних ghost'ов
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n))
        {
            if (!forest.treeRef(tag().tree()).isGhost() && forest.treeRef(neighbour(n).tree()).isGhost()) //есть ghost-сосед
            {
                auto& rnode = nodeRef(neighbour(n));
                if (rnode.hasChildren())
                    rnode.coarsen();
            }
        }
    }
    return 0;
}

//TODO: сделать корректную версию для осесимметричных координат
void TreeNode::updateEigenObjects() //создание или обновление Eigen матриц и т.д. после изменения сетки
{
    static int variables_num = POLY_COEFF_NUM; //число неизвестных (= число столбцов матрицы)
    const int equations_num = neighbours12Num(); //число уравнений (= число строк матрицы)
    coeffs = Eigen::VectorXd(variables_num);
    J = Eigen::MatrixXd(equations_num, variables_num);
    rs = Eigen::VectorXd(equations_num);
    decomp = Eigen::HouseholderQR<Eigen::MatrixXd>(equations_num, variables_num);

    const double x = box().center().x; //координаты центра данной ячейки
    const double y = box().center().y;
    const double h = box().size();
    int row = 0;
    for (auto& rn12tag : neighbours12)
    {
        if (!rn12tag.isNull())
        {
            auto& nnode = nodeRef(rn12tag);
            const double dx = nnode.box().center().x - x; //компоненты отрезка до центра соседа
            const double dy = nnode.box().center().y - y;
            const double nh = nnode.box().size(); //длина стороны соседа
            int col = 0;
            J(row, col++) = dx;
            J(row, col++) = dy;
            J(row, col++) = 0.5 * (dx * dx + 1.0 / 12.0 * (nh * nh - h * h));
            J(row, col++) = dx * dy;
            J(row, col++) = 0.5 * (dy * dy + 1.0 / 12.0 * (nh * nh - h * h));
            row++;
        }
    }
    decomp.compute(J);

    //линейные функции (TODO: переписать без массивов - посчитать, сколько точек, создать объекты Eigen и пробежаться по точкам)
    static const int PTS_MAX = 5; //максимальное число точек, что могут попасть в подшаблон
    Point ncenters[PTS_MAX]{};
    CellData nds[PTS_MAX];
    coeffsl = Eigen::VectorXd(2);
    for (auto q : Quadrants)
    {
        int pts = 0; //число ячеек, вошедших в подшлаблон
        for (auto& rn12tag : neighbours12)
        {
            if (!rn12tag.isNull())
            {
                auto& rnnode = nodeRef(rn12tag);
                if (isNodeInSubstencil(q, rnnode))
                {
                    ncenters[pts] = rnnode.box().center();
                    pts++;
                }
            }
        }
        const int nq = static_cast<int>(q);
        Jl[nq] = Eigen::MatrixXd(pts, 2);
        rsl[nq] = Eigen::VectorXd(pts);
        decompl[nq] = Eigen::HouseholderQR<Eigen::MatrixXd>(pts, 2);
        for (int p = 0; p < pts; p++)
        {
            double dx = ncenters[p].x - x;
            double dy = ncenters[p].y - y;
            Jl[nq](p, 0) = dx;
            Jl[nq](p, 1) = dy;
        }
        decompl[nq].compute(Jl[nq]);
    }
    return;
}

static const double ALPHA0 = 0.5; //коэффициент для параболоида P0
static const double ALPHA1 = 0.125; //для остальных (линейных) 4 функций
void TreeNode::calcPolynomialCWENO(rkStep rk) //вычисление коэффициентов 2D CWENO полинома
{
    static const int variables_num = POLY_COEFF_NUM; //число неизвестных (= число столбцов матрицы)
    if (forest.treeRef(tag().tree()).isGhost()) //нода в ghost-дереве
    {
        for (auto eq : Equations)
            for (int i = 0; i < POLY_COEFF_NUM; i++)
                polyCoeffs[static_cast<int>(eq)][i] = 0.0;
        return;
    }

    const int equations_num = neighbours12Num(); //число уравнений (= число строк матрицы)
    auto& rdata = dataRef();

    double Popt_coeffs[POLY_COEFF_NUM]{}; //коэффициенты для оптимального (осциллирующего) полинома
    Quadrant cq = static_cast<Quadrant>(tag().id() % QUADRANTS_NUM); //квадрант текущей ноды среди прочих siblings
    for (auto eq : Equations)
    {
        //оптимальный полином
        int row = 0;
        for (auto& rn12tag : neighbours12)
        {
            if (!rn12tag.isNull())
            {
                auto& rnnode = nodeRef(rn12tag);
                auto& rndata = rnnode.dataRef();
                rs(row) = rndata.Qref(rk)(eq) - rdata.Qref(rk)(eq);
                row++;
            }
        }
        coeffs = decomp.solve(rs);
        for (int c = 0; c < POLY_COEFF_NUM; c++)
            Popt_coeffs[c] = coeffs(c);

        //линейные функции (TODO: переписать без массивов - посчитать, сколько точек, создать объекты Eigen и пробежаться по точкам)
        static const int PTS_MAX = 5; //максимальное число точек, что могут попасть в подшаблон
        CellData nds[PTS_MAX];
        double linear_coeffs[QUADRANTS_NUM][2];
        for (auto q : Quadrants)
        {
            int pts = 0; //число ячеек, вошедших в подшлаблон
            for (auto& rn12tag : neighbours12)
            {
                if (!rn12tag.isNull())
                {
                    auto& rnnode = nodeRef(rn12tag);
                    if (isNodeInSubstencil(q, rnnode))
                    {
                        nds[pts] = rnnode.data(); //копия данных
                        pts++;
                    }
                }
            }
            const int nq = static_cast<int>(q);
            for (int p = 0; p < pts; p++)
            {
                rsl[nq](p) = nds[p].Qref(rk)(eq) - rdata.Qref(rk)(eq);
            }
            coeffsl = decompl[nq].solve(rsl[nq]);
            linear_coeffs[nq][0] = coeffsl(0);
            linear_coeffs[nq][1] = coeffsl(1);
        } //quadrants loop

        //параболоид P0: Popt = ALPHA0 P0 + sum(ALPHA1 * linear funcs)
        double P0_coeffs[POLY_COEFF_NUM]{};
        for (int c = 0; c < POLY_COEFF_NUM; c++)
            P0_coeffs[c] = Popt_coeffs[c] / ALPHA0;
        for (auto q : Quadrants)
        {
            const int nq = static_cast<int>(q);
            P0_coeffs[0] -= ALPHA1 * linear_coeffs[nq][0] / ALPHA0;
            P0_coeffs[1] -= ALPHA1 * linear_coeffs[nq][1] / ALPHA0;
        }

        //WENO3 веса
        double h = box().size(); //длина стороны ячейки
        double beta[QUADRANTS_NUM + 1]{}; //smoothness indicators
        double omega[QUADRANTS_NUM + 1]{}; //weights
        double omega_sum = 0.0;
        for (auto q : Quadrants)
        {
            const int nq = static_cast<int>(q);
            beta[nq] = h * h * (linear_coeffs[nq][0] * linear_coeffs[nq][0] + linear_coeffs[nq][1] * linear_coeffs[nq][1]);
            omega[nq] = ALPHA1 / (h + beta[nq]) / (h + beta[nq]); //h вместо epsilon (Semplice2015)
            omega_sum += omega[nq];
        }
        beta[QUADRANTS_NUM] = h * h * (P0_coeffs[0] * P0_coeffs[0] + P0_coeffs[1] * P0_coeffs[1]
            + h * h * (13.0 / 12.0 * P0_coeffs[2] * P0_coeffs[2] + 7.0 / 6.0 * P0_coeffs[3] * P0_coeffs[3] + 13.0 / 12.0 * P0_coeffs[4] * P0_coeffs[4]));
        omega[QUADRANTS_NUM] = ALPHA0 / (h + beta[QUADRANTS_NUM]) / (h + beta[QUADRANTS_NUM]);
        omega_sum += omega[QUADRANTS_NUM];
        double alpha[QUADRANTS_NUM + 1]{}; //normalized weights
        for (auto q : Quadrants)
            alpha[static_cast<int>(q)] = omega[static_cast<int>(q)] / omega_sum;
        alpha[QUADRANTS_NUM] = omega[QUADRANTS_NUM] / omega_sum;

        //финальная реконструкция
        const int neq = static_cast<int>(eq);
        for (int c = 0; c < POLY_COEFF_NUM; c++)
            polyCoeffs[neq][c] = alpha[QUADRANTS_NUM] * P0_coeffs[c];
        for (auto q : Quadrants)
        {
            const int nq = static_cast<int>(q);
            polyCoeffs[neq][0] += alpha[nq] * linear_coeffs[nq][0];
            polyCoeffs[neq][1] += alpha[nq] * linear_coeffs[nq][1];
        }
    } //eq loop
    return;
}

//inspectors -----------------------
bool TreeNode::isDeleted() const { return is_deleted; }
bool TreeNode::isLeaf() const { return is_leaf; }
bool TreeNode::hasChildren() const { return childrenId_ != null; } //есть ли дети
bool TreeNode::hasGrandChildren() const { return has_grandchildren; } //есть ли внуки
bool TreeNode::hasGrandChildren(Neighbour n) const //есть ли внуки с определенной стороны
{
    if (!hasChildren())
        return false;
    switch (n)
    {
    case Neighbour::top:
        if (childRef(Quadrant::top_left).hasChildren() ||
            childRef(Quadrant::top_right).hasChildren())
        {
            return true;
        }
        break;
    case Neighbour::right:
        if (childRef(Quadrant::top_right).hasChildren() ||
            childRef(Quadrant::bottom_right).hasChildren())
        {
            return true;
        }
        break;
    case Neighbour::bottom:
        if (childRef(Quadrant::bottom_right).hasChildren() ||
            childRef(Quadrant::bottom_left).hasChildren())
        {
            return true;
        }
        break;
    case Neighbour::left:
        if (childRef(Quadrant::bottom_left).hasChildren() ||
            childRef(Quadrant::top_left).hasChildren())
        {
            return true;
        }
        break;
    default:
        cout << "TreeNode.hasGrandChildren(Neighbour) error: invalid quadrant" << endl;
        break;
    }
    return false;
}
bool TreeNode::hasNeighbour(Neighbour n) const //есть ли сосед по направлению
{
    if (!neighbours[static_cast<int>(n)].isNull())
        return true;
    return false;
}
bool TreeNode::hasNeighbour12(Neighbour12 n12) const
{
    if (!neighbours12[static_cast<int>(n12)].isNull())
        return true;
    return false;
}
bool TreeNode::hasEdge(Edge etype) const //есть ли ребро
{
    if (edges[static_cast<int>(etype)] != null)
        return true;
    return false;
}
int TreeNode::neighbours12Num() const //число соседей12 TODO: кешировать(?)
{
    int ret = 0;
    for (auto n12 : Neighbours12)
        if (hasNeighbour12(n12))
            ret++;
    return ret;
}
bool TreeNode::isNodeInSubstencil(Quadrant q, const TreeNode& rnnode) const //попадает ли ячейка в подшаблон для линейной функции
{
    switch (q)
    {
    case Quadrant::top_left:
        if (rnnode.box().center().x - 0.5 * rnnode.box().size() - DOUBLE_EPS12 < box().center().x &&
            rnnode.box().center().y + 0.5 * rnnode.box().size() + DOUBLE_EPS12 > box().center().y)
            return true;
        break;
    case Quadrant::top_right:
        if (rnnode.box().center().x + 0.5 * rnnode.box().size() + DOUBLE_EPS12 > box().center().x &&
            rnnode.box().center().y + 0.5 * rnnode.box().size() + DOUBLE_EPS12 > box().center().y)
            return true;
        break;
    case Quadrant::bottom_right:
        if (rnnode.box().center().x + 0.5 * rnnode.box().size() + DOUBLE_EPS12 > box().center().x &&
            rnnode.box().center().y - 0.5 * rnnode.box().size() - DOUBLE_EPS12 < box().center().y)
            return true;
        break;
    case Quadrant::bottom_left:
        if (rnnode.box().center().x - 0.5 * rnnode.box().size() - DOUBLE_EPS12 < box().center().x &&
            rnnode.box().center().y - 0.5 * rnnode.box().size() - DOUBLE_EPS12 < box().center().y)
            return true;
        break;
    }
    return false;
}
//uncategorized -----------------------
TreeNode& TreeNode::nodeRef(const NodeTag& tag) //ссылка на ноду по тэгу
{
    try {
        TreeNode& rnode = forest.treeRef(tag.tree()).nodeRef(tag.depth(), tag.id());
        if(rnode.isDeleted())
        {
            throw std::invalid_argument("node is deleted");
        }
        return rnode;
    }
    catch (const std::invalid_argument& e)
    {
        cout << "TreeNode.nodeRef() error: " << e.what() << ", tag: " << tag << endl;
    }
    return nodeRef(NodeTag(0, 1, 0)); //to suppress warning
}

double TreeNode::magGradRho() const //примерный градиент плотности
{
    if (!is_leaf || is_deleted)
        return 0.0;

    const CellData& rd = dataRefConst();
    if (rd.rho() < DOUBLE_EPS12)
        return 0.0;
    double drhosum = 0.0;
    int n12num = 0;
    for (auto& rn12 : neighbours12)
    {
        if (!rn12.isNull() && !forest.treeRef(rn12.tree()).isGhost())
        {
            auto& rnnode = nodeRef(rn12);
            const auto& rndata = rnnode.dataRefConst(); //можно по ссылке, т.к. в neighbour12 лежат только листья
            drhosum += fabs(rd.rho() - rndata.rho()) / distance(box().center(), rnnode.box().center());
            n12num++;
        }
    }
    return drhosum / n12num / rd.rho();
}

ConservativeVector TreeNode::evalPolynomialAt(Point p, rkStep rk) //реконструированное полиномом CWENO значение TODO: сделать const
{
    ConservativeVector ret;
    ConservativeVector& rQ = dataRef().Qref(rk);
    double dx = p.x - box().center().x;
    double dy = p.y - box().center().y;
    double h = box().size();
    for (auto eq : Equations)
    {
        ret.set(eq, rQ(eq) + polyCoeff(eq, 0) * dx + polyCoeff(eq, 1) * dy
            + 0.5 * polyCoeff(eq, 2) * (dx * dx - 1.0 / 12.0 * h * h) + 0.5 * polyCoeff(eq, 3) * (dy * dy - 1.0 / 12.0 * h * h)
            + polyCoeff(eq, 4) * dx * dy);
    }
    return ret;
}

//output -----------------------
std::string TreeNode::dump() const //дамп ноды в строку
{
    std::stringstream buffer;
    buffer << "(" << std::setw(3) << (tag_.tree() == null ? " - " : std::to_string(tag_.tree())) << ", " << tag_.depth() << ", " << std::setw(3) << (tag_.id() == null ? " - " : std::to_string(tag_.id())) << ", [" << std::setw(3) << (parentId_ == null ? " - " : std::to_string(parentId_)) << "]): ";
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n))
            buffer << static_cast<int>(n) << ": " << neighbour(n) << " | ";
        else
            buffer << "               | ";
    }
    if (is_deleted)
        buffer << " deleted";
    buffer << endl;

    /*buffer << "                     [";
    for (auto n : Neighbours)
    {
        if (neighbours[n].tree_id != null && neighbours[n].id != null)
        {
            auto& nnode = nodeRef(neighbours[n]);
            Neighbour dir = (Neighbour)((int)(n + 2) % 4); //противоположное направление
            buffer << Dirs[dir] << ": " << std::setw(3) << nnode.neighbours[dir].tree_id << ", " << nnode.neighbours[dir].depth << ", " << std::setw(3) << nnode.neighbours[dir].id << " | ";
        }
        else
            buffer << "               | ";
    }
    buffer << "]" << endl;*/

    return buffer.str();
}
static const double Dx[] = { 0, 1, 0, -1 }; //для вычисления векторов
static const double Dy[] = { 1, 0, -1, 0 };
std::string TreeNode::dumpNeighbourVector(Neighbour n) const //дамп соседа в виде вектора
{
    if (!hasNeighbour(n))
        return "(no such neighbour)";

    std::stringstream buffer;
    double x = box().center().x + 0.8 * Dx[static_cast<int>(n)] * box().size() / 2.0; //стартовая точка вектора
    double y = box().center().y + 0.8 * Dy[static_cast<int>(n)] * box().size() / 2.0;
    auto& rnnode = nodeRef(neighbour(n));
    double dx = rnnode.box().center().x - x; //вектор
    double dy = rnnode.box().center().y - y;

    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;
    return buffer.str();
}
static const double Dx12[] = { -0.5, 0.5, 1,   1,    1,  1, 0.5, -0.5, -1,   -1,  -1, -1 };
static const double Dy12[] = {    1,   1, 1, 0.5, -0.5, -1,  -1,   -1, -1, -0.5, 0.5,  1 };
std::string TreeNode::dumpNeighbour12Vector(Neighbour12 n12) const
{
    if (!hasNeighbour12(n12))
        return "(no such neighbour12)";

    std::stringstream buffer;
    double x = box().center().x + 0.9 * Dx12[static_cast<int>(n12)] * box().size() / 2.0; //стартовая точка вектора
    double y = box().center().y + 0.9 * Dy12[static_cast<int>(n12)] * box().size() / 2.0;
    TreeNode& rnnode = nodeRef(neighbour12(n12));
    double dx = 0.85 * (rnnode.box().center().x - x); //вектор
    double dy = 0.85 * (rnnode.box().center().y - y);

    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;
    return buffer.str();
}
static const double Dxe[] = { -0.5, 0.5,   1,    1, 0.5, -0.5,   -1,  -1 }; //для вычисления стартовых точек
static const double Dye[] = {    1,   1, 0.5, -0.5,  -1,   -1, -0.5, 0.5 };
std::string TreeNode::dumpEdgeVector(Edge etype) const
{
    if (!hasEdge(etype))
        return "(no such edge)";

    std::stringstream buffer;
    NodeEdge& redge = forest.edgeRef(edge(etype));
    double x = box().center().x + 0.5 * Dxe[static_cast<int>(etype)] * box().size() / 2.0; //стартовая точка вектора
    double y = box().center().y + 0.5 * Dye[static_cast<int>(etype)] * box().size() / 2.0;
    Point ecenter = redge.middle();
    double dx = 0.85 * (ecenter.x - x); //вектор
    double dy = 0.85 * (ecenter.y - y);

    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;
    return buffer.str();
}
//------------------------------------------------------------------- ChildrenRefs -------------------------------------------------------------------
TreeNode& ChildrenRefs::operator()(Quadrant q)
{
    switch (q)
    {
    case Quadrant::top_left: return rTL;
    case Quadrant::top_right: return rTR;
    case Quadrant::bottom_right: return rBR;
    case Quadrant::bottom_left: return rBL;
    default:
        cout << "ChildrenRefs() error: invalid quadrant" << endl;
        return rTL; //to suppress warning
        break;
    }
}
