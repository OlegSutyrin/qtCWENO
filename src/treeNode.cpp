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
const NodeTag& ChildrenTags::operator()(Quadrant q) const { return tags[static_cast<int>(q)]; } //() operator: ��� �� ���������

//------------------------------------------------------------------- TreeNode -------------------------------------------------------------------
const NodeTag TreeNode::tag() const { return tag_; } //��������� ���� 
cellDataId TreeNode::dataId() const { return dataId_; } //��������� dataId
CellBox TreeNode::box() const {	return box_; } //��������� box'a
//QuadTree& TreeNode::treeRef() { return forest.treeRef(tag().tree()); } // ������ �� ������, ���������� ����
//const QuadTree& TreeNode::treeRefConst() const { return forest.treeRefConst(tag().tree()); } //(const) ������ �� ������, ���������� ����
const NodeTag TreeNode::neighbour(Neighbour n) const { return neighbours[static_cast<int>(n)]; }
const NodeTag TreeNode::neighbour12(Neighbour12 n12) const { return neighbours12[static_cast<int>(n12)]; }

CellData& TreeNode::dataRef() //������ �� ������
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
CellData TreeNode::data() const //����� ������
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
        else //����� ������� ������ � �����
        {
            CellData d; //���������������� ������
            for (auto q : Quadrants) //����������� �������������� ���� ������ � �����
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

TreeNode& TreeNode::childRef(Quadrant q) //������ �� ������� �� ���������
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
const TreeNode& TreeNode::childRef(Quadrant q) const //const ������
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
const ChildrenTags TreeNode::childrenTags() const //���� ���� �����
{
    return ChildrenTags(childRef(Quadrant::top_left).tag(), childRef(Quadrant::top_right).tag(), childRef(Quadrant::bottom_right).tag(), childRef(Quadrant::bottom_left).tag());
}

TreeNode& TreeNode::getChildOrSelfByCoords(Point p) //������ �� ������� �� �����������
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
void TreeNode::setTag(const NodeTag& t) { tag_ = t; } //������� ����
void TreeNode::markDeleted() { is_deleted = true; } //������� ���������
void TreeNode::setDataId(cellDataId id) { dataId_ = id; } //������� dataId
void TreeNode::setBox(const CellBox& b) { box_ = b; } //������� box'�
void TreeNode::setNeighbour(Neighbour n, NodeTag ntag) //������� ������ ������ � ����� � ������ �������
{ 
    neighbours[static_cast<int>(n)] = ntag; //������ ����
    if (hasChildren()) //�����
    {
        switch (n) 
        {
        case Neighbour::top: //������ ��� ������ ������
            childRef(Quadrant::top_left).setNeighbour(Neighbour::top, ntag);
            childRef(Quadrant::top_right).setNeighbour(Neighbour::top, ntag);
            break;
        case Neighbour::right: //������
            childRef(Quadrant::top_right).setNeighbour(Neighbour::right, ntag);
            childRef(Quadrant::bottom_right).setNeighbour(Neighbour::right, ntag);
            break;
        case Neighbour::bottom: //�����
            childRef(Quadrant::bottom_right).setNeighbour(Neighbour::bottom, ntag);
            childRef(Quadrant::bottom_left).setNeighbour(Neighbour::bottom, ntag);
            break;
        case Neighbour::left: //�����
            childRef(Quadrant::top_left).setNeighbour(Neighbour::left, ntag);
            childRef(Quadrant::bottom_left).setNeighbour(Neighbour::left, ntag);
            break;
        }

    }
}
void TreeNode::setChildrenNeighbours(Neighbour n, ChildrenTags tags) //�������� ������ � ������� ��� �����
{
    if (!hasChildren()) //���� ����� �� ������ ���������� ��� �������
    {
        cout << "TreeNode.setChildrenNeighbours() error: no children, tag: " << tag() << endl;
        return;
    }
    switch (n)
    {
    case Neighbour::top: //������ ������ �� ������ ������
        childRef(Quadrant::top_left).setNeighbour(Neighbour::top, tags(Quadrant::bottom_left));
        childRef(Quadrant::top_right).setNeighbour(Neighbour::top, tags(Quadrant::bottom_right));
        break;
    case Neighbour::right: //������
        childRef(Quadrant::top_right).setNeighbour(Neighbour::right, tags(Quadrant::top_left));
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::right, tags(Quadrant::bottom_left));
        break;
    case Neighbour::bottom: //�����
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::bottom, tags(Quadrant::top_right));
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::bottom, tags(Quadrant::top_left));
        break;
    case Neighbour::left: //�����
        childRef(Quadrant::top_left).setNeighbour(Neighbour::left, tags(Quadrant::top_right));
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::left, tags(Quadrant::bottom_right));
        break;
    }
}

void TreeNode::setNeighbour12(Neighbour12 n12, NodeTag ntag) { neighbours12[static_cast<int>(n12)] = ntag; }

void TreeNode::setChildrenOrSelfNeighbours12(Neighbour n, ChildrenTags tags) //�������� ������ � �������12 ��� ���� ��� �����
{
    if (!hasChildren()) 
    {
        switch (n)
        {
        case Neighbour::top: //������ ������ �� ������ ������
            setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_left));
            setNeighbour12(Neighbour12::top2, tags(Quadrant::bottom_right));
            break;
        case Neighbour::right: //������
            setNeighbour12(Neighbour12::right1, tags(Quadrant::top_left));
            setNeighbour12(Neighbour12::right2, tags(Quadrant::bottom_left));
            break;
        case Neighbour::bottom: //�����
            setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_right));
            setNeighbour12(Neighbour12::bottom2, tags(Quadrant::top_left));
            break;
        case Neighbour::left: //�����
            setNeighbour12(Neighbour12::left1, tags(Quadrant::bottom_right));
            setNeighbour12(Neighbour12::left2, tags(Quadrant::top_right));
            break;
        }
        return;
    }
    switch (n)
    {
    case Neighbour::top: //������ ������ �� ������ ������
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_left)); //������������ ����� �� �������
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top_right, tags(Quadrant::bottom_right)); //�� ���������
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top1, tags(Quadrant::bottom_right));
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top_left, tags(Quadrant::bottom_left));
        break;
    case Neighbour::right: //������
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::right1, tags(Quadrant::top_left));
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::bottom_right, tags(Quadrant::bottom_left));
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::right1, tags(Quadrant::bottom_left));
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::top_right, tags(Quadrant::top_left));
        break;
    case Neighbour::bottom: //�����
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_right));
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom_left, tags(Quadrant::top_left));
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom1, tags(Quadrant::top_left));
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom_right, tags(Quadrant::top_right));
        break;
    case Neighbour::left: //�����
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::left1, tags(Quadrant::top_right));
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::bottom_left, tags(Quadrant::bottom_right));
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::left1, tags(Quadrant::bottom_right));
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::top_left, tags(Quadrant::top_right));
        break;
    }
}
void TreeNode::setChildrenCommonNeighbour12(Neighbour n, NodeTag ntag) //�������� ������ �� ����� ������12 ��� �����
{
    if (!hasChildren()) //���� ����� �� ������ ���������� ��� �������
    {
        cout << "TreeNode.setChildrenCommonNeighbour12() error: no children, tag: " << tag() << endl;
        return;
    }
    switch (n)
    {
    case Neighbour::top: //������ ������ �� ������ ������
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::top1, ntag); //������������ ����� �� �������
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::top1, ntag);
        break;
    case Neighbour::right: //������
        childRef(Quadrant::top_right).setNeighbour12(Neighbour12::right1, ntag);
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::right1, ntag);
        break;
    case Neighbour::bottom: //�����
        childRef(Quadrant::bottom_right).setNeighbour12(Neighbour12::bottom1, ntag);
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::bottom1, ntag);
        break;
    case Neighbour::left: //�����
        childRef(Quadrant::top_left).setNeighbour12(Neighbour12::left1, ntag);
        childRef(Quadrant::bottom_left).setNeighbour12(Neighbour12::left1, ntag);
        break;
    }
}

void TreeNode::setData(const CellData& data) //������ ������
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
void TreeNode::setGrandParency(bool status) { has_grandchildren = status; }  //������� � ������� ������

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

/*void TreeNode::setChildrenCommonNeighbour(Neighbour n, NodeTag ntag) //�������� ������ �� ����� ������ ��� �����
{
    if (!hasChildren()) //���� ����� �� ������ ���������� ��� �������
    {
        cout << "TreeNode.setChildrenCommonNeighbour() error: no children, tag: " << tag() << endl;
        return;
    }
    switch (n)
    {
    case Neighbour::top: //������ ������ �� ������ ������
        childRef(Quadrant::top_left).setNeighbour(Neighbour::top, ntag);
        childRef(Quadrant::top_right).setNeighbour(Neighbour::top, ntag);
        break;
    case Neighbour::right: //�� ������ ������
        childRef(Quadrant::top_right).setNeighbour(Neighbour::right, ntag);
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::right, ntag);
        break;
    case Neighbour::bottom: //�����
        childRef(Quadrant::bottom_right).setNeighbour(Neighbour::bottom, ntag);
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::bottom, ntag);
        break;
    case Neighbour::left: //�����
        childRef(Quadrant::top_left).setNeighbour(Neighbour::left, ntag);
        childRef(Quadrant::bottom_left).setNeighbour(Neighbour::left, ntag);
        break;
    }
}*/

int TreeNode::markToRefine() //������� ������ � ���������
{
    //��������� ������ ������ ���������
    if (is_deleted)
        return INT_ERROR_CODE_CANT_REFINE_DELETED;

    //������ ��� ������ �������
    if (tag().depth() >= config.max_depth)
        return INT_ERROR_CODE_CANT_REFINE_DEEPEST_LEVEL;

    //������ ��� ���������
    if (!is_leaf)
        return INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED;

    forest.addNodeToRefine(tag());

    //������������
    if(!forest.treeRef(tag().tree()).isGhost()) //��� ghost-�������� ������������ �� �����
    {
        for (auto& rntag : neighbours)
        {
            if (!rntag.isNull())
            {
                //(?)�������� �� ghost-������ ������ ��� �� �����, �.�. ��� ������ ������������� ����
                if (rntag.depth() < tag().depth()) //����� �� 2+ ������ ������ ������� ����� ������ ����
                    nodeRef(rntag).markToRefine();
                //������������ ������� �������� ghost'��
                if (forest.treeRef(rntag.tree()).isGhost() && rntag.depth() <= tag().depth()) //ghost-����� ������� ������� ����� ������ ����
                    nodeRef(rntag).markToRefine();
            }
        }
    }
    return 0;
}

int TreeNode::refine() //��������� ������
{
    //��������� ������ ������ ���������
    if (is_deleted)
        return INT_ERROR_CODE_CANT_REFINE_DELETED;

    //������ ��� ������ �������
    if (tag().depth() >= config.max_depth)
        return INT_ERROR_CODE_CANT_REFINE_DEEPEST_LEVEL;

    //������ ��� ���������
    if (!is_leaf)
        return INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED;

    auto& rtree = forest.treeRef(tag().tree());
    CellData cdata = this->data(); //�� ��������, ����� ��������� ������������

    //��������� ������ ��� ����� ������� (���� �����)
    if (rtree.depth() < tag().depth() + 1)
        rtree.initNewLevel();

    //������ ���-�����
    childrenId_ = rtree.getVacantNodeId(tag().depth() + 1); //��������� id ���������� ����� ���� � ��� ������ � ��� ����
    for (auto q : Quadrants)
    {
        auto& rchild = rtree.nodeRef(tag().depth() + 1, childrenId_ + static_cast<int>(q)); //TODO: ������� ��������� ����� ����� ����� ������������� � ������� �������
        auto cbox = box().quarterBox(q);
        //������ ����� �����������, ��� ����������� ���� ���������������� ���������� �� ���������
        rchild = TreeNode(
            NodeTag(tag().tree(), tag().depth() + 1, childrenId_ + static_cast<int>(q)), 
            tag().id(), 
            rtree.getVacantDataId(), 
            cbox);

        //�������������� ���������� ������ �� �����
        if (config.coord_type == CoordType::axisymmetric) //���� r � cellData'� ����� ��� �������������� ���������
            cdata.setY(cbox.center().y);
        rchild.setData(cdata);
    }

    is_leaf = false; //���� ����� �� ����
    rtree.vacateData(dataId_);
    dataId_ = null; //������ ��� ���� ������ ����� ���, ���� ������ ��� �����
    rtree.incrementNodesCounter(tag().depth() + 1, QUADRANTS_NUM); //�������� ��� �� ��������� ������ ����� �� 4 ������
    rtree.incrementLeavesCounter(tag().depth() + 1, QUADRANTS_NUM); //������� �� ��������� ������ ����� �� 4 ������
    rtree.incrementLeavesCounter(tag().depth(), -1); //������� �� ������� ������ ����� �� 1 ������

    //���������� ������ � ������ ��������
    if (parentId_ != null)
        rtree.nodeRef(tag().depth() - 1, parentId_).setGrandParency(true);

    //���� � ������ �� �����
    ChildrenRefs crefs(childRef(Quadrant::top_left), childRef(Quadrant::top_right), childRef(Quadrant::bottom_right), childRef(Quadrant::bottom_left));
    ChildrenTags ctags(crefs(Quadrant::top_left).tag(), crefs(Quadrant::top_right).tag(), crefs(Quadrant::bottom_right).tag(), crefs(Quadrant::bottom_left).tag());

    //���������� ������ � �������
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
        if (hasNeighbour(n)) //���� �����
        {
            auto& rnnode = nodeRef(neighbour(n));
            if (rnnode.hasChildren()) //� ������ ���� ����
            {
                rnnode.setChildrenNeighbours(opposite(n), ctags); //���������� ��� ����� ������
                ChildrenTags ntags = rnnode.childrenTags();
                setChildrenNeighbours(n, ntags); //��� ����� ���� ����
            }
            else //����� �� �������� ������ ������-�����
            { 
                setNeighbour(n, neighbour(n)); //��� ����� ���� ����
            }
        }
    }

    //���������� ������ � �������12
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
    for (auto n : Neighbours) //������ �� ��������
    {
        if (hasNeighbour(n)) //���� �����
        {
            auto& rnnode = nodeRef(neighbour(n));
            rnnode.setChildrenOrSelfNeighbours12(opposite(n), ctags); //���������� ��� ������ ��� ��� �����
            if (rnnode.hasChildren()) //� ������ ���� ����
            {
                ChildrenTags ntags = rnnode.childrenTags();
                setChildrenOrSelfNeighbours12(n, ntags); //��� ����� ���� ����
            }
            else
            {
                setChildrenCommonNeighbour12(n, neighbour(n)); //��� ����� ���� ����
            }
        }
    }    
    for (auto n12 : Neighbours12Diagonal) //������ �� ��������
    {
        if (hasNeighbour12(n12)) //���� �����
        {
            auto& rnnode = nodeRef(neighbour12(n12));
            rnnode.setNeighbour12(opposite(n12), ctags(toQuadrant(n12)));
            crefs(toQuadrant(n12)).setNeighbour12(n12, neighbour12(n12));
        }
    }
    /*
    //��������� �����: �������� ���������� �����
    nodeEdge nedge = { cTL.tag(), cTR.tag(), ORIENTATION_VERTICAL }; //����� �����
    edgeId eid = forest.addEdge(nedge); //�������� � ������
    cTL.edges[NEIGHBOUR_RIGHT * 2] = eid; //�������� � ������
    cTR.edges[NEIGHBOUR_LEFT * 2] = eid;
    nedge = { cTR.tag(), cBR.tag(), ORIENTATION_HORIZONTAL };
    eid = forest.addEdge(nedge);
    cTR.edges[NEIGHBOUR_BOTTOM * 2] = eid;
    cBR.edges[NEIGHBOUR_TOP * 2] = eid;
    nedge = { cBL.tag(), cBR.tag(), ORIENTATION_VERTICAL };
    eid = forest.addEdge(nedge);
    cBR.edges[NEIGHBOUR_LEFT * 2] = eid;
    cBL.edges[NEIGHBOUR_RIGHT * 2] = eid;
    nedge = { cTL.tag(), cBL.tag(), ORIENTATION_HORIZONTAL };
    eid = forest.addEdge(nedge);
    cTL.edges[NEIGHBOUR_BOTTOM * 2] = eid;
    cBL.edges[NEIGHBOUR_TOP * 2] = eid;
    //������� �����
    for (auto n : Neighbours)
    {
        if (neighbours[n].id != null)
        {
            auto& nnode = getNode(neighbours[n]);
            ushorty eindex = n * 2; //������ ����� (��� ����� ������)
            ushorty eindex_opposite = ((int)(n + 2) % 4) * 2; //������ ���������������� ����� (��� ������)
            if (nnode.is_leaf) //� ������ ��� �����, ����� ��������� �����
            {
                switch (n)
                {
                case NEIGHBOUR_TOP:
                    forest.updateEdge(edges[eindex], neighbours[n], cTL.tag()); //������ ������ ���������� ����� ������ �������� �������
                    cTL.edges[eindex] = edges[eindex]; //������ � ����� � �������
                    nnode.edges[eindex_opposite + 1] = edges[eindex]; //������ ����� � �������� ������
                    cTR.edges[eindex] = forest.addEdge({ neighbours[n], cTR.tag(), ORIENTATION_HORIZONTAL }); //�������� ������� ����� � ������ � ��� � �������
                    nnode.edges[eindex_opposite] = cTR.edges[eindex]; //������ � ������ ����� � �������� ������
                    break;
                case NEIGHBOUR_RIGHT:
                    forest.updateEdge(edges[eindex], cTR.tag(), neighbours[n]);
                    cTR.edges[eindex] = edges[eindex];
                    nnode.edges[eindex_opposite + 1] = edges[eindex];
                    cBR.edges[eindex] = forest.addEdge({ cBR.tag(), neighbours[n], ORIENTATION_VERTICAL });
                    nnode.edges[eindex_opposite] = cBR.edges[eindex];
                    break;
                case NEIGHBOUR_BOTTOM:
                    forest.updateEdge(edges[eindex], cBR.tag(), neighbours[n]);
                    cBR.edges[eindex] = edges[eindex];
                    nnode.edges[eindex_opposite + 1] = edges[eindex];
                    cBL.edges[eindex] = forest.addEdge({ cBL.tag(), neighbours[n], ORIENTATION_HORIZONTAL });
                    nnode.edges[eindex_opposite] = cBL.edges[eindex];
                    break;
                case NEIGHBOUR_LEFT:
                    forest.updateEdge(edges[eindex], neighbours[n], cBL.tag());
                    cBL.edges[eindex] = edges[eindex];
                    nnode.edges[eindex_opposite + 1] = edges[eindex];
                    cTL.edges[eindex] = forest.addEdge({ neighbours[n], cTL.tag(), ORIENTATION_VERTICAL });
                    nnode.edges[eindex_opposite] = cTL.edges[eindex];
                    break;
                }
            }
            else //���� ����, ����� ��� ���������, ����� ������ �������� ������ � ������ � ������ �� ����� �����
            {
                switch (n)
                {
                case NEIGHBOUR_TOP:
                    forest.updateEdge(edges[eindex], {}, cTL.tag()); //��������� ���������� �����
                    cTL.edges[eindex] = edges[eindex]; //������ � ����� � �������
                    forest.updateEdge(edges[eindex + 1], {}, cTR.tag()); //������ �����
                    cTR.edges[eindex] = edges[eindex + 1];
                    break;
                case NEIGHBOUR_RIGHT:
                    forest.updateEdge(edges[eindex], cTR.tag(), {});
                    cTR.edges[eindex] = edges[eindex];
                    forest.updateEdge(edges[eindex + 1], cBR.tag(), {});
                    cBR.edges[eindex] = edges[eindex + 1];
                    break;
                case NEIGHBOUR_BOTTOM:
                    forest.updateEdge(edges[eindex], cBR.tag(), {});
                    cBR.edges[eindex] = edges[eindex];
                    forest.updateEdge(edges[eindex + 1], cBL.tag(), {});
                    cBL.edges[eindex] = edges[eindex + 1];
                    break;
                case NEIGHBOUR_LEFT:
                    forest.updateEdge(edges[eindex], {}, cBL.tag());
                    cBL.edges[eindex] = edges[eindex];
                    forest.updateEdge(edges[eindex + 1], {}, cTL.tag());
                    cTL.edges[eindex] = edges[eindex + 1];
                    break;
                }
            }
        }
    }
    //������� ������� � ������ ������ (�������������?)
    for (size_t i = 0; i < 8; i++)
        edges[i] = null;
    */
    return 0;
}

void TreeNode::gatherDataFromChildren() //���� ������ �� ����� ��� �����������
{
    try {
        if (dataId_ == null)
            throw std::invalid_argument("node is deleted");
        else
        {
            auto& rd = dataRef();
            rd.clear(); //���������
            //���������� �������������� ������� �� �����
            for (auto q : Quadrants)
            {
                CellData cd = childRef(q).data();
                rd.add(cd);
            }
            rd.divide(static_cast<double>(QUADRANTS_NUM));
            if (config.coord_type == CoordType::axisymmetric) //�������� r � cellData ��� ��������������� ���������
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
int TreeNode::coarsen() //������� ������
{
    //��������� ������, ����� ��� ������� ������ ���������
    if (is_deleted)
        return INT_ERROR_CODE_CANT_COARSEN_DELETED;
    if (is_leaf)
        return INT_ERROR_CODE_CANT_COARSEN_LEAF;
    if (has_grandchildren)
        return INT_ERROR_CODE_CANT_COARSEN_GRANDPARENT;

    //�������� �� �����c
    if (!forest.treeRef(tag().tree()).isGhost()) //� ghost-�������� ������ �� �����
    {
        for (auto n : Neighbours)
        {
            if (hasNeighbour(n) && nodeRef(neighbour(n)).hasGrandChildren(opposite(n))) //� ������ ���� ����� � ������� �������
            {
                    return INT_ERROR_CODE_CANT_COARSEN_UNBALANCED;
            }
        }
    }

    auto& rtree = forest.treeRef(tag().tree());
    //��������� ������ ������ ������ ��� ������� ���� � ���� ������ � ��
    dataId_ = rtree.getVacantDataId();
    gatherDataFromChildren();

    //���������� ������ � �������� �������
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n))
        {
            nodeRef(neighbour(n)).setNeighbour(opposite(n), tag());
        }
    }

    /*
    auto& cTL = getChild(TOP_LEFT);
    auto& cTR = getChild(TOP_RIGHT);
    auto& cBR = getChild(BOTTOM_RIGHT);
    auto& cBL = getChild(BOTTOM_LEFT);
    auto index = Neighbour12::top1;
    if (cTL.hasNeighbour12(index))
    {
        if (cTL.neighbours12[index] == cTR.neighbours12[index]) //����� ������� ����� ������ ������
        {
            auto& nTL = getNode(cTL.neighbours12[index]);
            nTL.updateNeighbour12(Opposites12[Nexts12[index]], tag()); //���������� ��� ������
            nTL.updateNeighbour12(Opposites12[index], {}); //������ "null"
            updateNeighbour12(index, cTL.neighbours12[index]); //��� ������ ������
            updateNeighbour12(Nexts12[index], {});
        }
        else //������ ���� �� �������, ��� ���� ������ ����
        {
            auto& nTL = getNode(cTL.neighbours12[index]);
            auto& nTR = getNode(cTR.neighbours12[index]);
            nTL.updateNeighbour12(Opposites12[Nexts12[index]], tag()); //��� �������
            nTR.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            updateNeighbour12(index, cTL.neighbours12[index]); //��� ������ ������
            updateNeighbour12(Nexts12[index], cTR.neighbours12[index]);
            nTL.updateNeighbour12(Opposites12[Prevs12[index]], {}); //���������
            nTR.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], {});
        }
    }
    index = Neighbour12::top_right;
    if (cTR.hasNeighbour12(index))
    {
        auto& nDTR = getNode(cTR.neighbours12[index]);
        nDTR.updateNeighbour12(Opposites12[index], tag());
        updateNeighbour12(index, nDTR.tag());
    }
    index = eighbour12::right1;
    if (cTR.hasNeighbour12(index))
    {
        if (cTR.neighbours12[index] == cBR.neighbours12[index])
        {
            auto& nTR = getNode(cTR.neighbours12[index]);
            nTR.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            nTR.updateNeighbour12(Opposites12[index], {});
            updateNeighbour12(index, cTR.neighbours12[index]);
            updateNeighbour12(Nexts12[index], {});
        }
        else
        {
            auto& nTR = getNode(cTR.neighbours12[index]);
            auto& nBR = getNode(cBR.neighbours12[index]);
            nTR.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            nBR.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            updateNeighbour12(index, cTR.neighbours12[index]);
            updateNeighbour12(Nexts12[index], cBR.neighbours12[index]);
            nTR.updateNeighbour12(Opposites12[Prevs12[index]], {});
            nBR.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], {});
        }
    }
    index = Neighbour12::bottom_right;
    if (cBR.hasNeighbour12(index))
    {
        auto& nDBR = getNode(cBR.neighbours12[index]);
        nDBR.updateNeighbour12(Opposites12[index], tag());
        updateNeighbour12(index, nDBR.tag());
    }
    index = eighbour12::bottom1;
    if (cBR.hasNeighbour12(index))
    {
        if (cBR.neighbours12[index] == cBL.neighbours12[index])
        {
            auto& nBR = getNode(cBR.neighbours12[index]);
            nBR.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            nBR.updateNeighbour12(Opposites12[index], {});
            updateNeighbour12(index, cBR.neighbours12[index]);
            updateNeighbour12(Nexts12[index], {});
        }
        else
        {
            auto& nBR = getNode(cBR.neighbours12[index]);
            auto& nBL = getNode(cBL.neighbours12[index]);
            nBR.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            nBL.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            updateNeighbour12(index, cBR.neighbours12[index]);
            updateNeighbour12(Nexts12[index], cBL.neighbours12[index]);
            nBR.updateNeighbour12(Opposites12[Prevs12[index]], {});
            nBL.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], {});
        }
    }
    index = Neighbour12::bottom_left;
    if (cBL.hasNeighbour12(index))
    {
        auto& nDBL = getNode(cBL.neighbours12[index]);
        nDBL.updateNeighbour12(Opposites12[index], tag());
        updateNeighbour12(index, nDBL.tag());
    }
    index = Neighbour12::left1;
    if (cBL.hasNeighbour12(index))
    {
        if (cBL.neighbours12[index] == cTL.neighbours12[index])
        {
            auto& nBL = getNode(cBL.neighbours12[index]);
            nBL.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            nBL.updateNeighbour12(Opposites12[index], {});
            updateNeighbour12(index, cBL.neighbours12[index]);
            updateNeighbour12(Nexts12[index], {});
        }
        else
        {
            auto& nBL = getNode(cBL.neighbours12[index]);
            auto& nTL = getNode(cTL.neighbours12[index]);
            nBL.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            nTL.updateNeighbour12(Opposites12[Nexts12[index]], tag());
            updateNeighbour12(index, cBL.neighbours12[index]);
            updateNeighbour12(Nexts12[index], cTL.neighbours12[index]);
            nBL.updateNeighbour12(Opposites12[Prevs12[index]], {});
            nTL.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], {});
        }
    }
    index = Neighbour12::top_left;
    if (cTL.hasNeighbour12(index))
    {
        auto& nDTL = getNode(cTL.neighbours12[index]);
        nDTL.updateNeighbour12(Opposites12[index], tag());
        updateNeighbour12(index, nDTL.tag());
    }


    //��������� �����
    forest.removeEdge(cTL.edges[NEIGHBOUR_RIGHT * 2]); //������ ����� = ������ ������ * 2
    forest.removeEdge(cTL.edges[NEIGHBOUR_BOTTOM * 2]);
    forest.removeEdge(cBR.edges[NEIGHBOUR_TOP * 2]);
    forest.removeEdge(cBR.edges[NEIGHBOUR_LEFT * 2]);
    //���������� ������� �����
    for (auto n : Neighbours)
    {
        if (neighbours[n].id != null)
        {
            auto& nnode = getNode(neighbours[n]);
            ushorty eindex = n * 2; //������ ����� (��� ������ ������)
            ushorty eindex_opposite = ((int)(n + 2) % 4) * 2; //������ ���������������� ����� (��� ������)
            if (nnode.is_leaf) //���� � ������ ��� ����� - ��������� � ��������� ������ � ����� � � �����
            {
                edgeId nedgeid = null; //id �����, � ������� ����� ������������ ��������� (������������ ����)
                switch (n)
                {
                case NEIGHBOUR_TOP:
                    nedgeid = cTL.edges[eindex];
                    forest.updateEdge(nedgeid, neighbours[n], ctag); //���������� ������ � ��������� ����� ������ ������� �������
                    forest.removeEdge(cTR.edges[eindex]); //�������� ������� �����
                    break;
                case NEIGHBOUR_RIGHT:
                    nedgeid = cTR.edges[eindex];
                    forest.updateEdge(nedgeid, ctag, neighbours[n]);
                    forest.removeEdge(cBR.edges[eindex]);
                    break;
                case NEIGHBOUR_BOTTOM:
                    nedgeid = cBR.edges[eindex];
                    forest.updateEdge(nedgeid, ctag, neighbours[n]);
                    forest.removeEdge(cBL.edges[eindex]);
                    break;
                case NEIGHBOUR_LEFT:
                    nedgeid = cBL.edges[eindex];
                    forest.updateEdge(nedgeid, neighbours[n], ctag);
                    forest.removeEdge(cTL.edges[eindex]);
                    break;
                }
                edges[eindex] = nedgeid; //������ � ����� � ������ ������
                edges[eindex + 1] = null; //������ �����
                nnode.edges[eindex_opposite] = nedgeid; //� �������� ������
                nnode.edges[eindex_opposite + 1] = null; //������ �����
            }
            else //� ������ ���� ���� - ������ ��������� ������ � ������ � � ������ 
            {
                switch (n)
                {
                case NEIGHBOUR_TOP:
                    edges[eindex] = cTL.edges[eindex]; //������ ����� ������� � ������ ������
                    edges[eindex + 1] = cTR.edges[eindex]; //������ �����
                    forest.updateEdge(edges[eindex], {}, ctag); //���������� ������ ���� � ����� �� ��� ������ ������
                    forest.updateEdge(edges[eindex + 1], {}, ctag); //������ �����
                    break;
                case NEIGHBOUR_RIGHT:
                    edges[eindex] = cTR.edges[eindex];
                    edges[eindex + 1] = cBR.edges[eindex];
                    forest.updateEdge(edges[eindex], ctag, {});
                    forest.updateEdge(edges[eindex + 1], ctag, {});
                    break;
                case NEIGHBOUR_BOTTOM:
                    edges[eindex] = cBR.edges[eindex];
                    edges[eindex + 1] = cBL.edges[eindex];
                    forest.updateEdge(edges[eindex], ctag, {});
                    forest.updateEdge(edges[eindex + 1], ctag, {});
                    break;
                case NEIGHBOUR_LEFT:
                    edges[eindex] = cBL.edges[eindex];
                    edges[eindex + 1] = cTL.edges[eindex];
                    forest.updateEdge(edges[eindex], {}, ctag);
                    forest.updateEdge(edges[eindex + 1], {}, ctag);
                    break;
                }
            }
        }
    }*/

    //�������� �����
    rtree.vacateNodeGroup(tag().depth()+1, childrenId_); //������� ������ ����� ��� ��� ����������
    for (auto q : Quadrants)
    {
        auto& rchild = childRef(q);
        rchild.markDeleted();
        rtree.vacateData(rchild.dataId()); //������� ����� � ������� tree.data ����������
    }
    rtree.incrementNodesCounter(tag().depth() + 1, -QUADRANTS_NUM); //�������� ��� �� ��������� ������ ����� �� 4 ������
    rtree.incrementLeavesCounter(tag().depth() + 1, -QUADRANTS_NUM); //������� �� ��������� ������ ����� �� 4 ������
    rtree.incrementLeavesCounter(tag().depth(), 1); //������� �� ������� ������ ����� �� 1 ������
        
    //������ ����� ������
    is_leaf = true;
    childrenId_ = null;

    //���������� ������ � ������ ��������
    if (parentId_ != null)
        nodeRef(NodeTag(tag().tree(), tag().depth() - 1, parentId_)).updateGrandParency();

    //�������� ���� ������ ��� �������������
    rtree.deleteLevelIfEmpty(tag().depth());

    //������������ ������� �������� ghost'��
    for (auto n : Neighbours)
    {
        if (hasNeighbour(n))
        {
            if (!forest.treeRef(tag().tree()).isGhost() && forest.treeRef(neighbour(n).tree()).isGhost()) //���� ghost-�����
            {
                auto& rnode = nodeRef(neighbour(n));
                if (rnode.hasChildren())
                    rnode.coarsen();
            }
        }
    }
    return 0;
}

//inspectors -----------------------
bool TreeNode::isDeleted() const { return is_deleted; }
bool TreeNode::isLeaf() const { return is_leaf; }
bool TreeNode::hasChildren() const { return childrenId_ != null; } //���� �� ����
bool TreeNode::hasGrandChildren() const { return has_grandchildren; } //���� �� �����
bool TreeNode::hasGrandChildren(Neighbour n) const //���� �� ����� � ������������ �������
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
bool TreeNode::hasNeighbour(Neighbour n) const //���� �� ����� �� �����������
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
//uncategorized -----------------------
TreeNode& TreeNode::nodeRef(const NodeTag& tag) //������ �� ���� �� ����
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
    return nodeRef({}); //to suppress warning
}

double TreeNode::magGradRho() const //��������� �������� ���������
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
            const auto& rndata = rnnode.dataRefConst(); //����� �� ������, �.�. � neighbour12 ����� ������ ������
            drhosum += fabs(rd.rho() - rndata.rho()) / distance(box().center(), rnnode.box().center());
            n12num++;
        }
    }
    return drhosum / n12num / rd.rho();
}

//output -----------------------
std::string TreeNode::dump() const //���� ���� � ������
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
            auto& nnode = getNode(neighbours[n]);
            Neighbour dir = (Neighbour)((int)(n + 2) % 4); //��������������� �����������
            buffer << Dirs[dir] << ": " << std::setw(3) << nnode.neighbours[dir].tree_id << ", " << nnode.neighbours[dir].depth << ", " << std::setw(3) << nnode.neighbours[dir].id << " | ";
        }
        else
            buffer << "               | ";
    }
    buffer << "]" << endl;*/

    return buffer.str();
}
const double Dx[] = { 0, 1, 0, -1 }; //��� ���������� ��������
const double Dy[] = { 1, 0, -1, 0 };
std::string TreeNode::dumpNeighbourVector(Neighbour n) const //���� ������ � ���� �������
{
    if (!hasNeighbour(n))
        return "(no such neighbour)";

    std::stringstream buffer;
    double x = box().center().x + 0.8 * Dx[static_cast<int>(n)] * box().size() / 2.0; //��������� ����� �������
    double y = box().center().y + 0.8 * Dy[static_cast<int>(n)] * box().size() / 2.0;
    auto& rnnode = nodeRef(neighbour(n));
    double dx = rnnode.box().center().x - x; //������
    double dy = rnnode.box().center().y - y;

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
