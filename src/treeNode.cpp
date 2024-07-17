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
NodeTag TreeNode::tag() const { return tag_; } //��������� ���� 
cellDataId TreeNode::dataId() const { return dataId_; } //��������� dataId
CellBox TreeNode::box() const {	return box_; } //��������� box'a
//QuadTree& TreeNode::treeRef() { return forest.treeRef(tag().tree()); } // ������ �� ������, ���������� ����
//const QuadTree& TreeNode::treeRefConst() const { return forest.treeRefConst(tag().tree()); } //(const) ������ �� ������, ���������� ����
NodeTag TreeNode::neighbour(Neighbour n) const { return neighbours[static_cast<int>(n)]; }
NodeTag TreeNode::neighbour12(Neighbour n12) const { return neighbours12[static_cast<int>(n12)]; }

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
        cout << "getDataRef error: " << e.what() << ", tag = " << tag() << endl;
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
        cout << "getDataRef error: " << e.what() << ", tag = " << tag() << endl;
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
                throw std::invalid_argument("null data reference in this node");
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
        cout << "getData error: " << e.what() << ", tag = " << tag() << endl;
        return forest.treeRef(0).data(0); //to suppress warning
    }
}

const int ERROR_NODE_NO_CHILDREN = -1;
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
            throw ERROR_NODE_NO_CHILDREN;
        }
    }
    catch (int err_code)
    {
        cout << "node.childRef error: " << err_code << " " << tag().tree() << " " << tag().depth() + 1 << " " << childrenId_ + static_cast<int>(q) << endl;
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
            throw ERROR_NODE_NO_CHILDREN;
        }
    }
    catch (int err_code)
    {
        cout << "node.childRef error: " << err_code << " " << tag().tree() << " " << tag().depth() + 1 << " " << childrenId_ + static_cast<int>(q) << endl;
        return *this; //to suppress warning
    }
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
void TreeNode::setDataId(cellDataId id) { dataId_ = id; } //������� dataId
void TreeNode::setBox(const CellBox& b) { box_ = b; } //������� box'�
void TreeNode::setNeighbour(Neighbour n, NodeTag t) { neighbours[static_cast<int>(n)] = t; }
void TreeNode::setNeighbour12(Neighbour12 n12, NodeTag t) { neighbours12[static_cast<int>(n12)] = t; }

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
    case Neighbour::right: //�� ������ ������
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

const int INT_ERROR_CODE_CANT_REFINE_DELETED = -1;
const int INT_ERROR_CODE_CANT_REFINE_DEEPEST_LEVEL = -2;
const int INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED = -3;
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
            //(?)�������� �� ghost-������ ������ ��� �� �����, �.�. ��� ������ ������������� ����
            if (!rntag.isNull() && rntag.depth() < tag().depth()) //����� �� 2+ ������ ������ ������� ����� ������ ����
                nodeRef(rntag).markToRefine();
            //������������ ������� �������� ghost'��
            if (forest.treeRef(rntag.tree()).isGhost() && rntag.depth() <= tag().depth()) //ghost-����� ������� ������� ����� ������ ����
                nodeRef(rntag).markToRefine();
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
    rtree.incrementLeavesCounter(tag().depth(), -1); //������� �� ������� ������ ����� ������ �� 1
    rtree.incrementLeavesCounter(tag().depth() + 1, QUADRANTS_NUM); //������� �� ��������� ������ ����� ������ �� 4

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
            if(rnnode.hasChildren()) //���������� ���� ����� �����, ������ ���� � ������ ������ ���� ����
                rnnode.setChildrenNeighbours(opposite(n), ctags);
        }
    }



    /*
    //���������� ������ � ������� � ������ ������������ � ������� �����
    cTL.neighbours12[N12_RIGHT_TOP] = cTR.tag(); //siblings
    cTL.neighbours12[N12_BOTTOM_RIGHT] = cBL.tag();
    cTL.neighbours12[N12_DIAG_RIGHT_BOTTOM] = cBR.tag();
    cTR.neighbours12[N12_LEFT_BOTTOM] = cTL.tag();
    cTR.neighbours12[N12_BOTTOM_RIGHT] = cBR.tag();
    cTR.neighbours12[N12_DIAG_BOTTOM_LEFT] = cBL.tag();
    cBR.neighbours12[N12_TOP_LEFT] = cTR.tag();
    cBR.neighbours12[N12_LEFT_BOTTOM] = cBL.tag();
    cBR.neighbours12[N12_DIAG_LEFT_TOP] = cTL.tag();
    cBL.neighbours12[N12_TOP_LEFT] = cTL.tag();
    cBL.neighbours12[N12_RIGHT_TOP] = cBR.tag();
    cBL.neighbours12[N12_DIAG_TOP_RIGHT] = cTR.tag();

    auto index = N12_TOP_LEFT;
    if (hasNeighbour12(index))
    {
        if (hasNeighbour12(Nexts12[index])) //�������� ������ ������, ��� ������
        {
            auto& nTL = getNode(neighbours12[index]);
            auto& nTR = getNode(neighbours12[Nexts12[index]]);
            nTL.updateNeighbour12(Opposites12[Nexts12[index]], cTL.tag()); //���������� ��� ������
            cTL.updateNeighbour12(index, nTL.tag()); //��� �������
            nTR.updateNeighbour12(Opposites12[Nexts12[index]], cTR.tag());
            cTR.updateNeighbour12(index, nTR.tag());
            nTL.updateNeighbour12(Opposites12[Prevs12[index]], cTR.tag()); //���������
            cTR.updateNeighbour12(Prevs12[index], nTL.tag());
            nTR.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], cTL.tag());
            cTL.updateNeighbour12(Nexts12[Nexts12[index]], nTR.tag());
        }
        else //����� ���� �� �������, ��� ������ ������ �� ���������
        {
            auto& nTL = getNode(neighbours12[index]);
            nTL.updateNeighbour12(Opposites12[index], cTL.tag()); //���������� ��� ������
            nTL.updateNeighbour12(Opposites12[Nexts12[index]], cTR.tag());
            cTL.updateNeighbour12(index, nTL.tag()); //��� �����
            cTR.updateNeighbour12(index, nTL.tag());
        }
    }
    index = N12_DIAG_TOP_RIGHT;
    if (hasNeighbour12(index))
    {
        auto& nDTR = getNode(neighbours12[index]);
        nDTR.updateNeighbour12(Opposites12[index], cTR.tag()); //���������� ��� ������
        cTR.updateNeighbour12(index, nDTR.tag()); //��� �������
    }
    index = N12_RIGHT_TOP;
    if (hasNeighbour12(index))
    {
        if (hasNeighbour12(Nexts12[index]))
        {
            auto& nTR = getNode(neighbours12[index]);
            auto& nBR = getNode(neighbours12[Nexts12[index]]);
            nTR.updateNeighbour12(Opposites12[Nexts12[index]], cTR.tag());
            cTR.updateNeighbour12(index, nTR.tag());
            nBR.updateNeighbour12(Opposites12[Nexts12[index]], cBR.tag());
            cBR.updateNeighbour12(index, nBR.tag());
            nTR.updateNeighbour12(Opposites12[Prevs12[index]], cBR.tag());
            cBR.updateNeighbour12(Prevs12[index], nTR.tag());
            nBR.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], cTR.tag());
            cTR.updateNeighbour12(Nexts12[Nexts12[index]], nBR.tag());
        }
        else
        {
            auto& nTR = getNode(neighbours12[index]);
            nTR.updateNeighbour12(Opposites12[index], cTR.tag());
            nTR.updateNeighbour12(Opposites12[Nexts12[index]], cBR.tag());
            cTR.updateNeighbour12(index, nTR.tag());
            cBR.updateNeighbour12(index, nTR.tag());
        }
    }
    index = N12_DIAG_RIGHT_BOTTOM;
    if (hasNeighbour12(index))
    {
        auto& nDBR = getNode(neighbours12[index]);
        nDBR.updateNeighbour12(Opposites12[index], cBR.tag());
        cBR.updateNeighbour12(index, nDBR.tag());
    }
    index = N12_BOTTOM_RIGHT;
    if (hasNeighbour12(index))
    {
        if (hasNeighbour12(Nexts12[index]))
        {
            auto& nBR = getNode(neighbours12[index]);
            auto& nBL = getNode(neighbours12[Nexts12[index]]);
            nBR.updateNeighbour12(Opposites12[Nexts12[index]], cBR.tag());
            cBR.updateNeighbour12(index, nBR.tag());
            nBL.updateNeighbour12(Opposites12[Nexts12[index]], cBL.tag());
            cBL.updateNeighbour12(index, nBL.tag());
            nBR.updateNeighbour12(Opposites12[Prevs12[index]], cBL.tag());
            cBL.updateNeighbour12(Prevs12[index], nBR.tag());
            nBL.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], cBR.tag());
            cBR.updateNeighbour12(Nexts12[Nexts12[index]], nBL.tag());
        }
        else
        {
            auto& nBR = getNode(neighbours12[index]);
            nBR.updateNeighbour12(Opposites12[index], cBR.tag());
            nBR.updateNeighbour12(Opposites12[Nexts12[index]], cBL.tag());
            cBR.updateNeighbour12(index, nBR.tag());
            cBL.updateNeighbour12(index, nBR.tag());
        }
    }
    index = N12_DIAG_BOTTOM_LEFT;
    if (hasNeighbour12(index))
    {
        auto& nDBL = getNode(neighbours12[index]);
        nDBL.updateNeighbour12(Opposites12[index], cBL.tag());
        cBL.updateNeighbour12(index, nDBL.tag());
    }
    index = N12_LEFT_BOTTOM;
    if (hasNeighbour12(index))
    {
        if (hasNeighbour12(Nexts12[index]))
        {
            auto& nBL = getNode(neighbours12[index]);
            auto& nTL = getNode(neighbours12[Nexts12[index]]);
            nBL.updateNeighbour12(Opposites12[Nexts12[index]], cBL.tag());
            cBL.updateNeighbour12(index, nBL.tag());
            nTL.updateNeighbour12(Opposites12[Nexts12[index]], cTL.tag());
            cTL.updateNeighbour12(index, nTL.tag());
            nBL.updateNeighbour12(Opposites12[Prevs12[index]], cTL.tag());
            cTL.updateNeighbour12(Prevs12[index], nBL.tag());
            nTL.updateNeighbour12(Opposites12[Nexts12[Nexts12[index]]], cBL.tag());
            cBL.updateNeighbour12(Nexts12[Nexts12[index]], nTL.tag());
        }
        else
        {
            auto& nBL = getNode(neighbours12[index]);
            nBL.updateNeighbour12(Opposites12[index], cBL.tag());
            nBL.updateNeighbour12(Opposites12[Nexts12[index]], cTL.tag());
            cBL.updateNeighbour12(index, nBL.tag());
            cTL.updateNeighbour12(index, nBL.tag());
        }
    }
    index = N12_DIAG_LEFT_TOP;
    if (hasNeighbour12(index))
    {
        auto& nDTL = getNode(neighbours12[index]);
        nDTL.updateNeighbour12(Opposites12[index], cTL.tag());
        cTL.updateNeighbour12(index, nDTL.tag());
    }

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

//inspectors -----------------------
bool TreeNode::isDeleted() const { return is_deleted; }
bool TreeNode::isLeaf() const { return is_leaf; }
bool TreeNode::hasChildren() const { return childrenId_ != null; } //���� �� ����
bool TreeNode::hasGrandChildren() const { return has_grandchildren; } //���� �� �����
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
    default: return rTL; //to suppress warning
    }
}
