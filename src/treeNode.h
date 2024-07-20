#ifndef qtCWENO_TreeNode_H //include guard
#define qtCWENO_TreeNode_H

#include "main.h"
#include "CellData.h"
#include "NodeEdge.h"
#include "NodeTag.h"
#include "QuadTree.h"

#include "Eigen\Dense"

const int INT_ERROR_CODE_CANT_REFINE_DELETED = -1;
const int INT_ERROR_CODE_CANT_REFINE_DEEPEST_LEVEL = -2;
const int INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED = -3;

//struct ���� �����
struct ChildrenTags
{
    std::array<NodeTag, QUADRANTS_NUM> tags = {};
    ChildrenTags(NodeTag tTL, NodeTag tTR, NodeTag tBR, NodeTag tBL) : tags({ tTL, tTR, tBR, tBL }) {};
    const NodeTag& operator()(Quadrant q) const; //() operator: ��� �� ���������
};

//���� ������
class TreeNode
{
    NodeTag tag_ = {}; //{} �������� ��������� ����������� NodeTag
    treeNodeId parentId_ = null; //����� �������� �� ���������� ������
    treeNodeId childrenId_ = null; //����� ������� child �� ��������� ������, ��������� ���� ����� �� ���
    bool has_grandchildren = false; //���� �� ����� (����� ��� ������������)
    bool is_leaf = true; //�������� �� ������ ������ (������ ��� ��� �������� ���������� ������) (���������, �� ����� ������, ��� ��������� ������� �����)
    bool is_deleted = false; //��� ������� ����� ����� �� ���������, �.�. �������� = ������� ����� �� ������ ���������
    cellDataId dataId_ = null; //����� ������ ��� ������ � ���������� ������� cell_data, ������ ��� �������
    CellBox box_ = {};
    std::array<NodeTag, NEIGHBOURS_NUM> neighbours = {}; //������
    std::array<NodeTag, MAX_NEIGHBOURS12_NUM> neighbours12 = {}; //������ � ������ ������������ � ���������� ���������� ����� ������
    //������ �����, ������������ ������� ������ (�� DIRECTIONS_NUM*2 ���� � ������ ���������� ���������� �����)
    
    std::array<nodeEdgeId, MAX_NEIGHBOURS12_NUM - 4> edges = {}; //��������� �� ������� �������, ������� � ����� �������� ������� �������
    //�������� �� ��������� ��� edges[8] ����� ��������� ����, �.�. ��� ���������� ������������ ��� nodeEdgeId
    double polyCoeffs[EQ_NUM][POLY_COEFF_NUM] = {}; //������������ ����������� CWENO ��� ������� ���������

    Eigen::VectorXd coeffs; //������ ����������� (px, py, pxx, pxy, pyy) ��� ������������ ��������
    Eigen::MatrixXd J; //������� �������
    Eigen::VectorXd rs; //������ ������ ������ �������
    Eigen::HouseholderQR<Eigen::MatrixXd> decomp; //��������
    Eigen::VectorXd coeffsl; //������ ����������� (px, py) ��� �������� �������
    Eigen::MatrixXd Jl[4]; //�� ����� ��� ������� ����������
    Eigen::VectorXd rsl[4];
    Eigen::HouseholderQR<Eigen::MatrixXd> decompl[4];

public:
    TreeNode() {}; //default constructor
    explicit TreeNode(NodeTag t, treeNodeId p, cellDataId d, CellBox b) : tag_(t), parentId_(p), dataId_(d), box_(b) {}; //����������� �� ��������� ������
    
    //accessors
    const NodeTag tag() const;
    cellDataId dataId() const;
    CellBox box() const;
    //QuadTree& treeRef(); //������ �� ������, ���������� ���� TODO:�����������, ������ �� �������� (circular dependency?)
    //const QuadTree& treeRefConst() const; //const ������
    const NodeTag neighbour(Neighbour n) const; 
    const NodeTag neighbour12(Neighbour12 n12) const; 
    CellData& dataRef(); //������ �� ������
    const CellData& dataRefConst() const; //const ������
    CellData data() const; //����� ������ (� ��������� ������ � �����)
    TreeNode& childRef(Quadrant q); //������ �� ������� �� ���������
    const TreeNode& childRef(Quadrant q) const; //const ������
    const ChildrenTags childrenTags() const; //���� ���� �����

    TreeNode& getChildOrSelfByCoords(Point p); //������ �� ������� (��� ����) �� �����������

    //mutators
    void setTag(const NodeTag& t);
    void markDeleted(); //������� ���������
    void setDataId(cellDataId id);
    void setBox(const CellBox& b);
    void setNeighbour(Neighbour n, NodeTag ntag); //������� ������ ������ � ����� � ������ �������
    void setChildrenNeighbours(Neighbour n, ChildrenTags tags); //�������� ������ � ������� ��� �����
    void clearNeighbours12(); //�������� ���� ������� � �������12 ����� ���������
    void setNeighbour12(Neighbour12 n12, NodeTag ntag);
    void setChildrenOrSelfNeighbour12(Neighbour n, NodeTag ntag); //�������� ������ � ������12 ��� ���� ��� ����� ����� ������� ������
    void setNeighbours12(Neighbour n, ChildrenTags tags); //�������� ������ � �������12 ��� ����
    void setChildrenOrSelfNeighbours12(Neighbour n, ChildrenTags tags); //�������� ������ � �������12 ��� ���� ��� �����
    void setChildrenCommonNeighbour12(Neighbour n, NodeTag ntag); //�������� ������ �� ����� ������12 ��� �����
    void setData(const CellData& data);
    void setGrandParency(bool status); //������� � ������� ������
    void updateGrandParency(); //���������� ������ � ������� ������ ����� ������� �������
    //void setChildrenCommonNeighbour(Neighbour n, NodeTag ntag); //�������� ������ �� ����� ������ ��� �����
    int markToRefine(); //������� ������ � ���������
    int refine(); //��������� ������
    void gatherDataFromChildren(); //���� ������ �� ����� ��� �����������
    int coarsen(); //������� ������

    //inspectors
    bool isDeleted() const;
    bool isLeaf() const;
    bool hasChildren() const; //���� �� ����
    bool hasGrandChildren() const; //���� �� �����
    bool hasGrandChildren(Neighbour n) const; //���� �� ����� � ������������ �������
    bool hasNeighbour(Neighbour n) const; //���� �� ����� �� �����������
    bool hasNeighbour12(Neighbour12 n12) const;

    //other
    static TreeNode& nodeRef(const NodeTag& tag); //������ �� ���� �� ���� (static - ����� ������� ��� ���� ���)
    double magGradRho() const; //��������� �������� ���������

    //output
    std::string dump() const; //���� ���� � ������
    std::string dumpNeighbourVector(Neighbour n) const; //���� ������ � ���� �������
    std::string dumpNeighbour12Vector(Neighbour12 n12) const; //���� ������12 � ���� �������

    friend class QuadTreeForest; //��� ���������� ������� TODO:�����������, ��� ����� ����������� ��������� ����� ��� ��������� ������������
    friend class QuadTree;

    /*
    NodeTag getNodeOrChildTag(int target_depth, Quadrant quadrant) const; //����� ��������� ���� ������� ������ ������ ������ (������� ������ �� ����� 1)
    bool hasEdge(ushorty etype) const; //���� �� �����
    nodeEdge& getEdge(ushorty etype) const; //������ �� �����
    int neighbours12Num() const; //����� �������
    //bool isBoundary(); //�������� �� ���������
    //bool isCorner(); //�������� �� �������
    TreeNode& getNeigbourOrSelf(Neighbour dir, int target_depth) const; //������ �� ������ (��� �� ����, ���� ��� ������)
    bool isNeighbour12InSubstencil(Neighbour12 n, Quadrant q) const; //�������� �� ����� � ��������� ��� �������� �������
    void updateEigenObjects(); //�������� ��� ���������� Eigen ������ � �.�. ����� ��������� �����
    void calcPolynomialCWENO(rkStep rk); //���������� ������������� 2D CWENO ��������
    CellData evalPolynomialAt(point p, rkStep rk = 0); //������������������ ��������� �������� 
    */
};

//struct ������ �� �����
struct ChildrenRefs
{
    TreeNode& rTL;
    TreeNode& rTR;
    TreeNode& rBR;
    TreeNode& rBL;
    ChildrenRefs(TreeNode& _rTL, TreeNode& _rTR, TreeNode& _rBR, TreeNode& _rBL) : rTL(_rTL), rTR(_rTR), rBR(_rBR), rBL(_rBL) {};
    TreeNode& operator()(Quadrant q); //() operator: ������ �� ���������
};

#endif
