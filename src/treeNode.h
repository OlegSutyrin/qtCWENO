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
    std::array<NodeTag, NEIGHBOURS_NUM> neighbours = {}; //������(������������ � �������� ��� ��������� ������)
    std::array<NodeTag, MAX_NEIGHBOURS12_NUM> neighbours12 = {}; //������ � ������ ���������� � ��������� ����� (������������ ��� CWENO)
    //������ �����, ������������ ������� ������ (�� EDGES_NUM=8 ���� � ������ ���������� ���������� �����)
    std::array<nodeEdgeId, EDGES_NUM> edges = { null, null, null, null, null, null, null, null }; //��������� �� ������� �������, ������� � ����� �������� ������� �������
    //null �� ��������� ��� edges[EDGES_NUM] ����� ��������� ����, �.�. nodeEdgeId �� ������� = 0
    double polyCoeffs[EQ_NUM][POLY_COEFF_NUM] = {}; //������������ ����������� CWENO ��� ������� ���������

    Eigen::VectorXd coeffs; //������ ����������� (px, py, pxx, pxy, pyy) ��� ������������ ��������
    Eigen::MatrixXd J; //������� �������
    Eigen::VectorXd rs; //������ ������ ������ �������
    Eigen::HouseholderQR<Eigen::MatrixXd> decomp; //��������
    Eigen::VectorXd coeffsl; //������ ����������� (px, py) ��� �������� �������
    Eigen::MatrixXd Jl[QUADRANTS_NUM]; //�� ����� ��� ������� ����������
    Eigen::VectorXd rsl[QUADRANTS_NUM];
    Eigen::HouseholderQR<Eigen::MatrixXd> decompl[QUADRANTS_NUM];

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
    const nodeEdgeId edge(Edge e) const; //id �����
    NodeEdge& edgeRef(Edge e); //������ �� �����
    double polyCoeff(Equation eq, int p) const; //����������� �������� CWENO

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
    void setEdge(Edge etype, nodeEdgeId eid); //������� �����
    void splitEdge(Neighbour n); //��������� ����� � �������� ������ � ���� � (����������) ������
    void joinEdge(Neighbour n); //������� ����� � �������� ������ � ���� � (����������) ������
    void updateChildrenEdges(Neighbour n); //�������� ����� � ����� � ��������� �����
    void gatherEdgesFromChildren(Neighbour n); //�������� ������� ����� ���� � �������� ����� � ��������� �����
    int markToRefine(); //������� ������ � ���������
    int refine(); //��������� ������
    void gatherDataFromChildren(); //���� ������ �� ����� ��� �����������
    int coarsen(); //������� ������
    void updateEigenObjects(); //�������� ��� ���������� Eigen ������ � �.�. ����� ��������� �����
    void calcPolynomialCWENO(rkStep rk); //���������� ������������� 2D CWENO ��������

    //inspectors
    bool isDeleted() const;
    bool isLeaf() const;
    bool isGhost() const; //����������� �� � ghost-������
    bool hasChildren() const; //���� �� ����
    bool hasGrandChildren() const; //���� �� �����
    bool hasGrandChildren(Neighbour n) const; //���� �� ����� � ������������ �������
    bool hasNeighbour(Neighbour n) const; //���� �� ����� �� �����������
    bool hasNeighbour12(Neighbour12 n12) const;
    bool hasEdge(Edge etype) const; //���� �� �����
    int neighbours12Num() const; //����� �������12
    bool isNodeInSubstencil(Quadrant q, const TreeNode& rnnode) const; //�������� �� ������ � ��������� ��� �������� �������

    //other
    static TreeNode& nodeRef(const NodeTag& tag); //������ �� ���� �� ���� (static - ����� ������� ��� ���� ���)
    double magGradRho() const; //��������� �������� ���������
    ConservativeVector evalPolynomialAt(Point p, rkStep rk = 0); //������������������ ��������� CWENO �������� TODO: ������� const
    
    //output
    std::string dump() const; //���� ���� � ������
    std::string dumpNeighbourVector(Neighbour n) const; //���� ������ � ���� �������
    std::string dumpNeighbour12Vector(Neighbour12 n12) const; //���� ������12 � ���� �������
    std::string dumpEdgeVector(Edge etype) const; //���� ����� ������ � ���� �������
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
