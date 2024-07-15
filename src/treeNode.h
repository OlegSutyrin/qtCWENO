#ifndef qtCWENO_TreeNode_H //include guard
#define qtCWENO_TreeNode_H

#include "main.h"
#include "CellData.h"
#include "NodeEdge.h"
#include "NodeTag.h"

#include "Eigen\Dense"

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
    //nodeEdgeId edges[DIRECTIONS_NUM * 2] = { null, null, null, null, null, null, null, null }; //��������� �� ������� �������, ������� � ����� �������� ������� �������
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
    NodeTag tag() const;
    cellDataId dataId() const;
    CellBox box() const;
    NodeTag getNeighbour(Neighbour n) const; 
    NodeTag getNeighbour12(Neighbour n12) const; 
    CellData& dataRef() const; //������ �� ������
    CellData data() const; //����� ������ (� ��������� ������ � �����)
    TreeNode& childRef(Quadrant q); //������ �� ������� �� ���������
    const TreeNode& childRef(Quadrant q) const; //const ������
    TreeNode& getChildByCoords(Point p); //������ �� ������� �� �����������

    //mutators
    void setTag(const NodeTag& t);
    void setDataId(cellDataId id);
    void setBox(const CellBox& b);
    void setNeighbour(Neighbour n, NodeTag t);
    void setNeighbour12(Neighbour12 n12, NodeTag t);
    void setData(const CellData& data); //������ ������

    //inspectors
    bool isDeleted() const;
    bool isLeaf() const;
    bool hasChildren() const; //���� �� ����
    bool hasGrandChildren() const; //���� �� �����

    //uncategorized
    static TreeNode& nodeRef(const NodeTag& tag); //������ �� ���� �� ���� (static - ����� ������� ��� ���� ���)
    double magGradRho() const; //��������� �������� ���������

    //output
    std::string dump() const; //���� ���� � ������

    friend class QuadTreeForest; //��� ���������� ������� TODO:�����������, ��� ����� ����������� ��������� ����� ��� ��������� ������������
    friend class QuadTree;


    /*static TreeNode& getNodeByCoords(point p); //������ �� ���� �� �����������
    NodeTag getNodeOrChildTag(int target_depth, Quadrant quadrant) const; //����� ��������� ���� ������� ������ ������ ������ (������� ������ �� ����� 1)
    bool hasEdge(ushorty etype) const; //���� �� �����
    nodeEdge& getEdge(ushorty etype) const; //������ �� �����
    double h() const; //����� ������� ������
    void gatherDataFromChildren(); //���� ������ �� ����� ��� �����������
    void updateGrandChildren(); //���������� ������ � ������
    //bool hasNeighbour(Neighbour n) const; //���� �� ����� �� �����������
    bool hasNeighbour12(Neighbour12 n12) const;
    int neighbours12Num() const; //����� �������
    //bool isBoundary(); //�������� �� ���������
    //bool isCorner(); //�������� �� �������
    TreeNode& getNeigbourOrSelf(Neighbour dir, int target_depth) const; //������ �� ������ (��� �� ����, ���� ��� ������)
    void updateNeighbour(Neighbour dir, NodeTag tag1, NodeTag tag2); //�������� ������ � ������ ����� ��� ���������
    void updateNeighbour(Neighbour dir, NodeTag tag); //�������� ������ � ������
    int markToRefine(); //������� ������ � ���������
    int refine(); //��������� ������
    int tryCoarsen(); //������� ������
    bool isNeighbour12InSubstencil(Neighbour12 n, Quadrant q) const; //�������� �� ����� � ��������� ��� �������� �������
    void updateEigenObjects(); //�������� ��� ���������� Eigen ������ � �.�. ����� ��������� �����
    void calcPolynomialCWENO(rkStep rk); //���������� ������������� 2D CWENO ��������
    CellData evalPolynomialAt(point p, rkStep rk = 0); //������������������ ��������� �������� 
    std::string dumpNeighbourVector(Neighbour n) const; //���� ������ � ���� �������
    std::string dumpNeighbour12Vector(Neighbour12 n) const;*/
};


#endif
