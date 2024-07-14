#ifndef qtCWENO_treeNode_H //include guard
#define qtCWENO_treeNode_H

#include "main.h"
#include "cellData.h"
#include "nodeEdge.h"
#include "nodeTag.h"

#include "Eigen\Dense"

//���� ������
class treeNode
{
    nodeTag tag_ = {}; //{} �������� ��������� ����������� nodeTag
    treeNodeId parentId_ = null; //����� �������� �� ���������� ������
    treeNodeId childrenId_ = null; //����� ������� child �� ��������� ������, ��������� ���� ����� �� ���
    bool has_grandchildren = false; //���� �� ����� (����� ��� ������������)
    bool is_leaf = true; //�������� �� ������ ������ (������ ��� ��� �������� ���������� ������) (���������, �� ����� ������, ��� ��������� ������� �����)
    bool is_deleted = false; //��� ������� ����� ����� �� ���������, �.�. �������� = ������� ����� �� ������ ���������
    cellDataId dataId_ = null; //����� ������ ��� ������ � ���������� ������� cell_data, ������ ��� �������
    cellBox box_ = {};
    std::array<nodeTag, NEIGHBOURS_NUM> neighbours = {}; //������
    std::array<nodeTag, MAX_NEIGHBOURS12_NUM> neighbours12 = {}; //������ � ������ ������������ � ���������� ���������� ����� ������
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
    treeNode() {}; //default constructor
    treeNode(nodeTag t, treeNodeId p, cellDataId d, cellBox b) : tag_(t), parentId_(p), dataId_(d), box_(b) {}; //����������� �� ��������� ������
    bool isDeleted() const;
    bool isLeaf() const;
    nodeTag tag() const; //��������� ���� 
    void setTag(nodeTag t); //������� ����
    cellDataId dataId() const;
    void setDataId(cellDataId id);
    cellBox box() const;
    void setBox(cellBox b);
    nodeTag getNeighbour(Neighbour n) const; 
    void setNeighbour(Neighbour n, nodeTag t);
    nodeTag getNeighbour12(Neighbour n12) const; 
    void setNeighbour12(Neighbour12 n12, nodeTag t);
    void setData(cellData data); //������ ������

    cellData& dataRef() const; //������ �� ������
    cellData data() const; //����� ������
    treeNode& childRef(Quadrant q); //������ �� ������� �� ���������
    const treeNode& childRef(Quadrant q) const; //(const) ������ �� �������
    treeNode& getChildByCoords(point p); //������ �� ������� �� �����������

    
    
    static treeNode& nodeRef(nodeTag tag); //������ �� ���� �� ���� (static - ����� ������� ��� ���� ���)

    double magGradRho() const; //��������� �������� ���������

    std::string dump() const; //���� ���� � ������



    friend class quadTree;
    friend class quadTreeForest; //��� ���������� ������� TODO:�����������, ��� ����� ����������� ��������� ����� ��� ��������� ������������



    /*static treeNode& getNodeByCoords(point p); //������ �� ���� �� �����������
    treeNode& getChildByCoords(point p) const; //������ �� ������� �� �����������
    treeNode& getChild(Quadrant q) const; //������ �� ������� �� ���������
    nodeTag getNodeOrChildTag(int target_depth, Quadrant quadrant) const; //����� ��������� ���� ������� ������ ������ ������ (������� ������ �� ����� 1)
    //bool hasEdge(ushorty etype) const; //���� �� �����
    //nodeEdge& getEdge(ushorty etype) const; //������ �� �����
    double h() const; //����� ������� ������
    bool hasChildren() const; //���� �� ����
    bool hasGrandChildren() const; //���� �� �����
    void gatherDataFromChildren(); //���� ������ �� ����� ��� �����������
    void updateGrandChildren(); //���������� ������ � ������
    //bool hasNeighbour(Neighbour n) const; //���� �� ����� �� �����������
    bool hasNeighbour12(Neighbour12 n12) const;
    int neighbours12Num() const; //����� �������
    //bool isBoundary(); //�������� �� ���������
    //bool isCorner(); //�������� �� �������
    treeNode& getNeigbourOrSelf(Neighbour dir, int target_depth) const; //������ �� ������ (��� �� ����, ���� ��� ������)
    void updateNeighbour(Neighbour dir, nodeTag tag1, nodeTag tag2); //�������� ������ � ������ ����� ��� ���������
    void updateNeighbour(Neighbour dir, nodeTag tag); //�������� ������ � ������
    int markToRefine(); //������� ������ � ���������
    int refine(); //��������� ������
    int tryCoarsen(); //������� ������
    double magGradRho() const; //��������� �������� ���������
    bool isNeighbour12InSubstencil(Neighbour12 n, Quadrant q) const; //�������� �� ����� � ��������� ��� �������� �������
    void updateEigenObjects(); //�������� ��� ���������� Eigen ������ � �.�. ����� ��������� �����
    void calcPolynomialCWENO(rkStep rk); //���������� ������������� 2D CWENO ��������
    cellData evalPolynomialAt(point p, rkStep rk = 0); //������������������ ��������� �������� 
    std::string dump() const; //���� ���� � ������
    std::string dumpNeighbourVector(Neighbour n) const; //���� ������ � ���� �������
    std::string dumpNeighbour12Vector(Neighbour12 n) const;*/
};


#endif
