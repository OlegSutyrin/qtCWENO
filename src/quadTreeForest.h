#ifndef qtCWENO_QuadTreeForest_H //include guard
#define qtCWENO_QuadTreeForest_H

#include <vector>
#include <fstream>

#include "QuadTree.h"
#include "NodeEdge.h"

class QuadTreeForest {
    std::vector<QuadTree> trees; //�������
    std::vector<std::vector<NodeTag>> toRefine; //����������� ������ ���, ���������� ������������
    std::vector<NodeEdge> edges; //������� �����
    std::vector<nodeEdgeId> vacant_edge_ids; //������ ��������� ���� � edges

public:
    //accessors
    QuadTree& treeRef(quadTreeId id); //������ �� ������ �� id
    const QuadTree& treeRefConst(quadTreeId id) const; //const ������
    QuadTree& getTreeByCoords(Point p); //����� ������ �� ����������� �����
    NodeEdge& edgeRef(nodeEdgeId eid); //������ �� ����� �� id

    //mutators
    void initialize(); //��������� ������ ��� ������� � ������
    void addTree(QuadTree tree); //���������� ������ � ������
    void addNodeToRefine(const NodeTag& t); //���������� ���� � ������ �� ���������
    nodeEdgeId addEdge(NodeEdge edge); //�������� ����� � ������
    nodeEdgeId addEdgeUnique(NodeEdge edge); //�������� ����� � ������ � ��������� ������������
    void updateEdge(nodeEdgeId eid, NodeTag n1, NodeTag n2); //���������� ������ �����
    void removeEdge(nodeEdgeId eid); //�������� ����� �� ������

    //inspectors
    size_t activeNodesNumber() const; //������� �������� (�����������) ���
    size_t leavesNumber() const; //������� �������

    //other
    dataExtrema getExtrema(); //���� ����������� ���� ������� ��� ������ � Tecplot
    nodeEdgeId getVacantEdgeId(); //��������� ������ ��������� ������ ��� ������� ����� � ������� edges
    void meshApplyRefineList(); //��������� ����� �� ������ toRefine � ������������
    void meshRefineInitial(); //��������� ��������� �����
    void meshCoarsenInitial(); //��������� ���������� ����� (��� �����)
    void meshUpdate(); //���������� �����
    void computeQuadraturePoints(); //������ ����� ���������� �� ���� ������
    void updateEigenObjects(); //�������� ��� ���������� Eigen �������� ��� ���� �����
    void initialCondition(); //��������� �������
    double CFLTimestepSize(); //�������� ���� �� �������
    void putQn(rkStep rk); //�������� Qn[rk_order] -> Qn[0]
    void computePolynomialCoeffs(rkStep rk); //���������� ������������� CWENO �������� �� ���� �������
    void computeFluxesCWENO(rkStep rk); //������ ������� �� ���� ������
    void advanceTime(); //������ ��� �� �������: ���������� Qn[rk_order]
    void boundaryConditions(rkStep rk); //��������� �������
    void boundaryConditionsAll(); //��������� ������� ��� ���� rk

    //output
    void exportScatter(std::string filename); //����� � ���� (Tecplot ASCII scatter)
    void exportNeighbours(std::string filename); //����� ������� � ���� (Tecplot ASCII vectors)
    void exportNeighbours12(std::string filename); //����� �������12 � ���� (Tecplot ASCII vectors)
    void exportEdges(std::string filename); //����� ����� � ���� (Tecplot ASCII scatter and vectors)
    void exportNodeEdges(std::string filename); //����� ����� ����� � ���� (Tecplot ASCII vectors)

    friend class QuadTree;
    friend class TreeNode;

    //������ ����� �� ��������
    /*template<typename Func>
    inline void forAllTrees(Func f, bool include_ghosts = SKIP_GHOSTS)
    {
        for (auto& tree : trees) //�� ������
        {
            if (!include_ghosts && tree.isGhost()) //������� ghost-��������
                continue;
            f(tree);
        }
    }

    //������ ����� �� �����
    template<typename Func>
    inline void forAllNodes(Func f, bool include_ghosts = SKIP_GHOSTS, bool include_branches = SKIP_BRANCHES)
    {
        for (auto& tree : trees) //�� ������
        {
            if (!include_ghosts && tree.isGhost()) //������� ghost-��������
                continue;
            tree.forAllNodes(f, include_branches);
        }
    }*/
};

//������������ �����-�����
const double RKcoeff[RK_ORDER_MAX][RK_ORDER_MAX][RK_ORDER_MAX + 1] = { //[order-1][step][Qn], [..][..][RK_ORDER_MAX] - �����������
    {
        {1.0, 0.0, 0.0, 1.0}, //1 �������: explicit Euler
        {},
        {}
    },
    {
        {1.0, 0.0, 0.0, 1.0}, //2 Heun's method
        {0.5, 0.5, 0.0, 0.5}, //(Butcher tableau: a21 = 1 | b1 = 1/2, b2 = 1/2)
        {}
    },
    {
        {1.0, 0.0, 0.0, 1.0}, //3 Strong stability preserving (TVD) Runge-Kutta (SSPRK3)
        {0.75, 0.25, 0.0, 0.25},
        {1.0 / 3.0, 0.0, 2.0 / 3.0, 2.0 / 3.0} //(Butcher tableau: a21 = 1 | a31 = 1/4, a32 = 1/4 | b1 = 1/6, b2 = 1/6, b3 = 2/3)
    }
};


#endif