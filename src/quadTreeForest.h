#ifndef qtCWENO_QuadTreeForest_H //include guard
#define qtCWENO_QuadTreeForest_H

#include <vector>
#include <fstream>

#include "QuadTree.h"
#include "NodeEdge.h"

class QuadTreeForest {
    std::vector<QuadTree> trees; //�������
    std::vector<std::vector<NodeTag>> toRefine; //����������� ������ ���, ���������� ������������

public: //������� ������� ���������� ��� ��������� ������ �� ������� � ������� (TODO: �����������, ��� ������� - � ����� �� - �� ����������)
    //std::vector<nodeEdge> edges; //������� �����
    //std::vector<nodeEdgeId> vacant_edge_ids; //������ ��������� ���� � edges

    void initialize(); //��������� ������ ��� ������� � ������
    void addTree(QuadTree tree); //���������� ������ � ������
    QuadTree& treeRef(quadTreeId id); //������ �� ������ �� id
    QuadTree& getTreeByCoords(Point p); //����� ������ �� ����������� �����
    //edgeId getVacantEdgeId(); //��������� ������ ��������� ������ ��� ������� ����� � ������� edges
    //edgeId addEdge(nodeEdge edge); //�������� ����� � ������
    //edgeId addEdgeUnique(nodeEdge edge); //�������� ����� � ������ � ��������� ������������
    //int updateEdge(edgeId eid, NodeTag n1, NodeTag n2); //���������� ������ �����
    //int removeEdge(edgeId eid); //�������� ����� �� ������
    dataExtrema getExtrema(); //���� ����������� ���� ������� ��� ������ � Tecplot
    size_t activeNodesNumber(); //������� �������� (�����������) ���
    size_t leavesNumber(); //������� �������

    void exportForestScatter(std::string filename); //����� � ���� (Tecplot ASCII scatter)

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

#endif