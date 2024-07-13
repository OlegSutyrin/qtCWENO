#ifndef qtCWENO_quadTreeForest_H //include guard
#define qtCWENO_quadTreeForest_H

#include <vector>
#include <fstream>

#include "quadTree.h"
#include "nodeEdge.h"

class quadTreeForest {
    std::vector<quadTree> trees; //�������
    std::vector<std::vector<nodeTag>> toRefine; //����������� ������ ���, ���������� ������������

public: //������� ������� ���������� ��� ��������� ������ �� ������� � ������� (TODO: �����������, ��� ������� - � ����� �� - �� ����������)
    //std::vector<nodeEdge> edges; //������� �����
    //std::vector<nodeEdgeId> vacant_edge_ids; //������ ��������� ���� � edges

    void initialize(); //��������� ������ ��� ������� � ������
    void addTree(quadTree tree); //���������� ������ � ������
    quadTree& getTreeByCoords(point p); //����� ������ �� ����������� �����
    //edgeId getVacantEdgeId(); //��������� ������ ��������� ������ ��� ������� ����� � ������� edges
    //edgeId addEdge(nodeEdge edge); //�������� ����� � ������
    //edgeId addEdgeUnique(nodeEdge edge); //�������� ����� � ������ � ��������� ������������
    //int updateEdge(edgeId eid, nodeTag n1, nodeTag n2); //���������� ������ �����
    //int removeEdge(edgeId eid); //�������� ����� �� ������
    dataExtrema getExtrema(); //���� ����������� ���� ������� ��� ������ � Tecplot
    size_t activeNodesNumber(); //������� �������� (�����������) ���
    size_t leavesNumber(); //������� �������

    void exportForestScatter(std::string filename); //����� � ���� (Tecplot ASCII scatter)

    friend class quadTree;
    friend class treeNode;


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