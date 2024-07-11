#ifndef qtCWENO_quadTree_H //include guard
#define qtCWENO_quadTree_H

#include "main.h"
#include "treeNode.h"

class quadTree
{
    const int FIRST_LEVEL = 1;
    const treeNodeId FIRST_ID = 0;
    quadTreeId id;
    bool is_ghost = false; //�������� �� ������ ghost'��
    int depth = 0;
    std::vector<std::vector<treeNode>> nodes; //����������� ��������� ������
    std::vector<std::vector<treeNodeId>> vacant_node_ids; //����������� ������ ��������� ���� � nodes
    std::vector<size_t> active_nodes_num; //����� �������� (�����������) ��� �� �������
    std::vector<size_t> leaf_nodes_num; //����� �������� ��� �� �������
    std::vector<cellData> data; //���������� ������ ���������� ������ ��� �������
    std::vector<cellDataId> vacant_data_ids; //������ ��������� ���� � data

public:
    quadTree(quadTreeId _id); //����������� ������ � ����� �����
    void initNewLevel(); //������������� ������ ������ ������ (cross-check with coarsenTreeNode)
    cellBox generateBox(cellBox global_box) const; //���������� bounding box ��� ������ �� ���������� ��������� ����� � ������� ������
    quadTreeId calcNeighbourTreeId(Neighbour Neighbour) const; //���������� id ��������� ������ �� id ��������
    //quadTreeId getNeighbour12TreeId(Neighbour12 Neighbour); //���������� id ��������� ������ �� id ��������
    treeNodeId getVacantNodeId(int depth); //��������� ������ ��������� ������ ��� �������� ����� � ������� nodes[depth]
    cellDataId getVacantDataId(); //��������� ������ ��������� ������ ��� �������� ����� � ������� data
    bool isGhost() const; //�������� �� ������ ghost'��
    bool isGhostCorner() const; //�������� �� ������ ������� ghost'��
    static bool isTreeGhost(quadTreeId id); //�������� �� ������ ghost'�� (�� id)
    const std::string dump() const; //���� ������ � ������
};

#endif
