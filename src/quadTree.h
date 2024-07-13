#ifndef qtCWENO_quadTree_H //include guard
#define qtCWENO_quadTree_H

#include "main.h"
#include "treeNode.h"

class quadTree
{
    const int FIRST_LEVEL = 1; //������ define (static �� ���������?)
    const treeNodeId FIRST_ID = 0;
    quadTreeId id_;
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
    treeNodeId id() const;
    const treeNode& root() const; //������ �� (������������) �������� ����
    bool isGhost() const; //�������� �� ������ ghost'��
    bool isGhostCorner() const; //�������� �� ������ ������� ghost'��
    static bool isTreeGhost(quadTreeId id); //�������� �� ������ ghost'�� (�� id)
    void initNewLevel(); //������������� ������ ������ ������ (cross-check with coarsenTreeNode)
    cellBox generateBox(cellBox global_box) const; //���������� bounding box ��� ������ �� ���������� ��������� ����� � ������� ������
    quadTreeId calcNeighbourTreeId(Neighbour Neighbour) const; //���������� id ��������� ������ �� id ��������
    //quadTreeId getNeighbour12TreeId(Neighbour12 Neighbour); //���������� id ��������� ������ �� id ��������
    treeNodeId getVacantNodeId(int depth); //��������� ������ ��������� ������ ��� �������� ����� � ������� nodes[depth]
    cellDataId getVacantDataId(); //��������� ������ ��������� ������ ��� �������� ����� � ������� data
    treeNode& getNodeByCoords(point p) const; //������ �� ���� �� �����������
    const std::string dump() const; //���� ������ � ������

    friend class quadTreeForest;
    friend class treeNode;

    //������ ����� �� �����
    /*template<typename Func>
    inline void forAllNodes(Func f, bool include_branches = SKIP_BRANCHES)
    {
        for (auto& nodes_level : nodes) //�� ������
        {
            for (auto& cnode : nodes_level) //�� ������
            {
                if (cnode.isDeleted() || !include_branches && !cnode.isLeaf())
                    continue;
                f(cnode);
            }
        }
    }*/
};

#endif
