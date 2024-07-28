#ifndef qtCWENO_QuadTree_H //include guard
#define qtCWENO_QuadTree_H

#include "main.h"
#include "TreeNode.h"

class QuadTree
{
    const int FIRST_LEVEL = 1; //������ define (static �� ���������?)
    const treeNodeId FIRST_ID = 0;
    quadTreeId id_;
    bool is_ghost = false; //�������� �� ������ ghost'��
    int depth_ = 0;

    std::vector<std::vector<TreeNode>> nodes; //����������� ��������� ������
    std::vector<std::vector<treeNodeId>> vacant_node_ids; //����������� ������ ��������� ���� � nodes
    std::vector<size_t> active_nodes_num; //����� �������� (�����������) ��� �� �������
    std::vector<size_t> leaf_nodes_num; //����� �������� ��� �� �������
    std::vector<CellData> data_; //���������� ������ ���������� ������ ��� �������
    std::vector<cellDataId> vacant_data_ids; //������ ��������� ���� � data

public:
    explicit QuadTree(quadTreeId _id); //����������� ������ � ����� �����

    //accessors
    treeNodeId id() const;
    int depth() const;
    TreeNode& rootRef(); //������ �� �������� ����
    const TreeNode& rootRefConst() const; //������ �� (const) �������� ����
    TreeNode& nodeRef(int depth, treeNodeId id); //������ �� ����
    CellData& dataRef(cellDataId id); //������ �� CellData �� id
    CellData data(cellDataId id) const; //����� CellData �� id
    TreeNode& getNodeByCoords(Point p) const; //������ �� ���� �� �����������

    //mutators
    void initNewLevel(); //������������� ������ ������ ������ (cross-check with coarsenTreeNode)
    void deleteLevelIfEmpty(int depth); //�������� ������, ���� �� ��� �� �������� �����
    void vacateNodeGroup(int depth, treeNodeId id); //������� ������ ��� ���������
    void vacateData(cellDataId did); //������� ������ ������ ���������
    void incrementCounterNodes(int dpth, int amount); //��������� �������� �������� ���
    void incrementCounterLeaves(int dpth, int amount); //��������� �������� �������

    //inspectors
    bool isGhost() const; //�������� �� ������ ghost'��
    bool isGhostCorner() const; //�������� �� ������ ������� ghost'��
    //static bool isTreeGhost(quadTreeId id); //�������� �� ������ ghost'�� (�� id)

    //other
    CellBox generateBox(const CellBox& global_box) const; //���������� bounding box ��� ������ �� ���������� ��������� ����� � ������� ������
    quadTreeId calcNeighbourTreeId(Neighbour n) const; //���������� id ��������� ������ �� id ��������
    quadTreeId calcNeighbour12TreeId(Neighbour12 n12) const; //���������� id ��������� ������ �� id ��������
    treeNodeId getVacantNodeId(int depth); //��������� ������ ��������� ������ ��� �������� ����� � ������� nodes[depth]
    cellDataId getVacantDataId(); //��������� ������ ��������� ������ ��� �������� ����� � ������� data


    //output
    const std::string dump() const; //���� ������ � ������

    friend class QuadTreeForest;
    friend class TreeNode;

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
