#include "main.h"
#include "cellBox.h"
#include "cellData.h"
#include "treeNode.h"
#include "quadTree.h"

#include <sstream>


quadTree::quadTree(quadTreeId _id)  //����������� ������ � ����� �����
{
    id = _id;
    //������� ������� �� ������������ ��� �������� ���������
    nodes.emplace_back();
    active_nodes_num.push_back(1); //�� ������� ������ �������� ����� ����� �������� ���
    leaf_nodes_num.push_back(1); //����� ����� �������
    vacant_node_ids.emplace_back();
    //������ ���� - ���� ����
    initNewLevel(); //��������� ����� ��� ������ ����
    nodes[FIRST_LEVEL].emplace_back(); //��������� ����� ��� ���� ���� (� ������������� ���������� �� ���������)
    treeNode& root = nodes[FIRST_LEVEL][FIRST_ID]; //�� ������
    root.setTag({ id, FIRST_LEVEL, FIRST_ID });
    cellBox root_box = generateBox(config.global_box);
    root.setBox(root_box);
    //root.is_leaf = true; (=true �� ���������)
    //root.is_deleted = false; (=false �� ���������)
    //root.parent = null; (=null �� ���������)
    active_nodes_num[FIRST_LEVEL] = 1; //������� �������� ���
    leaf_nodes_num[FIRST_LEVEL] = 1; //������� �������� ���

    //nodeTag ctag = { root.tree_id, root.depth, root.id }; //��� �������� ����
    //������ � �����
    for (auto n : Neighbours)
    {
        auto nid = calcNeighbourTreeId(n);
        if (nid == null)
            is_ghost = true;
        else
        {
            //������ � ������ (������ ��������� ������)
            root.setNeighbour(n, { nid, FIRST_LEVEL, FIRST_ID });
            //�������� ����� ����� ��������
            /*nodeEdge edge;
            switch (n)
            {
            case NEIGHBOUR_TOP:
                edge = { root.neighbours[n], ctag, ORIENTATION_HORIZONTAL };
                break;
            case NEIGHBOUR_RIGHT:
                edge = { ctag, root.neighbours[n], ORIENTATION_VERTICAL };
                break;
            case NEIGHBOUR_BOTTOM:
                edge = { ctag, root.neighbours[n], ORIENTATION_HORIZONTAL };
                break;
            case NEIGHBOUR_LEFT:
                edge = { root.neighbours[n], ctag, ORIENTATION_VERTICAL };
                break;
            }*/
            //���������� ����� � ������ � ��������� ������������
            //root.edges[2 * n] = forest.addEdgeUnique(edge); //��������� �� 2 ��-�� ��������� �����
            //root.edges[2 * n + 1] = null; //"������ ��������" ������������� �������
            //cout << "created edge " << forest.edges[root.edges[2 * n]] << endl;
        }
    }

    //������ � ������ ������������
    /*for (auto n : Neighbours12)
    {
        auto nid = getNeighbour12TreeId(n);
        root.neighbours12[n] = { nid, FIRST_LEVEL, FIRST_ID };
    }*/

    //��������� ������ ��� ���������� ������
    root.setDataId(getVacantDataId());
    if (config.coord_type == CoordType::axisymmetric) //�������� r � cellData ��� ��������������� ���������
    {
        data[root.dataId()].setY(root_box.center().y);
    }
}

void quadTree::initNewLevel()  //������������� ������ ������ ������ (cross-check with coarsenTreeNode)
{
    std::vector<treeNode> nodes_level;
    nodes.push_back(nodes_level); //����� ���� ���
    active_nodes_num.push_back(0); //����� ������ � ����� �������� ��� �� ����
    leaf_nodes_num.push_back(0); //����� ������ � ����� ������� �� ����
    std::vector<treeNodeId> vacant_node_ids_level;
    vacant_node_ids.push_back(vacant_node_ids_level); //����� ���� ������� � ��������� �����
    depth++; //������� ������ +1
    return;
}

cellBox quadTree::generateBox(cellBox global_box) const //���������� bounding box ��� ������ �� ���������� ��������� ����� � ������� ������
{
    double dx = (global_box.bottom_right().x - global_box.bottom_left().x) / config.Nx;
    double dy = (global_box.top_left().y - global_box.bottom_left().y) / config.Ny;
    size_t discrete_coord_x = (size_t)(id / config.Ny) % config.Nx;
    size_t discrete_coord_y = (size_t)id % config.Ny;
    return cellBox({ global_box.bottom_left().x + discrete_coord_x * dx,global_box.bottom_left().y + discrete_coord_y * dy }, 
        { global_box.bottom_left().x + (discrete_coord_x + 1) * dx,global_box.bottom_left().y + (discrete_coord_y + 1) * dy });
}

quadTreeId quadTree::calcNeighbourTreeId(Neighbour Neighbour) const //���������� id ��������� ������ �� id ��������
{
    quadTreeId ret = null;
    size_t discrete_coord_x = (size_t)(id / config.Ny) % config.Nx;
    size_t discrete_coord_y = id % config.Ny;
    switch (Neighbour)
    {
    case Neighbour::top:
        if (discrete_coord_y < config.Ny - 1)
            ret = id + 1;
        break;
    case Neighbour::right:
        if (discrete_coord_x < config.Nx - 1)
            ret = id + config.Ny;
        break;
    case Neighbour::bottom:
        if (discrete_coord_y > 0)
            ret = id - 1;
        break;
    case Neighbour::left:
        if (discrete_coord_x > 0)
            ret = id - config.Ny;
        break;
    }
    return ret;
}

cellDataId quadTree::getVacantDataId() //��������� ������ ��������� ������ ��� �������� ����� � ������� data
{
    cellDataId ret;
    //���� ��� ��������� ���� - ������ � ����� �������
    if (vacant_data_ids.empty())
    {
        ret = data.size(); //�� 1 ������ ������ ���������� ��������
        data.emplace_back(); //��������� ������ (������ ���������� ����� .size())
    }
    else
    {
        ret = vacant_data_ids.back(); //����� ��������� (�� ������� ������������) ��������� ������ � ������� data
        vacant_data_ids.pop_back(); //������ �� ������ ��������� �����
    }
    return ret;
}