#include "main.h"
#include "CellBox.h"
#include "CellData.h"
#include "TreeNode.h"
#include "QuadTree.h"

#include <sstream>


QuadTree::QuadTree(quadTreeId _id) //����������� ������ � ����� �����
{
    id_ = _id;
    //������� ������� �� ������������ ��� �������� ���������
    nodes.emplace_back();
    active_nodes_num.push_back(1); //�� ������� ������ �������� ����� ����� �������� ���
    leaf_nodes_num.push_back(1); //����� ����� �������
    vacant_node_ids.emplace_back();
    //������ ���� - ���� ����
    initNewLevel(); //��������� ����� ��� ������ ����
    nodes[FIRST_LEVEL].emplace_back(); //��������� ����� ��� ���� ���� (� ������������� ���������� �� ���������)
    TreeNode& root = nodes[FIRST_LEVEL][FIRST_ID]; //�� ������
    root.setTag(NodeTag(id_, FIRST_LEVEL, FIRST_ID));
    CellBox root_box = generateBox(config.global_box);
    root.setBox(root_box);
    //root.is_leaf = true; (=true �� ���������)
    //root.is_deleted = false; (=false �� ���������)
    //root.parent = null; (=null �� ���������)
    active_nodes_num[FIRST_LEVEL] = 1; //������� �������� ���
    leaf_nodes_num[FIRST_LEVEL] = 1; //������� �������� ���

    //NodeTag ctag = { root.tree_id, root.depth, root.id }; //��� �������� ����
    //������ � �����
    for (auto n : Neighbours)
    {
        auto nid = calcNeighbourTreeId(n);
        if (nid == null)
            is_ghost = true;
        else
        {
            //������ � ������ (������ ��������� ������)
            root.setNeighbour(n, NodeTag(nid, FIRST_LEVEL, FIRST_ID));
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
    for (auto n : Neighbours12)
    {
        auto nid = calcNeighbour12TreeId(n);
        root.setNeighbour12(n, NodeTag(nid, FIRST_LEVEL, FIRST_ID));
    }

    //��������� ������ ��� ���������� ������
    root.setDataId(getVacantDataId());
    if (config.coord_type == CoordType::axisymmetric) //�������� r � CellData ��� ��������������� ���������
    {
        data_[root.dataId()].setY(root_box.center().y);
    }
}

//accessors -----------------------
treeNodeId QuadTree::id() const { return id_; }
int QuadTree::depth() const { return depth_; }
const TreeNode& QuadTree::root() const { return nodes[FIRST_LEVEL][FIRST_ID]; } //������ �� (������������) �������� ����

TreeNode& QuadTree::nodeRef(int d, treeNodeId id)
{
    try {
        if (d <= depth_ && id < nodes[d].size()) //���� ����� ���� (�.�. ���������)
        {
            return nodes[d][id];
        }
        else
        {
            throw std::invalid_argument("no such node in the tree");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "tree.nodeRef() error: " << e.what() << "tag: (" << id_ << ", " << d << ", " << id << ")" << endl;
        return nodes[FIRST_LEVEL][FIRST_ID]; //to suppress warning
    }
}

CellData& QuadTree::dataRef(cellDataId id) //������ �� CellData �� id
{
    try {
        if (id < data_.size()) //���� ����� ������ (�.�. ���������)
        {
            return data_[id];
        }
        else
        {
            throw std::invalid_argument("invalid data id");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "tree.dataRef() error: " << e.what() << ", data id: " << id << endl;
        return data_[0]; //to suppress warning
    }
}
CellData QuadTree::data(cellDataId id) const //����� CellData �� id
{
    try {
        if (id < data_.size()) //���� ����� ������ (�.�. ���������)
        {
            return data_[id];
        }
        else
        {
            throw std::invalid_argument("invalid data id");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "tree.data() error: " << e.what() << ", data id: " << id << endl;
        return data_[0]; //to suppress warning
    }
}

TreeNode& QuadTree::getNodeByCoords(Point p) const //������ �� ���� �� �����������
{
    auto& tree = forest.getTreeByCoords(p);
    return tree.nodes[FIRST_LEVEL][FIRST_ID].getChildOrSelfByCoords(p);
}

//mutators -----------------------
void QuadTree::initNewLevel() //������������� ������ ������ ������
{
    std::vector<TreeNode> nodes_level;
    nodes.push_back(nodes_level); //����� ���� ���
    active_nodes_num.push_back(0); //����� ������ � ����� �������� ��� �� ����
    leaf_nodes_num.push_back(0); //����� ������ � ����� ������� �� ����
    std::vector<treeNodeId> vacant_node_ids_level;
    vacant_node_ids.push_back(vacant_node_ids_level); //����� ���� ������� � ��������� �����
    depth_++; //������� ������ +1
    return;
}
void QuadTree::deleteLevelIfEmpty(int depth)
{
    if (depth != depth_) //������� ������� �� ��������� ����
        return;
    if (active_nodes_num[depth] == 0)
    {
        nodes[depth].clear(); //�������� (���������) ���
        nodes[depth].shrink_to_fit(); //������������ ������
        nodes.pop_back(); //�������� ����
        active_nodes_num.pop_back(); //�������� ������ � ����� �������� ��� �� ����
        leaf_nodes_num.pop_back(); //�������� ������ � ����� ������� �� ����
        vacant_node_ids[depth].clear(); //�������� ���� ������� � ��������� ����� �� ����
        vacant_node_ids[depth].shrink_to_fit(); //������������ ������ 
        vacant_node_ids.pop_back(); //�������� ���� 
        depth--; //������� ������ -1
    }
}

void QuadTree::vacateNodeGroup(int depth, treeNodeId id) { vacant_node_ids[depth].push_back(id); } //������� ������ ��� ���������
void QuadTree::vacateData(cellDataId did) { vacant_data_ids.push_back(did); } //������� ������ ������ ���������
void QuadTree::incrementNodesCounter(int dpth, int amount) //��������� �������� �������� ���
{
    active_nodes_num[dpth] += amount;
    active_nodes_num[0] += amount; //����� ����� ���
}
void QuadTree::incrementLeavesCounter(int dpth, int amount) //��������� �������� �������
{
    leaf_nodes_num[dpth] += amount;
    leaf_nodes_num[0] += amount; //����� ����� �������
}

//inspectors -----------------------
bool QuadTree::isGhost() const { return is_ghost; }

//other -----------------------
CellBox QuadTree::generateBox(const CellBox& global_box) const //���������� bounding box ��� ������ �� ���������� ��������� ����� � ������� ������
{
    double dx = (global_box.bottom_right().x - global_box.bottom_left().x) / config.Nx;
    double dy = (global_box.top_left().y - global_box.bottom_left().y) / config.Ny;
    size_t discrete_coord_x = (size_t)(id_ / config.Ny) % config.Nx;
    size_t discrete_coord_y = (size_t)id_ % config.Ny;
    return CellBox({ global_box.bottom_left().x + discrete_coord_x * dx,global_box.bottom_left().y + discrete_coord_y * dy }, 
        { global_box.bottom_left().x + (discrete_coord_x + 1) * dx,global_box.bottom_left().y + (discrete_coord_y + 1) * dy });
}

quadTreeId QuadTree::calcNeighbourTreeId(Neighbour n) const //���������� id ��������� ������ �� id ��������
{
    quadTreeId ret = null;
    size_t discrete_coord_x = (size_t)(id_ / config.Ny) % config.Nx;
    size_t discrete_coord_y = id_ % config.Ny;
    switch (n)
    {
    case Neighbour::top:
        if (discrete_coord_y < config.Ny - 1)
            ret = id_ + 1;
        break;
    case Neighbour::right:
        if (discrete_coord_x < config.Nx - 1)
            ret = id_ + config.Ny;
        break;
    case Neighbour::bottom:
        if (discrete_coord_y > 0)
            ret = id_ - 1;
        break;
    case Neighbour::left:
        if (discrete_coord_x > 0)
            ret = id_ - config.Ny;
        break;
    }
    return ret;
}
quadTreeId QuadTree::calcNeighbour12TreeId(Neighbour12 n12) const
{
    quadTreeId ret = null;
    size_t discrete_coord_x = (size_t)(id_ / config.Ny) % config.Nx;
    size_t discrete_coord_y = id_ % config.Ny;
    switch (n12)
    {
    case Neighbour12::top1: //������������ ������� ������������ ��� top1 = <tag>, top2 = null
        if (discrete_coord_y < config.Ny - 1)
            ret = id_ + 1;
        break;
    case Neighbour12::top_right:
        if (discrete_coord_y < config.Ny - 1 && discrete_coord_x < config.Nx - 1)
            ret = id_ + 1 + config.Ny;
        break;
    case Neighbour12::right1: //������������ ������
        if (discrete_coord_x < config.Nx - 1)
            ret = id_ + config.Ny;
        break;
    case Neighbour12::bottom_right:
        if (discrete_coord_y > 0 && discrete_coord_x < config.Nx - 1)
            ret = id_ - 1 + config.Ny;
        break;
    case Neighbour12::bottom1: //������������ ������
        if (discrete_coord_y > 0)
            ret = id_ - 1;
        break;
    case Neighbour12::bottom_left:
        if (discrete_coord_y > 0 && discrete_coord_x > 0)
            ret = id_ - 1 - config.Ny;
        break;
    case Neighbour12::left1: //������������ �����
        if (discrete_coord_x > 0)
            ret = id_ - config.Ny;
        break;
    case Neighbour12::top_left:
        if (discrete_coord_y < config.Ny - 1 && discrete_coord_x > 0)
            ret = id_ + 1 - config.Ny;
        break;
    default:
        ret = null;
    }
    return ret;
}

treeNodeId QuadTree::getVacantNodeId(int depth) //��������� ������ ��������� ������ ��� �������� ����� � ������� nodes[depth]
{
    treeNodeId ret;
    //���� ��� ��������� ���� - ������ � ����� �������
    if (vacant_node_ids[depth].empty())
    {
        ret = nodes[depth].size(); //�� 1 ������ ������ ���������� ��������
        for (auto q : Quadrants) //����� ������� ����� �������?
            nodes[depth].emplace_back(); //��������� ������ (� �������������?), ������ ���������� ����� .size() ����
    }
    else
    {
        ret = vacant_node_ids[depth].back(); //����� ��������� (�� ������� ������������) ��������� ������ � ������� nodes[depth]
        vacant_node_ids[depth].pop_back(); //������ ���� ����� �� ������ ��������� �����
    }
    return ret;
}

cellDataId QuadTree::getVacantDataId() //��������� ������ ��������� ������ ��� �������� ����� � ������� data
{
    cellDataId ret;
    //���� ��� ��������� ���� - ������ � ����� �������
    if (vacant_data_ids.empty())
    {
        ret = data_.size(); //�� 1 ������ ������ ���������� ��������
        data_.emplace_back(); //��������� ������ (������ ���������� ����� .size())
    }
    else
    {
        ret = vacant_data_ids.back(); //����� ��������� (�� ������� ������������) ��������� ������ � ������� data
        vacant_data_ids.pop_back(); //������ �� ������ ��������� �����
    }
    return ret;
}

//output -----------------------
const std::string QuadTree::dump() const //���� ������ � ������ (��� ���)
{
    std::stringstream ret;
    ret << "Tree " << id() << ", " << "depth " << depth() << " ( ";
    for (auto& n : active_nodes_num)
    {
        ret << n << " ";
    }
    ret << "):" << endl;
    return ret.str();
}
