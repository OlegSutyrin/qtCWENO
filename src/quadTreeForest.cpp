#include "main.h"
#include "output.h"
#include "CellBox.h"
#include "NodeEdge.h"
#include "TreeNode.h"
#include "QuadTree.h"
#include "QuadTreeForest.h"

const int ERROR_BAD_TREE_ID = -1;
const int ERROR_TREE_BAD_COORDS = -2;
const int ERROR_COULDNT_FIND_TREE = -3;
void QuadTreeForest::initialize() //выделение памяти под деревья
{ 
    trees.reserve(config.Nx * config.Ny); //массив деревьев
    forest.toRefine.resize(config.max_depth + 1); //уровни списка ячеек для дробления (+1 для нулевого уровня) (resize увеличивает вектор и инициализирует элементы по умлочанию)

}
void QuadTreeForest::addTree(QuadTree tree) { trees.push_back(tree); } //доабвление дерева в список
QuadTree& QuadTreeForest::getTreeByCoords(Point p) //поиск дерева по координатам точки
{
    try {
        if (config.global_box.isPointInside(p))
        {
            for (auto& rtree : trees)
            {
                if (rtree.root().box().isPointInside(p))
                    return rtree;
            }
            throw ERROR_COULDNT_FIND_TREE;
        }
        else
        {
            throw ERROR_TREE_BAD_COORDS;
        }
    }
    catch (int err_code)
    {
        cout << "getTreeByCoords error: " << err_code << " coords " << p << endl;
        return trees[0]; //to suppress warning
    }
}

QuadTree& QuadTreeForest::treeRef(quadTreeId id) //ссылка на дерево по id
{
    try {
        if (id < trees.size()) //есть такое дерево
        {
            return trees[id];
        }
        else
        {
            throw ERROR_BAD_TREE_ID;
        }
    }
    catch (int err_code)
    {
        cout << "forest.tree() error: " << err_code << "tree id: " << id << endl;
        return trees[0]; //to suppress warning
    }
}

const int INDEX_MACH = TECPLOT_FIELDS_NUMBER - 3;
const int INDEX_LEVEL = TECPLOT_FIELDS_NUMBER - 2;
const int INDEX_MAGGRADRHO = TECPLOT_FIELDS_NUMBER - 1;
dataExtrema QuadTreeForest::getExtrema() //сбор экстремумов всех величин для вывода в Tecplot
{
    dataExtrema ret;
    double mins[TECPLOT_FIELDS_NUMBER], maxs[TECPLOT_FIELDS_NUMBER];
    for (int i = 0; i < TECPLOT_FIELDS_NUMBER; i++)
    {
        mins[i] = std::numeric_limits<double>::max();
        maxs[i] = -std::numeric_limits<double>::max();
    }

    double tree_h = (config.global_box.bottom_right().x - config.global_box.bottom_left().x) / config.Nx;
    mins[0] = config.global_box.bottom_left().x + tree_h; //поправка на один слой ghost-деревьев
    maxs[0] = config.global_box.bottom_right().x - tree_h;
    mins[1] = config.global_box.bottom_left().y + tree_h;
    maxs[1] = config.global_box.top_left().y - tree_h;
    mins[INDEX_LEVEL] = 1;
    maxs[INDEX_LEVEL] = config.max_depth;
    for (auto& tree : trees)
    {
        if (tree.isGhost()) //пропуск ghost-деревьев
            continue;
        for (auto& nodes_level : tree.nodes)
        {
            for (auto& cnode : nodes_level)
            {
                if (!cnode.isDeleted() && cnode.isLeaf())
                {
                    CellData& d = tree.dataRef(cnode.dataId());
                    int index = 2 + static_cast<int>(Equation::density); double tmp = d.rho();
                    if (mins[index] > tmp) { mins[index] = tmp; }
                    if (maxs[index] < tmp) { maxs[index] = tmp; }
                    index = 2 + static_cast<int>(Equation::momentum_x); tmp = d.u();
                    if (mins[index] > tmp) { mins[index] = tmp; }
                    if (maxs[index] < tmp) { maxs[index] = tmp; }
                    index = 2 + static_cast<int>(Equation::momentum_y); tmp = d.v();
                    if (mins[index] > tmp) { mins[index] = tmp; }
                    if (maxs[index] < tmp) { maxs[index] = tmp; }
                    index = 2 + static_cast<int>(Equation::energy); tmp = d.p();
                    if (mins[index] > tmp) { mins[index] = tmp; }
                    if (maxs[index] < tmp) { maxs[index] = tmp; }
                    index = INDEX_MACH; tmp = d.Mach();
                    if (mins[index] > tmp) { mins[index] = tmp; }
                    if (maxs[index] < tmp) { maxs[index] = tmp; }
                    index = INDEX_MAGGRADRHO;
                    tmp = cnode.magGradRho();
                    if (mins[index] > tmp) { mins[index] = tmp; }
                    if (maxs[index] < tmp) { maxs[index] = tmp; }
                }
            }
        }
    }
    //очистка от NaN
    for (int i = 0; i < TECPLOT_FIELDS_NUMBER; i++)
    {
        if (isnan(mins[i]) || isnan(maxs[i]) || mins[i] > maxs[i])
        {
            mins[i] = -999.0;
            maxs[i] = 999.0;
        }
    }

    for (int i = 0; i < TECPLOT_FIELDS_NUMBER; i++)
    {
        ret.minima[i] = mins[i];
        ret.maxima[i] = maxs[i];
    }
    return ret;
}

size_t QuadTreeForest::activeNodesNumber() //подсчет активных (неудаленных) нод
{
    size_t ret = 0;
    for (auto& rtree : forest.trees)
    {
        if (!rtree.isGhost())
            ret += rtree.active_nodes_num[0]; //в [0] хранится суммарное число по всем уровням
    }
    return ret;
}

size_t QuadTreeForest::leavesNumber() //подсчет листьев
{
    size_t ret = 0;
    for (auto& rtree : forest.trees)
    {
        if (!rtree.isGhost())
            ret += rtree.leaf_nodes_num[0]; //в [0] хранится суммарное число по всем уровням
    }
    return ret;
}

void QuadTreeForest::exportForestScatter(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"X\", \"Y\", \"Level\", \"Size\", \"rho\", \"u\", \"v\", \"p\", \"tree\", \"magGradRho\"" << endl;
    file_output << "ZONE T = \"Forest scatter\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& rtree : forest.trees)
    {
        for (auto& rnodes_level : rtree.nodes)
        {
            for (auto& rnode : rnodes_level)
            {
                if (!rnode.isDeleted() && rnode.isLeaf())
                {
                    CellData& rd = rnode.dataRef();
                    file_output << rnode.box().center().x << " " << rnode.box().center().y << " " << rnode.tag().depth() << " " << rnode.box().size() << " "
                        << " " << NaNcleared(rd.rho()) << " " << NaNcleared(rd.u()) << " " << NaNcleared(rd.v()) << " " << NaNcleared(rd.p()) << " " << rtree.id() << " " << NaNcleared(rnode.magGradRho());
                    file_output << endl;
                }
            }
        }
    }
    file_output.close();
    return;
}
