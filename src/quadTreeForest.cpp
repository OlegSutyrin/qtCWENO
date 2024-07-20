#include "main.h"
#include "output.h"
#include "CellBox.h"
#include "NodeEdge.h"
#include "TreeNode.h"
#include "QuadTree.h"
#include "QuadTreeForest.h"

//accessors -----------------------
QuadTree& QuadTreeForest::treeRef(quadTreeId id) //ссылка на дерево по id
{
    try {
        if (id < trees.size()) //есть такое дерево
        {
            return trees[id];
        }
        else
        {
            throw std::invalid_argument("invalid tree id");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "forest.treeRef() error: " << e.what() << ", tree id: " << id << endl;
        return trees[0]; //to suppress warning
    }
}
const QuadTree& QuadTreeForest::treeRefConst(quadTreeId id) const //const версия
{
    try {
        if (id < trees.size()) //есть такое дерево
        {
            return trees[id];
        }
        else
        {
            throw std::invalid_argument("invalid tree id");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "forest.treeRef() error: " << e.what() << ", tree id: " << id << endl;
        return trees[0]; //to suppress warning
    }
}
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
            throw std::invalid_argument("no tree with such point inside");
        }
        else
        {
            throw std::invalid_argument("point outise global box");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "forest.getTreeByCoords() error: " << e.what() << ", point: " << p << endl;
        return trees[0]; //to suppress warning
    }
}

//mutators -----------------------
void QuadTreeForest::initialize() //выделение памяти под деревья
{ 
    trees.reserve(config.Nx * config.Ny); //массив деревьев
    forest.toRefine.resize(config.max_depth + 1); //уровни списка ячеек для дробления (+1 для нулевого уровня) (resize увеличивает вектор и инициализирует элементы по умлочанию)

}
void QuadTreeForest::addTree(QuadTree tree) { trees.push_back(tree); } //добавление дерева в список
void QuadTreeForest::addNodeToRefine(const NodeTag& t) { forest.toRefine[t.depth()].push_back(t); } //добавление ноды в список на дробление

//inspectors -----------------------
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

//other -----------------------
const int INDEX_MACH = TECPLOT_FIELDS_NUMBER - 3;
const int INDEX_LEVEL = TECPLOT_FIELDS_NUMBER - 2;
const int INDEX_MAGGRADRHO = TECPLOT_FIELDS_NUMBER - 1;
dataExtrema QuadTreeForest::getExtrema() //сбор экстремумов всех величин для вывода в Tecplot
{
    dataExtrema ret{};
    std::array<double, TECPLOT_FIELDS_NUMBER> mins;
    std::array<double, TECPLOT_FIELDS_NUMBER> maxs;
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

void QuadTreeForest::meshApplyRefineList() //дробление ячеек из списка toRefine и балансировка
{
    for (auto depth = 1; depth < config.max_depth; depth++)
    {
        while (!forest.toRefine[depth].empty())
        {
            auto& rtag = forest.toRefine[depth].back();
            forest.toRefine[depth].pop_back();
            int err = TreeNode::nodeRef(rtag).refine();
            if (err != 0 && err != INT_ERROR_CODE_CANT_REFINE_ALREADY_REFINED)
                cout << "meshApplyRefineList refining error: " << err << ", tag: " << rtag << endl;
        }
    }
    return;
}

void QuadTreeForest::meshRefineInitial() //начальное дробление сетки
{
    if (config.meshRefineAll)
    {
        for (auto depth = 1; depth < config.max_depth; depth++)
        {
            cout << "Refining all, level " << depth << ": " << forest.leavesNumber();
            for (auto& rtree : trees) //проход по всем деревьям по ссылке
            {
                if (rtree.isGhost()) //пропуск ghost-деревьев
                    continue;
                for (auto& rnode : rtree.nodes[depth])
                    rnode.refine();
            }
            cout << " ---> " << forest.leavesNumber() << " leaves" << endl;
        }
        return;
    }
    //дробление ячеек вплоть до последнего слоя
    for (auto depth = 1; depth < config.max_depth; depth++)
    {
        cout << "Refining level " << depth << ": " << forest.leavesNumber();
        for (auto& tree : trees) //проход по всем деревьям по ссылке
        {
            if (tree.isGhost()) //пропуск ghost-деревьев
                continue;
            if (depth <= tree.depth()) //есть ячейки этого уровня
            {
                for (auto& rnode : tree.nodes[depth]) //проход по всем нодам на уровне
                {
                    bool to_refine = false;

                    //около УВ
                    if (rnode.box().intersectLineStraight(config.shock_position_x, Orientation::vertical) ||
                        rnode.box().bottom_left().isCloseToStraightLine(config.shock_position_x, Orientation::vertical) ||
                        rnode.box().bottom_right().isCloseToStraightLine(config.shock_position_x, Orientation::vertical))
                    {
                        to_refine = true;
                    }

                    //около границ слоя
                    if (config.problem == "layer")
                    {
                        //правая граница
                        if (rnode.box().top_left().y >= config.layer_bottom - config.meshInitialRefinePadding &&
                            rnode.box().bottom_left().y <= config.layer_top + config.meshInitialRefinePadding &&
                            (rnode.box().intersectLineStraight(config.layer_right, Orientation::vertical) ||
                                rnode.box().bottom_left().isCloseToStraightLine(config.layer_right, Orientation::vertical) ||
                                rnode.box().top_right().isCloseToStraightLine(config.layer_right, Orientation::vertical)))
                        {
                            to_refine = true;
                        }

                        //нижняя и верхняя границы
                        if (rnode.box().bottom_left().x <= config.layer_right + config.meshInitialRefinePadding &&
                            (rnode.box().intersectLineStraight(config.layer_bottom, Orientation::horizontal) ||
                                rnode.box().bottom_left().isCloseToStraightLine(config.layer_bottom, Orientation::horizontal) ||
                                rnode.box().top_right().isCloseToStraightLine(config.layer_bottom, Orientation::horizontal) ||
                                rnode.box().intersectLineStraight(config.layer_top, Orientation::horizontal) ||
                                rnode.box().bottom_left().isCloseToStraightLine(config.layer_top, Orientation::horizontal) ||
                                rnode.box().top_right().isCloseToStraightLine(config.layer_top, Orientation::horizontal)))
                        {
                            to_refine = true;
                        }
                    }

                    //около границы пузыря
                    if (config.problem == "bubble")
                    {
                        if (rnode.box().intersectLineEllipse(config.bubble_axle_x, config.bubble_axle_y))
                            to_refine = true;
                    }

                    //пометка к разделению
                    if (to_refine)
                        rnode.markToRefine();
                }
            }
        }
        //дробление и балансировка
        meshApplyRefineList();
        cout << " ---> " << forest.leavesNumber() << " leaves" << endl;
    }
    //сomputeQuadraturePoints(); //расчет точек квадратуры в ребрах
    //updateEigenObjects(); //перерасчет матриц Eigen в ячейках
    return;
}

void QuadTreeForest::meshCoarsenInitial() //начальное склеивание сетки (для теста)
{
    //склеивание от мелких к крупным
    for (auto depth = config.max_depth - 1; depth > 0; depth--)
    {
        cout << "Coarsening level " << depth << ": " << forest.leavesNumber();
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //пропуск ghost-деревьев
                continue;
            if (depth <= rtree.depth())
            {
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (!rnode.isDeleted() && rnode.hasChildren())
                    {
                        //около границы пузыря
                        if (config.problem == "bubble")
                        {
                            if (!rnode.box().intersectLineEllipse(config.bubble_axle_x, config.bubble_axle_y))
                                rnode.coarsen();
                        }
                    }
                }
            }
        }
        cout << " ---> " << forest.leavesNumber() << " leaves" << endl;
    }
}

void QuadTreeForest::meshUpdate() //обновление сетки
{
    //дробление от крупных к мелким
    for (auto depth = 1; depth < config.max_depth; depth++)
    {
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //пропуск ghost-деревьев
                continue;
            if (depth <= rtree.depth()) //есть ячейки этого уровня
            {
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (!rnode.isDeleted() && rnode.magGradRho() > globals.refineLvls[depth])
                        rnode.markToRefine();
                }
            }
        }
        //дробление и балансировка
        meshApplyRefineList();
    }

    //склеивание от мелких к крупным
    for (auto depth = config.max_depth - 1; depth > 0; depth--)
    {
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //пропуск ghost-деревьев
                continue;
            if (depth <= rtree.depth()) //есть ячейки этого уровня
            {
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (!rnode.isDeleted() && rnode.hasChildren() && !rnode.hasGrandChildren())
                    {
                        double maxGradRho = 0;
                        for (auto q : Quadrants)
                            maxGradRho = std::max(maxGradRho, rnode.childRef(q).magGradRho());
                        if (maxGradRho < globals.refineLvls[depth])
                            rnode.coarsen();
                    }
                }
            }
        }
    }
    //сomputeQuadraturePoints(); //перерасчет точек квадратуры в ребрах
    //updateEigenObjects(); //перерасчет матриц Eigen в ячейках
    return;
}

void QuadTreeForest::initialCondition() //начальные условия
{
    double double_gamma = config.gamma;
    double double_shock_Mach_number = config.Mach;
    double double_Atwood_number = config.Atwood;
    double omega = (1.0 + double_Atwood_number) / (1.0 - double_Atwood_number);
    double p0 = 1.0;
    double r0 = 1.0;
    double u0 = double_shock_Mach_number * sqrt(double_gamma); //скорость исходного скачка
    //параметры газа за исходным скачком
    double p1 = p0 * ((1.0 - double_gamma) / (double_gamma + 1.0) + 2.0 * double_gamma / (double_gamma + 1.0) * double_shock_Mach_number * double_shock_Mach_number);
    double r1 = r0 / ((double_gamma - 1.0) / (double_gamma + 1.0) + 2.0 / (double_gamma + 1.0) / double_shock_Mach_number / double_shock_Mach_number);
    double u1 = u0 * ((double_gamma - 1.0) / (double_gamma + 1.0) + 2.0 / (double_gamma + 1.0) / double_shock_Mach_number / double_shock_Mach_number);

    double double_shock_x_position = config.shock_position_x;
    double double_layer_x_right = config.layer_right;
    double double_layer_y_bottom = config.layer_bottom;
    double double_layer_y_top = config.layer_top;

    for (auto depth = 1; depth <= config.max_depth; depth++)
    {
        for (auto& rtree : trees)
        {
            if (depth <= rtree.depth()) //есть ячейки этого уровня
            {
                //проход по всем нодам на уровне
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (rnode.hasChildren()) //не-листья пропускаются
                        continue;
                    double x = rnode.box().center().x;
                    double y = rnode.box().center().y;
                    ConservativeVector Q = {};
                    if (config.problem == "layer")
                    {
                        if (x >= double_shock_x_position)
                        {
                            Q.set(Equation::density, r1);
                            Q.set(Equation::momentum_x, r1 * u1);
                            Q.set(Equation::momentum_y, 0); //v = 0
                            Q.set(Equation::energy, p1 / (config.gamma - 1.0) + 0.5 * r1 * (u1 * u1)); //v = 0
                        }
                        else
                        {
                            if (x <= double_layer_x_right && y >= double_layer_y_bottom && y <= double_layer_y_top)
                                Q.set(Equation::density, omega * r0);
                            else
                                Q.set(Equation::density, r0);
                            Q.set(Equation::momentum_x, Q(Equation::density) * u0);
                            Q.set(Equation::momentum_y, 0);  //v = 0
                            Q.set(Equation::energy, p0 / (config.gamma - 1.0) + 0.5 * Q(Equation::density) * (u0 * u0)); //v = 0
                        }
                    }
                    else if (config.problem == "bubble")
                    {
                        if (x <= double_shock_x_position)
                        {
                            Q.set(Equation::density, r1);
                            Q.set(Equation::momentum_x, r1 * (u0 - u1));
                            Q.set(Equation::momentum_y, 0);  //v = 0
                            Q.set(Equation::energy, p1 / (config.gamma - 1.0) + 0.5 * r1 * (u0 - u1) * (u0 - u1)); //v = 0
                        }
                        else
                        {
                            if (x * x / config.bubble_axle_x / config.bubble_axle_x + y * y / config.bubble_axle_y / config.bubble_axle_y <= 1.0)
                                Q.set(Equation::density, omega * r0);
                            else
                                Q.set(Equation::density, r0);
                            Q.set(Equation::momentum_x, 0); //u = 0
                            Q.set(Equation::momentum_y, 0); //v = 0
                            Q.set(Equation::energy, p0 / (config.gamma - 1.0)); //u,v = 0
                        }
                    }

                    CellData& rd = rnode.dataRef();
                    //умножение на r для осесимметричных координат
                    if (config.coord_type == CoordType::axisymmetric)
                    {
                        rd.setY(y);
                        for (auto eq : Equations)
                        {
                            Q.set(eq, Q(eq) * y);
                        }
                    }
                    //cout << "(" << x << ", " << y << "): " << Q << endl;
                    rd.setQ0(Q); //запись в CellData
                }
            }
        }
    }
    return;
}

//output -----------------------
void QuadTreeForest::exportScatter(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"X\", \"Y\", \"Level\", \"Size\", \"rho\", \"u\", \"v\", \"p\", \"tree\", \"magGradRho\"" << endl;
    file_output << "ZONE T = \"Forest scatter\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& rtree : trees)
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

void QuadTreeForest::exportNeighbours(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;
    file_output << "ZONE T = \"Forest neighbours\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& rtree : trees)
    {
        for (auto& rnodes_level : rtree.nodes)
        {
            for (auto& rnode : rnodes_level)
            {
                if (!rnode.isDeleted() && rnode.isLeaf())
                {
                    for (auto n : Neighbours)
                    {
                        if (rnode.hasNeighbour(n))
                            file_output << rnode.dumpNeighbourVector(n);
                    }
                }
            }
        }
    }
    file_output.close();
    return;
}
void QuadTreeForest::exportNeighbours12(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;
    file_output << "ZONE T = \"Forest neighbours12\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& rtree : trees)
    {
        for (auto& rnodes_level : rtree.nodes)
        {
            for (auto& rnode : rnodes_level)
            {
                if (!rnode.isDeleted() && rnode.isLeaf())
                {
                    for (auto n12 : Neighbours12)
                    {
                        if (rnode.hasNeighbour12(n12))
                            file_output << rnode.dumpNeighbour12Vector(n12);
                    }
                }
            }
        }
    }
    file_output.close();
    return;

}




