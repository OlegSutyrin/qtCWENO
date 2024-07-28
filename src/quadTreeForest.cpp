#include "main.h"
#include "output.h"
#include "CellBox.h"
#include "NodeEdge.h"
#include "TreeNode.h"
#include "QuadTree.h"
#include "QuadTreeForest.h"

//accessors -----------------------
QuadTree& QuadTreeForest::treeRef(quadTreeId id) //������ �� ������ �� id
{
    try {
        if (id < trees.size()) //���� ����� ������
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
const QuadTree& QuadTreeForest::treeRefConst(quadTreeId id) const //const ������
{
    try {
        if (id < trees.size()) //���� ����� ������
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
QuadTree& QuadTreeForest::getTreeByCoords(Point p) //����� ������ �� ����������� �����
{
    try {
        if (config.global_box.isPointInside(p))
        {
            for (auto& rtree : trees)
            {
                if (rtree.rootRefConst().box().isPointInside(p))
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
NodeEdge& QuadTreeForest::edgeRef(nodeEdgeId eid) //������ �� ����� �� id
{ 
    try {
        if (edges.size() > eid && !edges[eid].isDeleted())
        {
            return edges[eid];
        }
        else
        {
            throw std::invalid_argument("edge is deleted or doesn't exist");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "forest.edgeRef() error: " << e.what() << ", eid: " << eid << endl;
        return edges[0]; //to suppress warning
    }
}

//mutators -----------------------
void QuadTreeForest::initialize() //��������� ������ ��� �������
{ 
    trees.reserve(config.Nx * config.Ny); //������ ��������
    forest.toRefine.resize(config.max_depth + 1); //������ ������ ����� ��� ��������� (+1 ��� �������� ������) (resize ����������� ������ � �������������� �������� �� ���������)

}
void QuadTreeForest::addTree(QuadTree tree) { trees.push_back(tree); } //���������� ������ � ������
void QuadTreeForest::addNodeToRefine(const NodeTag& t) { forest.toRefine[t.depth()].push_back(t); } //���������� ���� � ������ �� ���������
nodeEdgeId QuadTreeForest::addEdge(NodeEdge _edge) //�������� ����� � ������
{
    nodeEdgeId eid = getVacantEdgeId();
    edges[eid] = _edge;
    return eid;
}

nodeEdgeId QuadTreeForest::addEdgeUnique(NodeEdge _edge) //�������� ����� � ������ � ��������� ������������
{
    //����� ������ ����� ����� ���������
    for (nodeEdgeId eid = 0; eid < edges.size(); eid++)
    {
        if (edges[eid] == _edge)
            return eid;
    }
    //���� �� ������� - ���������� ������
    return addEdge(_edge);
}

void QuadTreeForest::updateEdge(nodeEdgeId eid, NodeTag n1, NodeTag n2) //���������� ������ �����
{
    try {
        if (edges.size() > eid && !edges[eid].isDeleted())
        {
            if (!n1.isNull()) //���� ������ �� null - ��������� ���
                edges[eid].setN1(n1);
            if (!n2.isNull())
                edges[eid].setN2(n2);
        }
        else
        {
            throw std::invalid_argument("edge is deleted or doesn't exist");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "forest.updateEdge() error: " << e.what() << " " << eid << endl;
    }
    return;
}
void QuadTreeForest::removeEdge(nodeEdgeId eid) //�������� ����� �� ������
{
    try {
        if (edges.size() > eid && !edges[eid].isDeleted())
        {
            edges[eid].markDeleted();
            vacant_edge_ids.push_back(eid);
        }
        else
        {
            throw std::invalid_argument("edge is deleted or doesn't exist");
        }
    }
    catch (const std::invalid_argument& e)
    {
        cout << "forest.removeEdge() error: " << e.what() << " " << eid << endl;
    }
    return;
}

//inspectors -----------------------
size_t QuadTreeForest::activeNodesNumber() const //������� �������� (�����������) ���
{
    size_t ret = 0;
    for (auto& rtree : forest.trees)
    {
        if (!rtree.isGhost())
            ret += rtree.active_nodes_num[0]; //� [0] �������� ��������� ����� �� ���� �������
    }
    return ret;
}
size_t QuadTreeForest::leavesNumber() const //������� �������
{
    size_t ret = 0;
    for (auto& rtree : forest.trees)
    {
        if (!rtree.isGhost())
            ret += rtree.leaf_nodes_num[0]; //� [0] �������� ��������� ����� �� ���� �������
    }
    return ret;
}

//other -----------------------
const int INDEX_MACH = TECPLOT_FIELDS_NUMBER - 3;
const int INDEX_LEVEL = TECPLOT_FIELDS_NUMBER - 2;
const int INDEX_MAGGRADRHO = TECPLOT_FIELDS_NUMBER - 1;
dataExtrema QuadTreeForest::getExtrema() //���� ����������� ���� ������� ��� ������ � Tecplot
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
    mins[0] = config.global_box.bottom_left().x + tree_h; //�������� �� ���� ���� ghost-��������
    maxs[0] = config.global_box.bottom_right().x - tree_h;
    mins[1] = config.global_box.bottom_left().y + tree_h;
    maxs[1] = config.global_box.top_left().y - tree_h;
    mins[INDEX_LEVEL] = 1;
    maxs[INDEX_LEVEL] = config.max_depth;
    for (auto& tree : trees)
    {
        if (tree.isGhost()) //������� ghost-��������
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
    //������� �� NaN
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

nodeEdgeId QuadTreeForest::getVacantEdgeId() //��������� ������ ��������� ������ ��� ������� ����� � ������� edges
{
    nodeEdgeId ret;
    //���� ��� ��������� ���� - ������ � ����� �������
    if (vacant_edge_ids.empty())
    {
        ret = edges.size(); //�� 1 ������ ������ ���������� ��������
        edges.emplace_back(); //��������� ������ (������ ���������� ����� .size())
    }
    else
    {
        ret = vacant_edge_ids.back(); //����� ��������� (�� ������� ������������) ��������� ������ � ������� edges
        vacant_edge_ids.pop_back(); //������ �� ������ ��������� �����
    }
    return ret;
}

void QuadTreeForest::meshApplyRefineList() //��������� ����� �� ������ toRefine � ������������
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

void QuadTreeForest::meshRefineInitial() //��������� ��������� �����
{
    if (config.meshRefineAll) //���� ����� ���������� ��
    {
        for (auto depth = 1; depth < config.max_depth; depth++)
        {
            cout << "Refining all, level " << depth << ": " << forest.leavesNumber();
            for (auto& rtree : trees) //������ �� ���� �������� �� ������
            {
                if (rtree.isGhost()) //������� ghost-��������
                    continue;
                for (auto& rnode : rtree.nodes[depth])
                    rnode.refine();
            }
            cout << " ---> " << forest.leavesNumber() << " leaves" << endl;
        }
        return;
    }
    //��������� �� ��������� ��������� ������
    for (auto depth = 1; depth < config.max_depth; depth++)
    {
        cout << "Refining level " << depth << ": " << forest.leavesNumber();
        for (auto& tree : trees) //������ �� ���� �������� �� ������
        {
            if (tree.isGhost()) //������� ghost-��������
                continue;
            if (depth <= tree.depth()) //���� ������ ����� ������
            {
                for (auto& rnode : tree.nodes[depth]) //������ �� ���� ����� �� ������
                {
                    bool to_refine = false;

                    //����� ��
                    if (rnode.box().intersectLineStraight(config.shock_position_x, Orientation::vertical) ||
                        rnode.box().bottom_left().isCloseToStraightLine(config.shock_position_x, Orientation::vertical) ||
                        rnode.box().bottom_right().isCloseToStraightLine(config.shock_position_x, Orientation::vertical))
                    {
                        to_refine = true;
                    }

                    //����� ������ ����
                    if (config.problem == "layer")
                    {
                        //������ �������
                        if (rnode.box().top_left().y >= config.layer_bottom - config.meshInitialRefinePadding &&
                            rnode.box().bottom_left().y <= config.layer_top + config.meshInitialRefinePadding &&
                            (rnode.box().intersectLineStraight(config.layer_right, Orientation::vertical) ||
                                rnode.box().bottom_left().isCloseToStraightLine(config.layer_right, Orientation::vertical) ||
                                rnode.box().top_right().isCloseToStraightLine(config.layer_right, Orientation::vertical)))
                        {
                            to_refine = true;
                        }

                        //������ � ������� �������
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

                    //����� ������� ������
                    if (config.problem == "bubble")
                    {
                        if (rnode.box().intersectLineEllipse(config.bubble_axle_x, config.bubble_axle_y))
                            to_refine = true;
                    }

                    //������� � ����������
                    if (to_refine)
                        rnode.markToRefine();
                }
            }
        }
        //��������� � ������������
        meshApplyRefineList();
        cout << " ---> " << forest.leavesNumber() << " leaves" << endl;
    }
    computeQuadraturePoints(); //������ ����� ���������� � ������
    updateEigenObjects(); //���������� ������ Eigen � �������
    return;
}

void QuadTreeForest::meshCoarsenInitial() //��������� ���������� ����� (��� �����)
{
    //���������� �� ������ � �������
    for (auto depth = config.max_depth - 1; depth > 0; depth--)
    {
        cout << "Coarsening level " << depth << ": " << forest.leavesNumber();
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //������� ghost-��������
                continue;
            if (depth <= rtree.depth())
            {
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (!rnode.isDeleted() && rnode.hasChildren())
                    {
                        //����� ������� ������
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
        //computeQuadraturePoints(); //������ ����� ���������� � ������
        //ExportForest();
    }
    computeQuadraturePoints(); //������ ����� ���������� � ������
    updateEigenObjects(); //���������� ������ Eigen � �������
}

void QuadTreeForest::meshUpdate() //���������� �����
{
    //��������� �� ������� � ������
    for (auto depth = 1; depth < config.max_depth; depth++)
    {
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //������� ghost-��������
                continue;
            if (depth <= rtree.depth()) //���� ������ ����� ������
            {
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (!rnode.isDeleted() && rnode.magGradRho() > globals.refineLvls[depth])
                        rnode.markToRefine();
                }
            }
        }
        //��������� � ������������
        meshApplyRefineList();
    }

    //���������� �� ������ � �������
    for (auto depth = config.max_depth - 1; depth > 0; depth--)
    {
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //������� ghost-��������
                continue;
            if (depth <= rtree.depth()) //���� ������ ����� ������
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
    computeQuadraturePoints(); //���������� ����� ���������� � ������
    updateEigenObjects(); //���������� ������ Eigen � �������
    return;
}

void QuadTreeForest::computeQuadraturePoints() //������ ����� ���������� �� ���� ������
{
    for (auto& redge : edges)
    {
        if (!redge.isDeleted())
        {
            redge.computeQuadraturePoints();
        }
    }
}

void QuadTreeForest::updateEigenObjects() //�������� ��� ���������� Eigen �������� ��� ���� �����
{
    for (auto& rtree : trees)
    {
        if (rtree.isGhost()) //������� ghost-��������
            continue;
        for (auto& rlevel : rtree.nodes)
        {
            for (auto& rnode : rlevel)
            {
                if (!rnode.isDeleted() && !rnode.hasChildren())
                {
                    rnode.updateEigenObjects();
                }
            }
        }
    }
    return;
}

void QuadTreeForest::initialCondition() //��������� �������
{
    double double_gamma = config.gamma;
    double double_shock_Mach_number = config.Mach;
    double double_Atwood_number = config.Atwood;
    double omega = (1.0 + double_Atwood_number) / (1.0 - double_Atwood_number);
    double p0 = 1.0;
    double r0 = 1.0;
    double u0 = double_shock_Mach_number * sqrt(double_gamma); //�������� ��������� ������
    //��������� ���� �� �������� �������
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
            if (depth <= rtree.depth()) //���� ������ ����� ������
            {
                //������ �� ���� ����� �� ������
                for (auto& rnode : rtree.nodes[depth])
                {
                    if (rnode.hasChildren()) //��-������ ������������
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
                    //��������� �� r ��� ��������������� ���������
                    if (config.coord_type == CoordType::axisymmetric)
                    {
                        rd.setY(y);
                        for (auto eq : Equations)
                        {
                            Q.set(eq, Q(eq) * y);
                        }
                    }
                    //cout << "(" << x << ", " << y << "): " << Q << endl;
                    rd.setQ(Q); //������ � CellData (rk = 0)
                }
            }
        }
    }
    return;
}

double QuadTreeForest::CFLTimestepSize()
{
    double ret = std::numeric_limits<double>::infinity();
    for (auto& rtree : trees)
    {
        if (rtree.isGhost()) //������� ghost-��������
            continue;
        //������ �� ���� ����� �� ������ ������
        for (auto& rnode : rtree.nodes[rtree.depth()])
        {
            CellData& d = rnode.dataRef();
            double a = std::sqrt(config.gamma * d.p() / d.rho());
            double dx = rnode.box().size();
            double vmax = std::max(fabs(d.u()), fabs(d.v()));
            double dt = dx / (vmax + a);
            ret = std::min(ret, dt);
        }
    }
    ret *= config.Courant;
    return ret;
}

void QuadTreeForest::putQn(rkStep rk) //�������� Qn[rk_order] -> Qn[0]
{
    for (auto& rtree : trees)
    {
        for (auto& rlevel : rtree.nodes)
        {
            for (auto& rnode : rlevel)
            {
                if (!rnode.isDeleted() && !rnode.hasChildren())
                {
                    rnode.dataRef().putQn(rk);
                }
            }
        }
    }
}

void QuadTreeForest::computePolynomialCoeffs(rkStep rk) //���������� ������������� CWENO �������� �� ���� �������
{
    for (auto& rtree : trees)
    {
        for (auto& rlevel : rtree.nodes)
        {
            for (auto& rnode : rlevel)
            {
                if (!rnode.isDeleted() && !rnode.hasChildren())
                {
                    rnode.calcPolynomialCWENO(rk);
                }
            }
        }
    }
}

void QuadTreeForest::computeFluxesCWENO(rkStep rk) //������ ������� �� ���� ������
{
    //int edges_total = forest.edges.size();
    //#pragma omp parallel for TODO: ����������, ������ �������� ���� t -> inf
    for (auto& redge : edges)
        //for (int id = 0; id < edges_total; id++)
    {
        //auto& edge = forest.edges[id];
        if (redge.isDeleted() || (treeRef(redge.n1().tree()).isGhost() && treeRef(redge.n2().tree()).isGhost())) //��������� ����� � ����� ����� ghost'��� ������������
            continue;
        if (config.flux == FluxType::LF)
            redge.computeFluxLF(rk);
        /*else if (config.flux == FluxType::Riemann_exact)
            redge.computeFluxRiemannExact(rk);
        else if (config.flux == FluxType::HLLC)
            redge.computeFluxHLLC(rk);*/
    }
    return;
}

void QuadTreeForest::advanceTime() //������ ��� �� �������: ���������� Qn[rk_order]
{
    for (rkStep rk = 0; rk < config.rk_order; rk++)
    {
        //CheckNaNQns(rk);
        computePolynomialCoeffs(rk);
        computeFluxesCWENO(rk);
        //CheckNaNFluxes(rk);
        for (auto& rtree : trees)
        {
            if (rtree.isGhost()) //������� ghost-��������
                continue;
            for (auto& rlevel : rtree.nodes)
            {
                for (auto& rnode : rlevel)
                {
                    if (!rnode.isDeleted() && !rnode.hasChildren())
                    {
                        auto& rd = rnode.dataRef();
                        double h = rnode.box().size();
                        double area = h * h; //������� ������ (�����, �������� �� dz)
                        if (config.coord_type == CoordType::axisymmetric) //� ��������������� ����������� ����������� �� r (V = r dr dz dtheta)
                            area *= rnode.box().center().y;
                        ConservativeVector Q{};
                        for (auto eq : Equations)
                        {
                            for (rkStep r = 0; r <= rk; r++) //��������� �� Qn
                                Q.set(eq, Q(eq) + RKcoeff[config.rk_order - 1][rk][r] * rd.Qref(r)(eq));

                            //FV-�����: ����� ���������� ������� �� ���� ������
                            /*for (ushorty etype = 0; etype < DIRECTIONS_NUM * 2; etype++)
                            {
                                double sgn; //���� ������: ����� � ����� "+", ������ � ������ "-"
                                if (etype < DIRECTIONS_NUM)
                                    sgn = -1.0;
                                else
                                    sgn = 1.0;
                                if (cnode.hasEdge(etype))
                                {
                                    auto& cedge = cnode.getEdge(etype);
                                    double dflux = RKcoeff[config.rk_order - 1][rk][RK_ORDER_MAX] * globals.dt / area * sgn * cedge.FQ[eq];
                                    if (config.coord_type == COORD_TYPE_AXISYMMETRIC) //��������������� ����������
                                    {
                                        if (cedge.orientation == ORIENTATION_VERTICAL) //����� ��� x 
                                        {
                                            d.Qn[rk + 1][eq] += dflux * d.y;
                                        }
                                        else //����� ��� y
                                        {
                                            if (etype == 0 || etype == 1) //������� �����: ������� ����� ������������ ������ ������
                                                d.Qn[rk + 1][eq] += dflux * (d.y + h / 2.0);
                                            else //������ �����: ������� ����� ������
                                                d.Qn[rk + 1][eq] += dflux * (d.y - h / 2.0);
                                        }
                                    }
                                    else
                                        d.Qn[rk + 1][eq] += dflux;
                                }
                            }
                            if (config.coord_type == COORD_TYPE_AXISYMMETRIC && eq == EQ_MOMENTUM_Y) //������������ ����
                            {
                                d.Qn[rk + 1][eq] += RKcoeff[config.rk_order - 1][rk][RK_ORDER_MAX] * d.p(rk);
                            }*/
                        }
                        rd.setQ(Q, rk + 1); //������ rd.Qn[rk + 1]
                    }
                }
            }
        }
        if (rk + 1 < config.rk_order) //��� ����, ����� ���������� (��� ���� ���������� � main)
            boundaryConditions(rk + 1);
    }
    return;
}

void QuadTreeForest::boundaryConditions(rkStep rk) //��������� �������
{
    for (auto& rtree : trees)
    {
        if (!rtree.isGhost()) //����� ������� ������������
            continue;
        if (rtree.isGhostCorner()) //������� �������
        {
            TreeNode& rroot = rtree.rootRef(); //�������������� ������� ������ �������� ������
            if (rroot.hasNeighbour12(Neighbour12::top_right)) //��� ������ - ������ �����
            {
                ConservativeVector Q = TreeNode::nodeRef(rroot.neighbour12(Neighbour12::top_right)).dataRef().Qref(rk); //�� ��������
                if (config.boundary_conditions[static_cast<int>(Directions::down)] == BCType::wall)
                    Q.flipVelocity(Orientation::vertical);
                if (config.boundary_conditions[static_cast<int>(Directions::left)] == BCType::wall)
                    Q.flipVelocity(Orientation::horizontal);
                rroot.dataRef().setQ(Q, rk);
                continue;
            }
            if (rroot.hasNeighbour12(Neighbour12::bottom_right))
            {
                ConservativeVector Q = TreeNode::nodeRef(rroot.neighbour12(Neighbour12::bottom_right)).dataRef().Qref(rk);
                if (config.boundary_conditions[static_cast<int>(Directions::left)] == BCType::wall)
                    Q.flipVelocity(Orientation::horizontal);
                if (config.boundary_conditions[static_cast<int>(Directions::up)] == BCType::wall)
                    Q.flipVelocity(Orientation::vertical);
                rroot.dataRef().setQ(Q, rk);
                continue;
            }
            if (rroot.hasNeighbour12(Neighbour12::bottom_left))
            {
                ConservativeVector Q = TreeNode::nodeRef(rroot.neighbour12(Neighbour12::bottom_left)).dataRef().Qref(rk);
                if (config.boundary_conditions[static_cast<int>(Directions::up)] == BCType::wall)
                    Q.flipVelocity(Orientation::vertical);
                if (config.boundary_conditions[static_cast<int>(Directions::right)] == BCType::wall)
                    Q.flipVelocity(Orientation::horizontal);
                rroot.dataRef().setQ(Q, rk);
                continue;
            }
            if (rroot.hasNeighbour12(Neighbour12::top_left))
            {
                ConservativeVector Q = TreeNode::nodeRef(rroot.neighbour12(Neighbour12::top_left)).dataRef().Qref(rk);
                if (config.boundary_conditions[static_cast<int>(Directions::right)] == BCType::wall)
                    Q.flipVelocity(Orientation::horizontal);
                if (config.boundary_conditions[static_cast<int>(Directions::down)] == BCType::wall)
                    Q.flipVelocity(Orientation::vertical);
                rroot.dataRef().setQ(Q, rk);
                continue;
            }
        }
        //���������, ��������� �������
        for (auto& rlevel : rtree.nodes)
        {
            for (auto& rnode : rlevel)
            {
                if (!rnode.isDeleted() && !rnode.hasChildren())
                {
                    for (auto n : Neighbours)
                    {
                        if (rnode.hasNeighbour(n) && !forest.treeRef(rnode.neighbour(n).tree()).isGhost()) //���� ����� �����
                        {
                            ConservativeVector Q = TreeNode::nodeRef(rnode.neighbour(n)).dataRef().Qref(rk); //�� ��������
                            switch (n)
                            {
                            case Neighbour::top:
                                if (config.coord_type == CoordType::cartesian && config.boundary_conditions[static_cast<int>(Directions::down)] == BCType::wall) //����� ����� ������ = ������ ������� ��������� �����
                                    Q.flipVelocity(Orientation::vertical);
                                break;
                            case Neighbour::right:
                                if (config.boundary_conditions[static_cast<int>(Directions::left)] == BCType::wall) //����� �������
                                    Q.flipVelocity(Orientation::horizontal);
                                break;
                            case Neighbour::bottom:
                                if (config.boundary_conditions[static_cast<int>(Directions::up)] == BCType::wall) //������� �������
                                    Q.flipVelocity(Orientation::vertical);
                                break;
                            case Neighbour::left:
                                if (config.boundary_conditions[static_cast<int>(Directions::right)] == BCType::wall) //������ �������
                                    Q.flipVelocity(Orientation::horizontal);
                                break;
                            }
                            rnode.dataRef().setQ(Q, rk); //������ ������� ������
                        }
                    }
                }
            }
        }
    }
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

void QuadTreeForest::exportEdges(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;

    //������� �������
    file_output << "ZONE T = \"Forest edge neighbours\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;
    for (auto& edge : forest.edges)
    {
        if (!edge.isDeleted())
        {
            file_output << edge.dumpNeigboursVectors();
            auto& rn1 = TreeNode::nodeRef(edge.n1());
            auto& rn2 = TreeNode::nodeRef(edge.n2());
            if (rn1.hasChildren() || rn2.hasChildren())
                cout << "suspicious edge! " << edge << endl;
        }
    }

    //����� ����������
    file_output << endl << "ZONE T = \"Edge quadrature points\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;
    for (auto& edge : forest.edges)
    {
        if (!edge.isDeleted())
        {
            file_output << edge.dumpQuadraturePoints();
        }
    }
    return;

}

void QuadTreeForest::exportNodeEdges(std::string filename) //����� ����� ����� � ���� (Tecplot ASCII vectors)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;
    file_output << "ZONE T = \"Node edges\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& rtree : trees)
    {
        for (auto& rnodes_level : rtree.nodes)
        {
            for (auto& rnode : rnodes_level)
            {
                if (!rnode.isDeleted() && rnode.isLeaf())
                {
                    for (auto e : Edges)
                    {
                        if (rnode.hasEdge(e))
                            file_output << rnode.dumpEdgeVector(e);
                    }
                }
            }
        }
    }
    file_output.close();
    return;
}




