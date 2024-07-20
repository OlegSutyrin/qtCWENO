#include "main.h"
#include "globalFuncs.h"
#include "CellBox.h"
#include "NodeTag.h"
#include "TreeNode.h"
#include "NodeEdge.h"

#include <sstream>

double NodeEdge::FQ(Equation eq) const { return FQ_[static_cast<int>(eq)]; } //компонента потока
double NodeEdge::length() const { return distance(qps[0], qps[1]) / ISQ3; } //вычисление длины ребра по точкам квадратуры TODO: хранить вершины и длину в самом ребре

void NodeEdge::computeQuadraturePoints() //расчет точек квадратуры на ребре
{
    auto& rnode1 = TreeNode::nodeRef(n1_);
    auto& rnode2 = TreeNode::nodeRef(n2_);
    if (rnode1.box().size() <= rnode2.box().size()) //если первый сосед меньший из двух (или они равны)
    {
        if (orientation_ == Orientation::vertical)
            qps = quadraturePoints(rnode1.box().bottom_right(), rnode1.box().top_right()); //это ребро - правая сторона левого соседа
        else
            qps = quadraturePoints(rnode1.box().bottom_left(), rnode1.box().bottom_right()); //нижняя сторона верхнего
    }
    else
    {
        if (orientation_ == Orientation::vertical)
            qps = quadraturePoints(rnode2.box().bottom_left(), rnode2.box().top_left()); //левая сторона правого
        else
            qps = quadraturePoints(rnode2.box().top_left(), rnode2.box().top_right()); //верхняя сторона нижнего
    }
}

std::ostream& operator<<(std::ostream& os, const NodeEdge& edge) //output overload
{
    os << edge.n1_;
    if (edge.orientation_ == Orientation::vertical)
        os << " | ";
    else
        os << " - ";
    os << edge.n1_;
    os << ", flux: (" << edge.FQ(Equation::density) << ", " << edge.FQ(Equation::momentum_x) << ", " << edge.FQ(Equation::momentum_y) << ", " << edge.FQ(Equation::energy) << ")";
    return os;
}

