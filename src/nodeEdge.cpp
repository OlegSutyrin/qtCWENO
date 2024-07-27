#include "main.h"
#include "globalFuncs.h"
#include "CellBox.h"
#include "NodeTag.h"
#include "TreeNode.h"
#include "NodeEdge.h"

#include <sstream>

const NodeTag NodeEdge::n1() const { return n1_; }
const NodeTag NodeEdge::n2() const { return n2_; }


//accessors
double NodeEdge::FQ(Equation eq) const { return FQ_[static_cast<int>(eq)]; } //компонента потока
Point NodeEdge::middle() { return ::middle(qps[0], qps[1]); } //(:: - глобальное namespace)
double NodeEdge::length() const { return distance(qps[0], qps[1]) / ISQ3; } //вычисление длины ребра по точкам квадратуры TODO: хранить вершины и длину в самом ребре

//mutators
void NodeEdge::markDeleted() { is_deleted = true; }
void NodeEdge::setN1(NodeTag ntag) { n1_ = ntag; }
void NodeEdge::setN2(NodeTag ntag) { n2_ = ntag; }

//inspectors
bool NodeEdge::isDeleted() { return is_deleted; }

//other
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

bool NodeEdge::operator==(const NodeEdge& rhs) //equality overload
{
    if (n1_ == rhs.n1_ && n2_ == rhs.n2_)
        return true;
    return false;
}

//output
std::string NodeEdge::dumpNeigboursVectors() //дамп соседей в виде векторов
{
    std::stringstream buffer;
    auto& rn1 = TreeNode::nodeRef(n1_);
    auto& rn2 = TreeNode::nodeRef(n2_);
    //стартовая точка вектора
    double x, y;
    if (rn1.tag().depth() > rn2.tag().depth()) //rn1 - меньшая ячейка
    {
        if (orientation_ == Orientation::vertical)
        {
            x = rn1.box().center().x + (rn2.box().center().x - rn1.box().center().x) / 3.0;
            y = rn1.box().center().y;
        }
        else
        {
            x = rn1.box().center().x;
            y = rn1.box().center().y + (rn2.box().center().y - rn1.box().center().y) / 3.0;
        }
    }
    else if (rn1.tag().depth() < rn2.tag().depth()) //rn2 - меньшая ячейка
    {
        if (orientation_ == Orientation::vertical)
        {
            x = rn2.box().center().x + (rn1.box().center().x - rn2.box().center().x) / 3.0;
            y = rn2.box().center().y;
        }
        else
        {
            x = rn2.box().center().x;
            y = rn2.box().center().y + (rn1.box().center().y - rn2.box().center().y) / 3.0;
        }
    }
    else //одинаковые ячейки
    {
        x = (rn1.box().center().x + rn2.box().center().x) / 2.0;
        y = (rn1.box().center().y + rn2.box().center().y) / 2.0;
    }

    //векторы
    double dx = 0.8 * (rn1.box().center().x - x);
    double dy = 0.8 * (rn1.box().center().y - y);
    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;
    dx = 0.8 * (rn2.box().center().x - x);
    dy = 0.8 * (rn2.box().center().y - y);
    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;

    return buffer.str();
}

std::string NodeEdge::dumpQuadraturePoints() //дамп точек квадратуры
{
    std::stringstream buffer;
    //координаты точек (и нули для dx, dy)
    buffer << std::to_string(qps[0].x) << " " << std::to_string(qps[0].y) << " 0 0 " << endl;
    buffer << std::to_string(qps[1].x) << " " << std::to_string(qps[1].y) << " 0 0 " << endl;

    return buffer.str();
}


std::ostream& operator<<(std::ostream& os, const NodeEdge& edge) //output overload
{
    os << edge.n1_;
    if (edge.orientation_ == Orientation::vertical)
        os << " | ";
    else
        os << " - ";
    os << edge.n2_;
    os << ", flux: (" << edge.FQ(Equation::density) << ", " << edge.FQ(Equation::momentum_x) << ", " << edge.FQ(Equation::momentum_y) << ", " << edge.FQ(Equation::energy) << ")";
    return os;
}

