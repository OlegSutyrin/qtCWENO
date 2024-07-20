#include "main.h"
#include "globalFuncs.h"
#include "CellBox.h"
#include "NodeTag.h"
#include "TreeNode.h"
#include "NodeEdge.h"

#include <sstream>

double NodeEdge::FQ(Equation eq) const { return FQ_[static_cast<int>(eq)]; } //���������� ������
double NodeEdge::length() const { return distance(qps[0], qps[1]) / ISQ3; } //���������� ����� ����� �� ������ ���������� TODO: ������� ������� � ����� � ����� �����

void NodeEdge::computeQuadraturePoints() //������ ����� ���������� �� �����
{
    auto& rnode1 = TreeNode::nodeRef(n1_);
    auto& rnode2 = TreeNode::nodeRef(n2_);
    if (rnode1.box().size() <= rnode2.box().size()) //���� ������ ����� ������� �� ���� (��� ��� �����)
    {
        if (orientation_ == Orientation::vertical)
            qps = quadraturePoints(rnode1.box().bottom_right(), rnode1.box().top_right()); //��� ����� - ������ ������� ������ ������
        else
            qps = quadraturePoints(rnode1.box().bottom_left(), rnode1.box().bottom_right()); //������ ������� ��������
    }
    else
    {
        if (orientation_ == Orientation::vertical)
            qps = quadraturePoints(rnode2.box().bottom_left(), rnode2.box().top_left()); //����� ������� �������
        else
            qps = quadraturePoints(rnode2.box().top_left(), rnode2.box().top_right()); //������� ������� �������
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

