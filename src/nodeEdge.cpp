#include "main.h"
#include "globalFuncs.h"
#include "CellBox.h"
#include "NodeTag.h"
#include "TreeNode.h"
#include "NodeEdge.h"

#include <sstream>

const NodeTag NodeEdge::n1() const { return n1_; }
const NodeTag NodeEdge::n2() const { return n2_; }
const Orientation NodeEdge::ori() const { return orientation_; }

//accessors
double NodeEdge::FQ(Equation eq) const { return FQ_[static_cast<int>(eq)]; } //���������� ������
Point NodeEdge::middle() { return ::middle(qps[0], qps[1]); } //(:: - ���������� namespace)
double NodeEdge::length() const { return distance(qps[0], qps[1]) / ISQ3; } //���������� ����� ����� �� ������ ���������� TODO: ������� ������� � ����� � ����� �����

//mutators
void NodeEdge::markDeleted() { is_deleted = true; }
void NodeEdge::setN1(NodeTag ntag) { n1_ = ntag; }
void NodeEdge::setN2(NodeTag ntag) { n2_ = ntag; }

//inspectors
bool NodeEdge::isDeleted() { return is_deleted; }

//other
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

void NodeEdge::computeFluxLF(rkStep rk) //������ ������ (Lax-Friedrich flux)
{
    auto& rn1 = TreeNode::nodeRef(n1_);
    auto& rn2 = TreeNode::nodeRef(n2_);
    double h = length(); //����� �����
    
    //TODO: ���������, ��� ����� �������� � r � �������������� �����������
    double y1 = 1.0, y2 = 1.0;
    if (config.coord_type == CoordType::axisymmetric)
    {
        y1 = rn1.box().center().y;
        y2 = rn2.box().center().y;
    }
    
    //������������������ ������ ������� � ������ ����������, [�����][�����]
    ConservativeVector Qs[2][2]{};
    CellData gdata{}; //ghost-������ ��� ��������� �������
    double gp�fs[EQ_NUM][POLY_COEFF_NUM]; //������������ ����������� ghost-������
    Point gcenter{};
    if (rn1.isGhost()) //������ �������� ������ - ghost
    {
        Qs[0][1] = rn2.evalPolynomialAt(qps[0], rk); //rn2 �����
        Qs[1][1] = rn2.evalPolynomialAt(qps[1], rk);
        //����� ������ ������ rn2
        gdata = rn2.data();
        ConservativeVector& gQ = gdata.Qref(rk);
        for (auto eq : Equations)
            for (int i = 0; i < POLY_COEFF_NUM; i++)
                gp�fs[static_cast<int>(eq)][i] = rn2.polyCoeff(eq, i);
        double h = rn2.box().size();
        if(orientation_ == Orientation::vertical) //ghost-������ �����
        { 
            gcenter.x = rn2.box().center().x - h;
            gcenter.y = rn2.box().center().y;
            if (config.boundary_conditions[static_cast<int>(Directions::left)] == BCType::wall) //������������
            {
                gQ.flipVelocity(Orientation::horizontal); //��������� ����� �������� ����� ��� x
                gp�fs[static_cast<int>(Equation::momentum_x)][0] *= -1; //px = -px, ���������� ��� ��� ����
                gp�fs[static_cast<int>(Equation::momentum_x)][2] *= -1; //pxx = -pxx
                gp�fs[static_cast<int>(Equation::momentum_x)][3] *= -1; //pxy = -pxy, ���������� ��� ��� ����
                for (auto eq : Equations) //��������� ���� ������� ������������ ������������ ���
                {
                    gp�fs[static_cast<int>(eq)][0] *= -1; //px = -px
                    gp�fs[static_cast<int>(eq)][3] *= -1; //pxy = -pxy
                }
            }
            else //d/dx=0
            {
                for (auto eq : Equations)
                {
                    gp�fs[static_cast<int>(eq)][0] = 0; //px = 0
                    gp�fs[static_cast<int>(eq)][2] = 0; //pxx = 0
                    gp�fs[static_cast<int>(eq)][3] = 0; //pxy = 0
                }
            }
        }
        else //ghost-������ ������
        {
            gcenter.x = rn2.box().center().x;
            gcenter.y = rn2.box().center().y + h;
            if (config.boundary_conditions[static_cast<int>(Directions::up)] == BCType::wall) //������������
            {
                gQ.flipVelocity(Orientation::vertical); //��������� ����� �������� ����� ��� y
                gp�fs[static_cast<int>(Equation::momentum_y)][1] *= -1; //py = -py, ���������� ��� ��� ����
                gp�fs[static_cast<int>(Equation::momentum_y)][3] *= -1; //pxy = -pxy, ���������� ��� ��� ����
                gp�fs[static_cast<int>(Equation::momentum_y)][4] *= -1; //pyy = -pyy
                for (auto eq : Equations) //��������� ���� ������� ������������ �������������� ���
                {
                    gp�fs[static_cast<int>(eq)][1] *= -1; //py = -py
                    gp�fs[static_cast<int>(eq)][3] *= -1; //pxy = -pxy
                }
            }
            else //d/dy=0
            {
                for (auto eq : Equations)
                {
                    gp�fs[static_cast<int>(eq)][1] = 0; //py = 0
                    gp�fs[static_cast<int>(eq)][3] = 0; //pxy = 0
                    gp�fs[static_cast<int>(eq)][4] = 0; //pyy = 0
                }
            }
        }
        for (int p = 0; p < 2; p++) //���� �� ������ ����������
        {
            double dx = qps[p].x - gcenter.x;
            double dy = qps[p].y - gcenter.y;
            for (auto eq : Equations)
            {
                Qs[p][0].set(eq, gQ(eq) + gp�fs[static_cast<int>(eq)][0] * dx + gp�fs[static_cast<int>(eq)][1] * dy
                    + 0.5 * gp�fs[static_cast<int>(eq)][2] * (dx * dx - 1.0 / 12.0 * h * h) + gp�fs[static_cast<int>(eq)][3] * dx * dy + 0.5 * gp�fs[static_cast<int>(eq)][4] * (dy * dy - 1.0 / 12.0 * h * h));
            }
        }
    }
    else if (rn2.isGhost()) //������ �������� ������ - ghost
    {
        Qs[0][0] = rn1.evalPolynomialAt(qps[0], rk); //rn1 �����
        Qs[1][0] = rn1.evalPolynomialAt(qps[1], rk);
        //����� ������ ������ rn1
        gdata = rn1.data();
        ConservativeVector& gQ = gdata.Qref(rk);
        for (auto eq : Equations)
            for (int i = 0; i < POLY_COEFF_NUM; i++)
                gp�fs[static_cast<int>(eq)][i] = rn1.polyCoeff(eq, i);
        double h = rn1.box().size();
        if (orientation_ == Orientation::vertical) //ghost-������ ������
        {
            gcenter.x = rn1.box().center().x + h;
            gcenter.y = rn1.box().center().y;
            if (config.boundary_conditions[static_cast<int>(Directions::right)] == BCType::wall) //������������
            {
                gQ.flipVelocity(Orientation::horizontal); //��������� ����� �������� ����� ��� x
                gp�fs[static_cast<int>(Equation::momentum_x)][0] *= -1; //px = -px, ���������� ��� ��� ����
                gp�fs[static_cast<int>(Equation::momentum_x)][2] *= -1; //pxx = -pxx
                gp�fs[static_cast<int>(Equation::momentum_x)][3] *= -1; //pxy = -pxy, ���������� ��� ��� ����
                for (auto eq : Equations) //��������� ���� ������� ������������ ������������ ���
                {
                    gp�fs[static_cast<int>(eq)][0] *= -1; //px = -px
                    gp�fs[static_cast<int>(eq)][3] *= -1; //pxy = -pxy
                }
            }
            else //d/dx=0
            {
                for (auto eq : Equations)
                {
                    gp�fs[static_cast<int>(eq)][0] = 0; //px = 0
                    gp�fs[static_cast<int>(eq)][2] = 0; //pxx = 0
                    gp�fs[static_cast<int>(eq)][3] = 0; //pxy = 0
                }
            }
        }
        else //ghost-������ �����
        {
            gcenter.x = rn1.box().center().x;
            gcenter.y = rn1.box().center().y - h;
            if (config.boundary_conditions[static_cast<int>(Directions::down)] == BCType::wall) //������������
            {
                gQ.flipVelocity(Orientation::vertical); //��������� ����� �������� ����� ��� y
                gp�fs[static_cast<int>(Equation::momentum_y)][1] *= -1; //py = -py, ���������� ��� ��� ����
                gp�fs[static_cast<int>(Equation::momentum_y)][3] *= -1; //pxy = -pxy, ���������� ��� ��� ����
                gp�fs[static_cast<int>(Equation::momentum_y)][4] *= -1; //pyy = -pyy
                for (auto eq : Equations) //��������� ���� ������� ������������ �������������� ���
                {
                    gp�fs[static_cast<int>(eq)][1] *= -1; //py = -py
                    gp�fs[static_cast<int>(eq)][3] *= -1; //pxy = -pxy
                }
            }
            else //d/dy=0
            {
                for (auto eq : Equations)
                {
                    gp�fs[static_cast<int>(eq)][1] = 0; //py = 0
                    gp�fs[static_cast<int>(eq)][3] = 0; //pxy = 0
                    gp�fs[static_cast<int>(eq)][4] = 0; //pyy = 0
                }
            }
        }
        for (int p = 0; p < 2; p++) //���� �� ������ ����������
        {
            double dx = qps[p].x - gcenter.x;
            double dy = qps[p].y - gcenter.y;
            for (auto eq : Equations)
            {
                Qs[p][1].set(eq, gQ(eq) + gp�fs[static_cast<int>(eq)][0] * dx + gp�fs[static_cast<int>(eq)][1] * dy
                    + 0.5 * gp�fs[static_cast<int>(eq)][2] * (dx * dx - 1.0 / 12.0 * h * h) + gp�fs[static_cast<int>(eq)][3] * dx * dy + 0.5 * gp�fs[static_cast<int>(eq)][4] * (dy * dy - 1.0 / 12.0 * h * h));
            }
        }
    }
    else //��� �������� ������ �����
    {
        Qs[0][0] = rn1.evalPolynomialAt(qps[0], rk);
        Qs[0][1] = rn2.evalPolynomialAt(qps[0], rk);
        Qs[1][0] = rn1.evalPolynomialAt(qps[1], rk);
        Qs[1][1] = rn2.evalPolynomialAt(qps[1], rk);
    }
    //ConservativeVector Qs[2][2] = { {rn1.evalPolynomialAt(qps[0], rk), rn2.evalPolynomialAt(qps[0], rk)},
        //{rn1.evalPolynomialAt(qps[1], rk), rn2.evalPolynomialAt(qps[1], rk)} };
    double tmpFQ[2] = { 0, 0 };
    if (orientation_ == Orientation::vertical) //����� ����� ��� X
    {
        double lmax = std::max(std::max(Qs[0][0].lambdaGLFX(y1), Qs[0][1].lambdaGLFX(y2)),
            std::max(Qs[1][0].lambdaGLFX(y1), Qs[1][1].lambdaGLFX(y2)));
        for (auto eq : Equations)
        {
            for (int p = 0; p < 2; p++)
            {
                tmpFQ[p] = 0.5 * (Qs[p][0].F(eq, y1) + Qs[p][1].F(eq, y2) + lmax * (Qs[p][0](eq) - Qs[p][1](eq))); //GLF2 flux
            }
            FQ_[static_cast<int>(eq)] = 0.5 * h * (tmpFQ[0] + tmpFQ[1]); //FQ = flux integral over edge (Gaussian quadrature)
        }
    }
    else  //����� ��� Y
    {
        double lmax = std::max(std::max(Qs[0][0].lambdaGLFY(y1), Qs[0][1].lambdaGLFY(y2)),
            std::max(Qs[1][0].lambdaGLFY(y1), Qs[1][1].lambdaGLFY(y2)));
        for (auto eq : Equations)
        {
            for (int p = 0; p < 2; p++)
            {
                tmpFQ[p] = 0.5 * (Qs[p][0].G(eq, y1) + Qs[p][1].G(eq, y2) + lmax * (Qs[p][1](eq) - Qs[p][0](eq))); //��������� ������� ������ ����
            }
            FQ_[static_cast<int>(eq)] = 0.5 * h * (tmpFQ[0] + tmpFQ[1]);
        }
    }
}

bool NodeEdge::operator==(const NodeEdge& rhs) //equality overload
{
    if (n1_ == rhs.n1_ && n2_ == rhs.n2_)
        return true;
    return false;
}

//output
std::string NodeEdge::dumpNeigboursVectors() //���� ������� � ���� ��������
{
    std::stringstream buffer;
    auto& rn1 = TreeNode::nodeRef(n1_);
    auto& rn2 = TreeNode::nodeRef(n2_);
    //��������� ����� �������
    double x, y;
    if (rn1.tag().depth() > rn2.tag().depth()) //rn1 - ������� ������
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
    else if (rn1.tag().depth() < rn2.tag().depth()) //rn2 - ������� ������
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
    else //���������� ������
    {
        x = (rn1.box().center().x + rn2.box().center().x) / 2.0;
        y = (rn1.box().center().y + rn2.box().center().y) / 2.0;
    }

    //�������
    double dx = 0.8 * (rn1.box().center().x - x);
    double dy = 0.8 * (rn1.box().center().y - y);
    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;
    dx = 0.8 * (rn2.box().center().x - x);
    dy = 0.8 * (rn2.box().center().y - y);
    buffer << std::to_string(x) << " " << std::to_string(y) << " " << std::to_string(dx) << " " << std::to_string(dy) << endl;

    return buffer.str();
}

std::string NodeEdge::dumpQuadraturePoints() //���� ����� ����������
{
    std::stringstream buffer;
    //���������� ����� (� ���� ��� dx, dy)
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

