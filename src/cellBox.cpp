#include <utility> //for std::pair

#include "main.h"
#include "CellBox.h"

//------------------------------------------------------------------- Point -------------------------------------------------------------------
bool Point::isCloseToStraightLine(double line_coord, Orientation ori) const //����� �� ����� ������ ������ ������� 2*(refine padding) ����� ������ �����
{
    if (ori == Orientation::horizontal)
    {
        if (fabs(y - line_coord) < config.meshInitialRefinePadding)
            return true;
    }
    else
    {
        if (fabs(x - line_coord) < config.meshInitialRefinePadding)
            return true;
    }
    return false;
}
std::ostream& operator<<(std::ostream& os, const Point& p) //output overload
{
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}

Point middle(Point const& p1, Point const& p2) //�������� ������� ����� �������
{
    return { (p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0 };
}
double distanceSquared(Point const& p1, Point const& p2) //������� ���������� ����� �������
{
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}
double distance(Point const& p1, Point const& p2) //���������� ����� �������
{
    return sqrt(distanceSquared(p1, p2));
}

std::pair<Point, Point> quadraturePoints(Point const& p1, Point const& p2) //��� ����� ��������� ���������� �� �������
{
    std::pair<Point, Point> ret;
    double dx = fabs(p2.x - p1.x);
    double dy = fabs(p2.y - p1.y);
    ret.first  = { middle(p1, p2).x - ISQ3 * 0.5 * dx, middle(p1, p2).y - ISQ3 * 0.5 * dy };
    ret.second = { middle(p1, p2).x + ISQ3 * 0.5 * dx, middle(p1, p2).y + ISQ3 * 0.5 * dy };
    return ret;
}

//------------------------------------------------------------------- CellBox -------------------------------------------------------------------
CellBox::CellBox(Point bl, Point tr) //����������� �� bottom_left � top_right ������
{
    vertices[static_cast<int>(Quadrant::top_left)] = { bl.x, tr.y };
    vertices[static_cast<int>(Quadrant::top_right)] = tr;
    vertices[static_cast<int>(Quadrant::bottom_right)] = { tr.x, bl.y };
    vertices[static_cast<int>(Quadrant::bottom_left)] = bl;
    updateCenter();
}

//point CellBox::getP(Quadrant q) const { return vertices[static_cast<int>(q)]; } //��������� ������ box'�
Point CellBox::top_left() const { return vertices[static_cast<int>(Quadrant::top_left)]; }
Point CellBox::top_right() const { return vertices[static_cast<int>(Quadrant::top_right)]; }
Point CellBox::bottom_right() const { return vertices[static_cast<int>(Quadrant::bottom_right)]; }
Point CellBox::bottom_left() const { return vertices[static_cast<int>(Quadrant::bottom_left)]; }
Point CellBox::center() const { return center_; } //��������� ������

//void CellBox::setP(Quadrant q, point pt) { vertices[static_cast<int>(q)] = pt; } //������� �������

double CellBox::size() const { return bottom_right().x - bottom_left().x; }

void CellBox::updateCenter() //���������� ������ �� ������� ������
{ 
    center_ = middle(vertices[static_cast<int>(Quadrant::bottom_left)], vertices[static_cast<int>(Quadrant::top_right)]); 
}

bool CellBox::isPointInside(Point p) const //�������� �� ����� � box (������� ����� � ������� �������)
{
    if (p.x >= vertices[static_cast<int>(Quadrant::bottom_left)].x && p.x < vertices[static_cast<int>(Quadrant::bottom_right)].x &&
        p.y >= vertices[static_cast<int>(Quadrant::bottom_left)].y && p.y < vertices[static_cast<int>(Quadrant::top_left)].y)
        return true;
    return false;
}

CellBox CellBox::quarterBox(Quadrant q) //��������� ����������� box'�
{
    switch (q)
    {
    case Quadrant::top_left:
        return CellBox(middle(top_left(), bottom_left()), middle(top_left(), top_right()));
        break;
    case Quadrant::top_right:
        return CellBox(center_, top_right());
        break;
    case Quadrant::bottom_right:
        return CellBox(middle(bottom_left(), bottom_right()), middle(top_right(), bottom_right()));
        break;
    case Quadrant::bottom_left:
        return CellBox(bottom_left(), center_);
        break;
    default: //to suppress warning
        cout << "quarterBox error: invalid quadrant!" << endl;
        return CellBox();
        break;
    }
}
