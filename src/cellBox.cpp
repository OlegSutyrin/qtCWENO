#include <utility> //for std::pair

#include "main.h"
#include "cellBox.h"

bool point::isCloseToStraightLine(double line_coord, Orientation ori) //����� �� ����� ������ ������ ������� 2*(refine padding) ����� ������ �����
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
std::ostream& operator<<(std::ostream& os, const point& p) //output overload
{
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}

point middle(point const& p1, point const& p2) //�������� ������� ����� �������
{
    return { (p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0 };
}
double distanceSquared(point const& p1, point const& p2) //������� ���������� ����� �������
{
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}
double distance(point const& p1, point const& p2) //���������� ����� �������
{
    return sqrt(distanceSquared(p1, p2));
}

std::pair<point, point> quadraturePoints(point const& p1, point const& p2) //��� ����� ��������� ���������� �� �������
{
    std::pair<point, point> ret;
    double dx = fabs(p2.x - p1.x);
    double dy = fabs(p2.y - p1.y);
    ret.first  = { middle(p1, p2).x - ISQ3 * 0.5 * dx, middle(p1, p2).y - ISQ3 * 0.5 * dy };
    ret.second = { middle(p1, p2).x + ISQ3 * 0.5 * dx, middle(p1, p2).y + ISQ3 * 0.5 * dy };
    return ret;
}

point cellBox::getP(Quadrant q) const { return vertices[static_cast<int>(q)]; } //��������� ������ box'�
point cellBox::top_left() const { return vertices[static_cast<int>(Quadrant::top_left)]; }
point cellBox::top_right() const { return vertices[static_cast<int>(Quadrant::top_right)]; }
point cellBox::bottom_right() const { return vertices[static_cast<int>(Quadrant::bottom_right)]; }
point cellBox::bottom_left() const { return vertices[static_cast<int>(Quadrant::bottom_left)]; }

void cellBox::setP(Quadrant q, point pt) //������� �������
{
    vertices[static_cast<int>(q)] = pt;
}
point cellBox::getCenter() const //��������� ������
{
    return center;
}

void cellBox::updateCenter()
{
    center = middle(vertices[static_cast<int>(Quadrant::bottom_left)], vertices[static_cast<int>(Quadrant::top_right)]);
}

bool  cellBox::isPointInside(point p) //�������� �� ����� � box (������� ����� � ������� �������)
{
    if (p.x >= vertices[static_cast<int>(Quadrant::bottom_left)].x && p.x < vertices[static_cast<int>(Quadrant::bottom_right)].x &&
        p.y >= vertices[static_cast<int>(Quadrant::bottom_left)].y && p.y < vertices[static_cast<int>(Quadrant::top_left)].y)
        return true;
    return false;
}

cellBox cellBox::quarterBox(Quadrant q) //��������� ����������� box'�
{
    cellBox ret;
    switch (q)
    {
    case Quadrant::top_left:
        ret.setP(Quadrant::top_left, top_left());
        ret.setP(Quadrant::top_right, middle(top_left(), top_right()));
        ret.setP(Quadrant::bottom_right, center);
        ret.setP(Quadrant::bottom_left, middle(top_left(), bottom_left()));
        break;
    case Quadrant::top_right:
        ret.setP(Quadrant::top_left, middle(top_left(), top_right()));
        ret.setP(Quadrant::top_right, top_right());
        ret.setP(Quadrant::bottom_right, middle(top_right(), bottom_right()));
        ret.setP(Quadrant::bottom_left, center);
        break;
    case Quadrant::bottom_right:
        ret.setP(Quadrant::top_left, center);
        ret.setP(Quadrant::top_right, middle(top_right(), bottom_right()));
        ret.setP(Quadrant::bottom_right, bottom_right());
        ret.setP(Quadrant::bottom_left, middle(bottom_left(), bottom_right()));
        break;
    case Quadrant::bottom_left:
        ret.setP(Quadrant::top_left, middle(bottom_left(), top_left()));
        ret.setP(Quadrant::top_right, center);
        ret.setP(Quadrant::bottom_right, middle(bottom_left(), bottom_right()));
        ret.setP(Quadrant::bottom_left, bottom_left());
        break;
    }
    ret.updateCenter();
    return ret;
}
