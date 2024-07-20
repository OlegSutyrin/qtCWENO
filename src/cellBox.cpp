#include <utility> //for std::pair

#include "main.h"
#include "CellBox.h"

//------------------------------------------------------------------- Point -------------------------------------------------------------------
bool Point::isCloseToStraightLine(double line_coord, Orientation ori) const //лежит ли точка внутри полосы шириной 2*(refine padding) около прямой линии
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

Point middle(Point const& p1, Point const& p2) //середина отрезка между точками
{
    return { (p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0 };
}
double distanceSquared(Point const& p1, Point const& p2) //квадрат расстояния между точками
{
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}
double distance(Point const& p1, Point const& p2) //расстояние между точками
{
    return sqrt(distanceSquared(p1, p2));
}

std::array<Point, QUADRATURE_POINTS_NUM> quadraturePoints(Point const& p1, Point const& p2) //две точки Гауссовой квадратуры на отрезке
{
    std::array<Point, QUADRATURE_POINTS_NUM> ret;
    double dx = fabs(p2.x - p1.x);
    double dy = fabs(p2.y - p1.y);
    ret[0] = {middle(p1, p2).x - ISQ3 * 0.5 * dx, middle(p1, p2).y - ISQ3 * 0.5 * dy};
    ret[1] = {middle(p1, p2).x + ISQ3 * 0.5 * dx, middle(p1, p2).y + ISQ3 * 0.5 * dy};
    return ret;
}

//------------------------------------------------------------------- CellBox -------------------------------------------------------------------
CellBox::CellBox(Point bl, Point tr) //конструктор по bottom_left и top_right точкам
{
    vertices[static_cast<int>(Quadrant::top_left)] = { bl.x, tr.y };
    vertices[static_cast<int>(Quadrant::top_right)] = tr;
    vertices[static_cast<int>(Quadrant::bottom_right)] = { tr.x, bl.y };
    vertices[static_cast<int>(Quadrant::bottom_left)] = bl;
    updateCenter();
}

//point CellBox::getP(Quadrant q) const { return vertices[static_cast<int>(q)]; } //получение вершин box'а
Point CellBox::top_left() const { return vertices[static_cast<int>(Quadrant::top_left)]; }
Point CellBox::top_right() const { return vertices[static_cast<int>(Quadrant::top_right)]; }
Point CellBox::bottom_right() const { return vertices[static_cast<int>(Quadrant::bottom_right)]; }
Point CellBox::bottom_left() const { return vertices[static_cast<int>(Quadrant::bottom_left)]; }
Point CellBox::center() const { return center_; } //получение центра

//void CellBox::setP(Quadrant q, point pt) { vertices[static_cast<int>(q)] = pt; } //задание вершины

double CellBox::size() const { return bottom_right().x - bottom_left().x; }

void CellBox::updateCenter() //вычисление центра по угловым точкам
{ 
    center_ = middle(vertices[static_cast<int>(Quadrant::bottom_left)], vertices[static_cast<int>(Quadrant::top_right)]); 
}

bool CellBox::isPointInside(Point p) const //попадает ли точка в box (включая левую и ниюжнюю границы)
{
    if (p.x >= vertices[static_cast<int>(Quadrant::bottom_left)].x && p.x < vertices[static_cast<int>(Quadrant::bottom_right)].x &&
        p.y >= vertices[static_cast<int>(Quadrant::bottom_left)].y && p.y < vertices[static_cast<int>(Quadrant::top_left)].y)
        return true;
    return false;
}

bool CellBox::intersectLineStraight(double line_coord, Orientation ori) const //пересекает ли ячейку прямая линия (горизонтальная или вертикальная)
{
    if (ori == Orientation::horizontal)
    {
        if ((bottom_left().y - line_coord) * (top_right().y - line_coord) < 0)
            return true;
    }
    else
    {
        if ((bottom_left().x - line_coord) * (top_right().x - line_coord) < 0)
            return true;
    }
    return false;
}

//inline makes compiler to substitute formula in place of function call
static inline double lsfEllipse(Point p, double axle_x, double axle_y) //level-set function for ellipse
{
    return p.x * p.x / axle_x / axle_x + p.y * p.y / axle_y / axle_y - 1.0;
}
bool CellBox::intersectLineEllipse(double axle_x, double axle_y) const //пересекает ли ячейку эллипс
{
    if (lsfEllipse(bottom_left(), axle_x, axle_y) * lsfEllipse(top_right(), axle_x, axle_y) < 0 ||
        lsfEllipse(top_left(), axle_x, axle_y) * lsfEllipse(bottom_right(), axle_x, axle_y) < 0)
    {
        return true;
    }
    return false;
}

CellBox CellBox::quarterBox(Quadrant q) //получение четвертинки box'а
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
