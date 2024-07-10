#include <utility> //for std::pair

#include "main.h"
#include "cellBox.h"

point middle(point const& p1, point const& p2) //середина отрезка между точками
{
    return { (p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0 };
}
double distanceSquared(point const& p1, point const& p2) //квадрат расстояния между точками
{
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}
double distance(point const& p1, point const& p2) //расстояние между точками
{
    return sqrt(distanceSquared(p1, p2));
}
bool point::isCloseToStraightLine(double line_coord, Orientation ori) //лежит ли точка внутри полосы шириной 2*(refine padding) около прямой линии
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
