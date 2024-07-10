#ifndef qtCWENO_cellBox_H //include guard
#define qtCWENO_cellBox_H

class point
{
public:
    double x, y;
    bool isCloseToStraightLine(double line_coord, Orientation ori); //лежит ли точка внутри полосы шириной 2*(refine padding) около прямой линии
};
point middle(point const& p1, point const& p2); //середина отрезка между точками
double distanceSquared(point const& p1, point const& p2); //квадрат расстояния между точками
double distance(point const& p1, point const& p2); //расстояние между точками

const double ISQ3 = 0.5773502691896258; //1/sqrt(3)
std::pair<point, point> quadraturePoints(point const& p1, point const& p2); //две точки Гауссовой квадратуры на отрезке

static constexpr double inf = std::numeric_limits<double>::infinity();
struct cellBox
{
    point top_left{ -inf,  inf };
    point top_right{ -inf,  -inf };
    point bottom_right{ inf, -inf };
    point bottom_left{ inf,  inf };
    point center = { 0, 0 };

    double size(); //длина стороны
    bool isPointInside(point p); //попадает ли точка в box
    cellBox quarterBox(Quadrant q); //четвертинка box'а
};


#endif