#ifndef qtCWENO_cellBox_H //include guard
#define qtCWENO_cellBox_H

class point
{
public:
    double x, y;
    bool isCloseToStraightLine(double line_coord, Orientation ori); //лежит ли точка внутри полосы шириной 2*(refine padding) около прямой линии
    friend std::ostream& operator<<(std::ostream& os, const point& p); //output overload
};
point middle(point const& p1, point const& p2); //середина отрезка между точками
double distanceSquared(point const& p1, point const& p2); //квадрат расстояния между точками
double distance(point const& p1, point const& p2); //расстояние между точками

const double ISQ3 = 0.5773502691896258; //1/sqrt(3)
std::pair<point, point> quadraturePoints(point const& p1, point const& p2); //две точки Гауссовой квадратуры на отрезке

static constexpr double inf = std::numeric_limits<double>::infinity();
class cellBox
{
    point vertices[QUADRANTS_NUM] = { { -inf, inf } ,{ -inf, -inf } ,{ inf, -inf } ,{ inf, inf } }; //top_left, top_right, bottom_right, bottom_left
    point center = { 0, 0 };

public:
    point getP(Quadrant q) const; //получение вершин box'а
    point top_left() const;
    point top_right() const;
    point bottom_right() const;
    point bottom_left() const;
    void setP(Quadrant q, point pt); //задание вершины
    point getCenter() const; //получение центра
    void updateCenter(); //обновление центра по вершинам
    bool isPointInside(point p); //попадает ли точка в box
    cellBox quarterBox(Quadrant q); //получение четвертинки box'а
};


#endif