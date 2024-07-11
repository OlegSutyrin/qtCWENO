#ifndef qtCWENO_cellBox_H //include guard
#define qtCWENO_cellBox_H

class point
{
public:
    double x, y;
    bool isCloseToStraightLine(double line_coord, Orientation ori) const; //лежит ли точка внутри полосы шириной 2*(refine padding) около прямой линии
    friend std::ostream& operator<<(std::ostream& os, const point& p); //output overload
};
point middle(point const& p1, point const& p2); //середина отрезка между точками
double distanceSquared(point const& p1, point const& p2); //квадрат расстояния между точками
double distance(point const& p1, point const& p2); //расстояние между точками

const double ISQ3 = 0.5773502691896258; //1/sqrt(3)
std::pair<point, point> quadraturePoints(point const& p1, point const& p2); //две точки Гауссовой квадратуры на отрезке

class cellBox
{
    point vertices[QUADRANTS_NUM] = { { -1.0, 1.0 } ,{ 1.0, 1.0 } ,{ 1.0, -1.0 } ,{ -1.0, -1.0 } }; //top_left, top_right, bottom_right, bottom_left
    point center_ = { 0, 0 };

public:
    cellBox() {}; //дефолтный конструктор
    cellBox(point bl, point tr); //конструктор по bottom_left и top_right точкам
    //point getP(Quadrant q) const;
    point top_left() const; //получение вершин box'а
    point top_right() const;
    point bottom_right() const;
    point bottom_left() const;
    point center() const; //получение центра
    //void setP(Quadrant q, point pt); //задание вершины
    void updateCenter(); //вычисление центра по угловым точкам
    bool isPointInside(point p) const; //попадает ли точка в box
    cellBox quarterBox(Quadrant q); //получение четвертинки box'а
};


#endif