#ifndef qtCWENO_CellBox_H //include guard
#define qtCWENO_CellBox_H

//struct точка на плоскости
struct Point
{
    double x, y;
    bool isCloseToStraightLine(double line_coord, Orientation ori) const; //лежит ли точка внутри полосы шириной 2*(refine padding) около кардинальной линии
    bool isInsideWedge(double angle_bottom, double angle_top); //попадает ли точка нутрь клина, углы в градусах

    friend std::ostream& operator<<(std::ostream& os, const Point& p); //output overload
};
Point middle(Point const& p1, Point const& p2); //середина отрезка между точками
double distanceSquared(Point const& p1, Point const& p2); //квадрат расстояния между точками
double distance(Point const& p1, Point const& p2); //расстояние между точками

const double ISQ3 = 0.5773502691896258; //1/sqrt(3)
std::array<Point, QUADRATURE_POINTS_NUM> quadraturePoints(Point const& p1, Point const& p2); //две точки Гауссовой квадратуры на отрезке

class CellBox
{
    Point vertices[QUADRANTS_NUM] = { { -1.0, 1.0 }, { 1.0, 1.0 }, { 1.0, -1.0 }, { -1.0, -1.0 } }; //top_left, top_right, bottom_right, bottom_left
    Point center_ = { 0, 0 };

public:
    CellBox() {}; //дефолтный конструктор
    CellBox(Point bl, Point tr); //конструктор по bottom_left и top_right точкам
    //point getP(Quadrant q) const;
    Point top_left() const; //получение вершин box'а
    Point top_right() const;
    Point bottom_right() const;
    Point bottom_left() const;
    Point center() const; //получение центра
    //void setP(Quadrant q, point pt); //задание вершины
    double size() const;
    void updateCenter(); //вычисление центра по угловым точкам
    bool isPointInside(Point p) const; //попадает ли точка в box
    bool intersectLineStraight(double line_coord, Orientation ori) const; //пересекает ли ячейку прямая линия (горизонтальная или вертикальная)
    bool intersectLineSlanted(double k, double b) const; //пересекает ли ячейку прямая линия с уравнением kx+b
    bool intersectLineEllipse(double axle_x, double axle_y) const; //пересекает ли ячейку эллипс
    CellBox quarterBox(Quadrant q); //получение четвертинки box'а
};


#endif