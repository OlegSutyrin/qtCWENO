#ifndef qtCWENO_CellBox_H //include guard
#define qtCWENO_CellBox_H

//struct ����� �� ���������
struct Point
{
    double x, y;
    bool isCloseToStraightLine(double line_coord, Orientation ori) const; //����� �� ����� ������ ������ ������� 2*(refine padding) ����� ������������ �����
    bool isInsideWedge(double angle_bottom, double angle_top); //�������� �� ����� ����� �����, ���� � ��������

    friend std::ostream& operator<<(std::ostream& os, const Point& p); //output overload
};
Point middle(Point const& p1, Point const& p2); //�������� ������� ����� �������
double distanceSquared(Point const& p1, Point const& p2); //������� ���������� ����� �������
double distance(Point const& p1, Point const& p2); //���������� ����� �������

const double ISQ3 = 0.5773502691896258; //1/sqrt(3)
std::array<Point, QUADRATURE_POINTS_NUM> quadraturePoints(Point const& p1, Point const& p2); //��� ����� ��������� ���������� �� �������

class CellBox
{
    Point vertices[QUADRANTS_NUM] = { { -1.0, 1.0 }, { 1.0, 1.0 }, { 1.0, -1.0 }, { -1.0, -1.0 } }; //top_left, top_right, bottom_right, bottom_left
    Point center_ = { 0, 0 };

public:
    CellBox() {}; //��������� �����������
    CellBox(Point bl, Point tr); //����������� �� bottom_left � top_right ������
    //point getP(Quadrant q) const;
    Point top_left() const; //��������� ������ box'�
    Point top_right() const;
    Point bottom_right() const;
    Point bottom_left() const;
    Point center() const; //��������� ������
    //void setP(Quadrant q, point pt); //������� �������
    double size() const;
    void updateCenter(); //���������� ������ �� ������� ������
    bool isPointInside(Point p) const; //�������� �� ����� � box
    bool intersectLineStraight(double line_coord, Orientation ori) const; //���������� �� ������ ������ ����� (�������������� ��� ������������)
    bool intersectLineSlanted(double k, double b) const; //���������� �� ������ ������ ����� � ���������� kx+b
    bool intersectLineEllipse(double axle_x, double axle_y) const; //���������� �� ������ ������
    CellBox quarterBox(Quadrant q); //��������� ����������� box'�
};


#endif