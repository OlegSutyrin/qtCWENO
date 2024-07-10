#ifndef qtCWENO_cellBox_H //include guard
#define qtCWENO_cellBox_H

class point
{
public:
    double x, y;
    bool isCloseToStraightLine(double line_coord, Orientation ori); //����� �� ����� ������ ������ ������� 2*(refine padding) ����� ������ �����
};
point middle(point const& p1, point const& p2); //�������� ������� ����� �������
double distanceSquared(point const& p1, point const& p2); //������� ���������� ����� �������
double distance(point const& p1, point const& p2); //���������� ����� �������

const double ISQ3 = 0.5773502691896258; //1/sqrt(3)
std::pair<point, point> quadraturePoints(point const& p1, point const& p2); //��� ����� ��������� ���������� �� �������

static constexpr double inf = std::numeric_limits<double>::infinity();
struct cellBox
{
    point top_left{ -inf,  inf };
    point top_right{ -inf,  -inf };
    point bottom_right{ inf, -inf };
    point bottom_left{ inf,  inf };
    point center = { 0, 0 };

    double size(); //����� �������
    bool isPointInside(point p); //�������� �� ����� � box
    cellBox quarterBox(Quadrant q); //����������� box'�
};


#endif