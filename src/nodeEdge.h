#ifndef qtCWENO_nodeEdge_H //include guard
#define qtCWENO_nodeEdge_H

#include "main.h"
#include "NodeTag.h"

//ребро - сторона (или еЄ половина) €чейки, дл€ хранени€ потоков
class NodeEdge {
    bool is_deleted = false;
    Orientation orientation_ = Orientation::horizontal;
    NodeTag n1_{}, n2_{}; //тэги €чеек-соседей, упор€дочены "верхн€€, нижн€€" или "лева€, права€"
    std::array<double, EQ_NUM> FQ_ = {}; //поток через сторону
    std::array<Point, QUADRATURE_POINTS_NUM> qps = {}; //точки дл€ квадратуры, упор€дочены "нижн€€, верхн€€" или "лева€, права€"

public:
    NodeEdge() {}; //default contructor
    NodeEdge(NodeTag n1, NodeTag n2, Orientation orientation) : n1_(n1), n2_(n2), orientation_(orientation) {}; //конструктор по двум €чейкам
    
    //accessors
    const NodeTag n1() const;
    const NodeTag n2() const;
    double FQ(Equation eq) const; //компонента потока
    Point middle(); //центр ребра
    double length() const; //вычисление длины ребра по точкам квадратуры TODO: хранить вершины и длину в самом ребре

    //mutators
    void markDeleted();
    void setN1(NodeTag ntag);
    void setN2(NodeTag ntag);
    
    //inspectors
    bool isDeleted();

    //other
    void computeQuadraturePoints(); //расчет точек квадратуры
    void computeFluxLF(rkStep rk); //расчет потока (Lax-Friedrich flux)
    void computeFluxHLLC(rkStep rk); //расчет потока (HLLC Rieman solver)
    void computeFluxRiemannExact(rkStep rk); //расчет потока (exact Rieman solver)

    bool operator==(const NodeEdge& rhs); //equality overload

    //output
    std::string dumpNeigboursVectors(); //дамп соседей в виде векторов
    std::string dumpQuadraturePoints(); //дамп точек квадратуры
    std::string dumpFlux(); //дамп значений потоков

    friend std::ostream& operator<<(std::ostream& os, const NodeEdge& edge); //output overload
};

#endif
