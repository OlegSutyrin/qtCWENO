#ifndef qtCWENO_TreeNode_H //include guard
#define qtCWENO_TreeNode_H

#include "main.h"
#include "CellData.h"
#include "NodeEdge.h"
#include "NodeTag.h"
#include "QuadTree.h"

#include "Eigen\Dense"

//struct теги детей
struct ChildrenTags
{
    std::array<NodeTag, QUADRANTS_NUM> tags = {};
    ChildrenTags(NodeTag tTL, NodeTag tTR, NodeTag tBR, NodeTag tBL) : tags({ tTL, tTR, tBR, tBL }) {};
    const NodeTag& operator()(Quadrant q) const; //() operator: тэг по квадранту
};

//узел дерева
class TreeNode
{
    NodeTag tag_ = {}; //{} вызывает дефолтный конструктор NodeTag
    treeNodeId parentId_ = null; //номер родителя на предыдущем уровне
    treeNodeId childrenId_ = null; //номер первого child на следующем уровне, остальные идут сразу за ним
    bool has_grandchildren = false; //есть ли внуки (нужно для балансировки)
    bool is_leaf = true; //является ли ячейка листом (только для них хранятся физические данные) (избыточно, но проще писать, чем проверять наличие детей)
    bool is_deleted = false; //для отличия живых ячеек от удаленных, т.к. удаление = пометка места на уровне свободным
    cellDataId dataId_ = null; //номер записи для данных в одномерном массиве cell_data, только для листьев
    CellBox box_ = {};
    std::array<NodeTag, NEIGHBOURS_NUM> neighbours = {}; //соседи
    std::array<NodeTag, MAX_NEIGHBOURS12_NUM> neighbours12 = {}; //соседи с учетом диагональных и возможного разделения ребра надвое
    //номера ребер, составляющих стороны ячейки (до DIRECTIONS_NUM*2 штук с учетом возможного разделения ребер)
    //nodeEdgeId edges[DIRECTIONS_NUM * 2] = { null, null, null, null, null, null, null, null }; //нумерация по часовой стрелке, начиная с левой половины верхней стороны
    //значения по умолчанию для edges[8] нужно укзаывать явно, т.к. нет дефолтного конструктора для nodeEdgeId
    double polyCoeffs[EQ_NUM][POLY_COEFF_NUM] = {}; //коэффициенты параболоида CWENO для каждого уравнения

    Eigen::VectorXd coeffs; //вектор неизвестных (px, py, pxx, pxy, pyy) для оптимального полинома
    Eigen::MatrixXd J; //матрица системы
    Eigen::VectorXd rs; //вектор правых частей системы
    Eigen::HouseholderQR<Eigen::MatrixXd> decomp; //решатель
    Eigen::VectorXd coeffsl; //вектор неизвестных (px, py) для линейных функций
    Eigen::MatrixXd Jl[4]; //по одной для каждого подшаблона
    Eigen::VectorXd rsl[4];
    Eigen::HouseholderQR<Eigen::MatrixXd> decompl[4];

public:
    TreeNode() {}; //default constructor
    explicit TreeNode(NodeTag t, treeNodeId p, cellDataId d, CellBox b) : tag_(t), parentId_(p), dataId_(d), box_(b) {}; //конструктор по частичным данным
    
    //accessors
    NodeTag tag() const;
    cellDataId dataId() const;
    CellBox box() const;
    //QuadTree& treeRef(); //ссылка на дерево, содержущее ноду TODO:разобраться, почему не работает (circular dependency?)
    //const QuadTree& treeRefConst() const; //const версия
    NodeTag neighbour(Neighbour n) const; 
    NodeTag neighbour12(Neighbour n12) const; 
    CellData& dataRef(); //ссылка на данные
    const CellData& dataRefConst() const; //const версия
    CellData data() const; //копия данных (с возможным сбором с детей)
    TreeNode& childRef(Quadrant q); //ссылка на ребенка по квадранту
    const TreeNode& childRef(Quadrant q) const; //const версия
    TreeNode& getChildOrSelfByCoords(Point p); //ссылка на ребенка (или себя) по координатам

    //mutators
    void setTag(const NodeTag& t);
    void setDataId(cellDataId id);
    void setBox(const CellBox& b);
    void setNeighbour(Neighbour n, NodeTag t);
    void setNeighbour12(Neighbour12 n12, NodeTag t);
    void setData(const CellData& data);
    void setGrandParency(bool status); //отметка о наличии внуков
    void setChildrenNeighbours(Neighbour n, ChildrenTags tags); //внесение данных о соседях для детей
    int markToRefine(); //пометка ячейки к дроблению
    int refine(); //дробление ячейки
    int tryCoarsen(); //склейка ячейки

    //inspectors
    bool isDeleted() const;
    bool isLeaf() const;
    bool hasChildren() const; //есть ли дети
    bool hasGrandChildren() const; //есть ли внуки
    bool hasNeighbour(Neighbour n) const; //есть ли сосед по направлению
    bool hasNeighbour12(Neighbour12 n12) const;

    //other
    static TreeNode& nodeRef(const NodeTag& tag); //ссылка на ноду по тэгу (static - общая функция для всех нод)
    double magGradRho() const; //примерный градиент плотности

    //output
    std::string dump() const; //дамп ноды в строку
    std::string dumpNeighbourVector(Neighbour n) const; //дамп соседа в виде вектора

    friend class QuadTreeForest; //для глобальных функций TODO:разобраться, как лучше реализовать вложенные циклы без нарушения инкапсуляции
    friend class QuadTree;

    /*
    NodeTag getNodeOrChildTag(int target_depth, Quadrant quadrant) const; //поиск граничной ноды нужного уровня внутри данной (глубина поиска не более 1)
    bool hasEdge(ushorty etype) const; //есть ли ребро
    nodeEdge& getEdge(ushorty etype) const; //ссылка на ребро
    double h() const; //длина стороны ячейки
    void gatherDataFromChildren(); //сбор данных из детей для объединения
    void updateGrandChildren(); //обновление данных о внуках
    //bool hasNeighbour(Neighbour n) const; //есть ли сосед по направлению
    bool hasNeighbour12(Neighbour12 n12) const;
    int neighbours12Num() const; //число соседей
    //bool isBoundary(); //является ли граничной
    //bool isCorner(); //является ли угловой
    TreeNode& getNeigbourOrSelf(Neighbour dir, int target_depth) const; //ссылка на соседа (или на себя, если нет соседа)
    void updateNeighbour(Neighbour dir, NodeTag tag1, NodeTag tag2); //внесение данных о соседе после его дробления
    void updateNeighbour(Neighbour dir, NodeTag tag); //внесение данных о соседе
    int markToRefine(); //пометка ячейки к дроблению
    int refine(); //дробление ячейки
    int tryCoarsen(); //склейка ячейки
    bool isNeighbour12InSubstencil(Neighbour12 n, Quadrant q) const; //попадает ли сосед в подшаблон для линейной функции
    void updateEigenObjects(); //создание или обновление Eigen матриц и т.д. после изменения сетки
    void calcPolynomialCWENO(rkStep rk); //вычисление коэффициентов 2D CWENO полинома
    CellData evalPolynomialAt(point p, rkStep rk = 0); //реконструированное полиномом значение 
    std::string dumpNeighbourVector(Neighbour n) const; //дамп соседа в виде вектора
    std::string dumpNeighbour12Vector(Neighbour12 n) const;*/
};

//struct ссылки на детей
struct ChildrenRefs
{
    TreeNode& rTL;
    TreeNode& rTR;
    TreeNode& rBR;
    TreeNode& rBL;
    ChildrenRefs(TreeNode& _rTL, TreeNode& _rTR, TreeNode& _rBR, TreeNode& _rBL) : rTL(_rTL), rTR(_rTR), rBR(_rBR), rBL(_rTL) {};
    TreeNode& operator()(Quadrant q); //() operator: ссылка по квадранту
};

#endif
