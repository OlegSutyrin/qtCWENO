#ifndef qtCWENO_treeNode_H //include guard
#define qtCWENO_treeNode_H

#include "main.h"
#include "cellData.h"
#include "nodeEdge.h"
#include "nodeTag.h"

#include "Eigen\Dense"

//узел дерева
class treeNode
{
    nodeTag tag_ = {}; //{} вызывает дефолтный конструктор nodeTag
    treeNodeId parentId_ = null; //номер родителя на предыдущем уровне
    treeNodeId childrenId_ = null; //номер первого child на следующем уровне, остальные идут сразу за ним
    bool has_grandchildren = false; //есть ли внуки (нужно для балансировки)
    bool is_leaf = true; //является ли ячейка листом (только для них хранятся физические данные) (избыточно, но проще писать, чем проверять наличие детей)
    bool is_deleted = false; //для отличия живых ячеек от удаленных, т.к. удаление = пометка места на уровне свободным
    cellDataId dataId_ = null; //номер записи для данных в одномерном массиве cell_data, только для листьев
    cellBox box_ = {};
    std::array<nodeTag, NEIGHBOURS_NUM> neighbours = {}; //соседи
    std::array<nodeTag, MAX_NEIGHBOURS12_NUM> neighbours12 = {}; //соседи с учетом диагональных и возможного разделения ребра надвое
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
    treeNode() {}; //default constructor
    treeNode(nodeTag t, treeNodeId p, cellDataId d, cellBox b) : tag_(t), parentId_(p), dataId_(d), box_(b) {}; //конструктор по частичным данным
    bool isDeleted() const;
    bool isLeaf() const;
    nodeTag tag() const; //получение тэга 
    void setTag(nodeTag t); //задание тэга
    cellDataId dataId() const;
    void setDataId(cellDataId id);
    cellBox box() const;
    void setBox(cellBox b);
    nodeTag getNeighbour(Neighbour n) const; 
    void setNeighbour(Neighbour n, nodeTag t);
    nodeTag getNeighbour12(Neighbour n12) const; 
    void setNeighbour12(Neighbour12 n12, nodeTag t);
    void setData(cellData data); //запись данных

    cellData& dataRef() const; //ссылка на данные
    cellData data() const; //копия данных
    treeNode& childRef(Quadrant q); //ссылка на ребенка по квадранту
    const treeNode& childRef(Quadrant q) const; //(const) ссылка на ребенка
    treeNode& getChildByCoords(point p); //ссылка на ребенка по координатам

    
    
    static treeNode& nodeRef(nodeTag tag); //ссылка на ноду по тэгу (static - общая функция для всех нод)

    double magGradRho() const; //примерный градиент плотности

    std::string dump() const; //дамп ноды в строку



    friend class quadTree;
    friend class quadTreeForest; //для глобальных функций TODO:разобраться, как лучше реализовать вложенные циклы без нарушения инкапсуляции



    /*static treeNode& getNodeByCoords(point p); //ссылка на ноду по координатам
    treeNode& getChildByCoords(point p) const; //ссылка на ребенка по координатам
    treeNode& getChild(Quadrant q) const; //ссылка на ребенка по квадранту
    nodeTag getNodeOrChildTag(int target_depth, Quadrant quadrant) const; //поиск граничной ноды нужного уровня внутри данной (глубина поиска не более 1)
    //bool hasEdge(ushorty etype) const; //есть ли ребро
    //nodeEdge& getEdge(ushorty etype) const; //ссылка на ребро
    double h() const; //длина стороны ячейки
    bool hasChildren() const; //есть ли дети
    bool hasGrandChildren() const; //есть ли внуки
    void gatherDataFromChildren(); //сбор данных из детей для объединения
    void updateGrandChildren(); //обновление данных о внуках
    //bool hasNeighbour(Neighbour n) const; //есть ли сосед по направлению
    bool hasNeighbour12(Neighbour12 n12) const;
    int neighbours12Num() const; //число соседей
    //bool isBoundary(); //является ли граничной
    //bool isCorner(); //является ли угловой
    treeNode& getNeigbourOrSelf(Neighbour dir, int target_depth) const; //ссылка на соседа (или на себя, если нет соседа)
    void updateNeighbour(Neighbour dir, nodeTag tag1, nodeTag tag2); //внесение данных о соседе после его дробления
    void updateNeighbour(Neighbour dir, nodeTag tag); //внесение данных о соседе
    int markToRefine(); //пометка ячейки к дроблению
    int refine(); //дробление ячейки
    int tryCoarsen(); //склейка ячейки
    double magGradRho() const; //примерный градиент плотности
    bool isNeighbour12InSubstencil(Neighbour12 n, Quadrant q) const; //попадает ли сосед в подшаблон для линейной функции
    void updateEigenObjects(); //создание или обновление Eigen матриц и т.д. после изменения сетки
    void calcPolynomialCWENO(rkStep rk); //вычисление коэффициентов 2D CWENO полинома
    cellData evalPolynomialAt(point p, rkStep rk = 0); //реконструированное полиномом значение 
    std::string dump() const; //дамп ноды в строку
    std::string dumpNeighbourVector(Neighbour n) const; //дамп соседа в виде вектора
    std::string dumpNeighbour12Vector(Neighbour12 n) const;*/
};


#endif
