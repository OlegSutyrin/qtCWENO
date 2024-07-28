#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include "json.hpp" //NLohmann json lib

#include "main.h"
#include "output.h"
#include "ProblemConfig.h"
#include "globalFuncs.h"
#include "CellBox.h"
#include "QuadTreeForest.h"

//функции преобразования enum class'ов
Neighbour opposite(Neighbour n) //противоположный сосед
{
    switch (n)
    {
    case Neighbour::top: return Neighbour::bottom;
    case Neighbour::right: return Neighbour::left;
    case Neighbour::bottom: return Neighbour::top;
    case Neighbour::left: return Neighbour::right;
    default: return Neighbour::top; //to suppress warning
    }
}
Neighbour12 opposite(Neighbour12 n12) //противоположный (по ребру или вершине) сосед
{
    switch (n12)
    {
    case Neighbour12::top1: return Neighbour12::bottom2;
    case Neighbour12::top2: return Neighbour12::bottom1;
    case Neighbour12::top_right: return Neighbour12::bottom_left;
    case Neighbour12::right1: return Neighbour12::left2;
    case Neighbour12::right2: return Neighbour12::left1;
    case Neighbour12::bottom_right: return Neighbour12::top_left;
    case Neighbour12::bottom1: return Neighbour12::top2;
    case Neighbour12::bottom2: return Neighbour12::top1;
    case Neighbour12::bottom_left: return Neighbour12::top_right;
    case Neighbour12::left1: return Neighbour12::right2;
    case Neighbour12::left2: return Neighbour12::right1;
    case Neighbour12::top_left: return Neighbour12::bottom_right;
    default: return Neighbour12::top1; //to suppress warning
    }
}
Quadrant toQuadrant(Neighbour12 n12) //квадрант по соседу12
{
    switch (n12)
    {
    case Neighbour12::top_right: return Quadrant::top_right;
    case Neighbour12::bottom_right: return Quadrant::bottom_right;
    case Neighbour12::bottom_left: return Quadrant::bottom_left;
    case Neighbour12::top_left: return Quadrant::top_left;
    default: return Quadrant::top_right; //to suppress warning
    }
}
Neighbour12 toNeighbour12(Neighbour n) //сосед12 по соседу
{
    switch (n)
    {
    case Neighbour::top: return Neighbour12::top1;
    case Neighbour::right: return Neighbour12::right1;
    case Neighbour::bottom: return Neighbour12::bottom1;
    case Neighbour::left: return Neighbour12::left1;
    default: return Neighbour12::top1; //to suppress warning
    }
}
Edge toEdge(Neighbour n) //ребро по соседу4
{
    switch (n)
    {
    case Neighbour::top: return Edge::top1;
    case Neighbour::right: return Edge::right1;
    case Neighbour::bottom: return Edge::bottom1;
    case Neighbour::left: return Edge::left1;
    default: return Edge::top1; //to suppress warning
    }
}
Edge next(Edge e) //следующее ребро
{
    return static_cast<Edge>(static_cast<int>(e) + 1);
}
Edge opposite(Edge e) //противоположное (по стороне) ребро
{
    switch (e)
    {
    case Edge::top1: return Edge::bottom2;
    case Edge::top2: return Edge::bottom1;
    case Edge::right1: return Edge::left2;
    case Edge::right2: return Edge::left1;
    case Edge::bottom1: return Edge::top2;
    case Edge::bottom2: return Edge::top1;
    case Edge::left1: return Edge::right2;
    case Edge::left2: return Edge::right1;
    default: return Edge::top1; //to suppress warning
    }
}

//глобальные объекты
ProblemConfig config; //конфиг задачи
Globals globals; //глобальные переменные
QuadTreeForest forest; //лес деревьев

//"static" limits function scope to this translation unit
static bool validateBoxAndTrees() //проверка соответствия числа деревьев global box'у
{
    double box_ratio = (config.global_box.bottom_right().x - config.global_box.bottom_left().x) / (config.global_box.top_left().y - config.global_box.bottom_left().y);
    double trees_ratio = (double)config.Nx / (double)config.Ny;
    if (fabs(box_ratio - trees_ratio) > DOUBLE_EPS12)
        return false;
    return true;
}

static void addGhostLayer() //добавление слоя ghost-деревьев
{
    double tree_h = (config.global_box.bottom_right().x - config.global_box.bottom_left().x) / config.Nx;
    Point new_bottom_left = { config.global_box.bottom_left().x - tree_h, config.global_box.bottom_left().y - tree_h };
    Point new_top_right = { config.global_box.top_right().x + tree_h, config.global_box.top_right().y + tree_h };
    config.global_box = CellBox(new_bottom_left, new_top_right);
    config.Nx += 2;
    config.Ny += 2;
    return;
}

bool checkExportTime() //проверка, не пора ли выводить в файл после данного шага по времени
{
    double n = 0.0;
    double f1 = modf(globals.time / config.exportPeriod, &n); //signed дробная часть
    double f2 = modf((globals.time + globals.dt) / config.exportPeriod, &n);
    if (f2 < f1) //переход через целое число 
        return true;
    else
        return false;
}


//void printNode(TreeNode node) { cout << node.dump(); }
//void printTree(QuadTree tree) { cout << tree.dump(); }

//------------------------------------- main ------------------------------------
int main(int argc, char** argv)
{
    cout << "Starting...\n";
    auto chrono_start = std::chrono::steady_clock::now();
    globals.file_output.open("output"); //output file

    //чтение конфига из файла
    std::string cfg_filename = (argc < 2 ? "config.json" : argv[1]); //если нет параметра командной строки, читать из config.json
    config = ProblemConfig(cfg_filename); //TODO: валидация конфига
    cout << config;
    globals.file_output << config;
    //std::cin.get(); //пауза

    //проверка соответствия global box и числа деревьев
    if (!validateBoxAndTrees())
    {
        cout << "bad trees configuration for global box!" << endl;
        globals.file_output << "bad trees configuration for global box!" << endl;
        std::cin.get(); //пауза
        return 0;
    }

    addGhostLayer(); //добавление слоя ghost-деревьев (изменение config'а)
    //задание уровней magGradRho
    cout << "refine levels: ";
    double lvl = config.meshRefinelevel0;
    for (size_t i = 0; i <= config.max_depth; i++)
    {
        globals.refineLvls.push_back(lvl);
        cout << "(" << i << ": " << globals.refineLvls[i] << ") ";
        lvl *= config.meshRefinefactor;
    }
    cout << endl;

    //создание леса
    forest.initialize(); //выделение памяти
    for (quadTreeId id = 0; id < config.Nx * config.Ny; id++) //генерация стартовых деревьев
        forest.addTree(QuadTree(id));

    //forest.forAllTrees(printTree, INCLUDE_GHOSTS);
    //forest.forAllNodes(printNode, INCLUDE_GHOSTS, INCLUDE_BRANCHES);

    forest.meshRefineInitial();
    forest.initialCondition();
    ExportForest();


    //--------------------- цикл по времени -----------------------------
    /*const int steps = 20;
    double dx = config.bubble_axle_x / (double)steps;
    //double dy = config.bubble_axle_y / (double)steps;
    for (int step = 0; step < steps; step++)
    {
        globals.time += 0.5 / (double)steps;
        config.bubble_axle_x -= dx;
        forest.meshCoarsenInitial();
        forest.meshRefineInitial();
        ExportForest();
    }
    for (int step = 0; step < steps; step++)
    {
        globals.time += 0.5 / (double)steps;
        config.bubble_axle_x += dx;
        forest.meshCoarsenInitial();
        forest.meshRefineInitial();
        ExportForest();
    }*/
    while (globals.time < config.end_time)
    {
        globals.dt = forest.CFLTimestepSize(); //TODO: адаптировать под FV + Riemann solver (брать скорости волн оттуда)?
        //globals.dt = AdjustTimeStepForExport();
        cout << "Step " << (globals.timestep_number + 1) << ": " << globals.time << " --> " << (globals.time + globals.dt) << " (dt = " << globals.dt << ")";
        globals.file_output << "Step " << (globals.timestep_number + 1) << ": " << globals.time << " --> " << (globals.time + globals.dt) << " (dt = " << globals.dt << ")";

        forest.advanceTime(); //шаг по времени
        forest.putQn(config.rk_order); //переброс Qn[rk_order] --> Qn[0] после последнего шага Рунге-Кутты
        forest.boundaryConditions(0);

        bool to_export = checkExportTime(); //нужно ли будет выводить после этого шага
        globals.time += globals.dt;
        globals.timestep_number++;

        if (!config.meshRefineAll && globals.timestep_number % config.meshRefinePeriod == 0)
        {
            forest.meshUpdate(); //TODO: разобраться, нужно ли делать перерасчет средних значений в соседних ячейках (Semplice2015, p.14)
            cout << " (remeshed to " << forest.activeNodesNumber() << " nodes, " << forest.leavesNumber() << " leaves)";
            globals.file_output << " (remeshed to " << forest.activeNodesNumber() << " nodes, " << forest.leavesNumber() << " leaves)";
            //ExportForest();
        }
        cout << endl;
        globals.file_output << endl;

        if (true || to_export)
            ExportForest();
    } //конец цикла по времени


    auto chrono_duration = std::chrono::steady_clock::now() - chrono_start;
    int duration_seconds = (int)round(std::chrono::duration<double, std::milli>(chrono_duration).count() / 1000);
    cout << "Job's done in " << std::to_string(duration_seconds) << " seconds." << endl;
    return 0;
}