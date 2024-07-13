#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include "json.hpp" //NLohmann json lib

#include "main.h"
#include "output.h"
#include "problemConfig.h"
#include "globalFuncs.h"
#include "cellBox.h"
#include "quadTreeForest.h"


//глобальные объекты
problemConfig config; //конфиг задачи
Globals globals; //глобальные переменные
quadTreeForest forest; //лес деревьев


bool validateBoxAndTrees() //проверка соответствия числа деревьев global box'у
{
    double box_ratio = (config.global_box.bottom_right().x - config.global_box.bottom_left().x) / (config.global_box.top_left().y - config.global_box.bottom_left().y);
    double trees_ratio = (double)config.Nx / (double)config.Ny;
    if (fabs(box_ratio - trees_ratio) > DOUBLE_EPS12)
        return false;
    return true;
}

void addGhostLayer() //добавление слоя ghost-деревьев
{
    double tree_h = (config.global_box.bottom_right().x - config.global_box.bottom_left().x) / config.Nx;
    point new_bottom_left = { config.global_box.bottom_left().x - tree_h, config.global_box.bottom_left().y - tree_h };
    point new_top_right = { config.global_box.top_right().x + tree_h, config.global_box.top_right().y + tree_h };
    config.global_box = cellBox(new_bottom_left, new_top_right);
    config.Nx += 2;
    config.Ny += 2;
    return;
}

//void printNode(treeNode node) { cout << node.dump(); }
//void printTree(quadTree tree) { cout << tree.dump(); }

//------------------------------------- main ------------------------------------
int main(int argc, char** argv)
{
    cout << "Starting...\n";
    auto chrono_start = std::chrono::steady_clock::now();
    globals.file_output.open("output"); //output file

    //чтение конфига из файла
    std::string cfg_filename = (argc < 2 ? "config.json" : argv[1]); //если нет параметра командной строки, читать из config.json
    config = problemConfig(cfg_filename); //TODO: валидация конфига
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
    forest.initialize(); //выделение памяти в лесу

    //генерация стартовых деревьев
    for (quadTreeId id = 0; id < config.Nx * config.Ny; id++)
        forest.addTree(quadTree(id));

    //forest.forAllTrees(printTree, INCLUDE_GHOSTS);
    //forest.forAllNodes(printNode, INCLUDE_GHOSTS, INCLUDE_BRANCHES);

    ExportForest();

    auto chrono_duration = std::chrono::steady_clock::now() - chrono_start;
    int duration_seconds = (int)round(std::chrono::duration<double, std::milli>(chrono_duration).count() / 1000);
    cout << "Job's done in " << std::to_string(duration_seconds) << " seconds." << endl;
    return 0;
}