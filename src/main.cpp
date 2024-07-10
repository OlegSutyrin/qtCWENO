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
//Globals globals; //глобальные переменные
//quadTreeForest forest; //лес деревьев


//------------------------------------- main ------------------------------------
int main(int argc, char** argv)
{
    cout << "Starting...\n";
    auto chrono_start = std::chrono::steady_clock::now();

    //std::ofstream file_output("output");

    //чтение конфига из файла
    std::ifstream input;
    std::string cfg_filename;
    if (argc < 2) //если нет параметра командной строки
        cfg_filename = "config.json";
    else
        cfg_filename = argv[1];
    config = problemConfig(cfg_filename);
    cout << config;
    std::cin.get(); //пауза




    auto chrono_duration = std::chrono::steady_clock::now() - chrono_start;
    int duration_seconds = (int)round(std::chrono::duration<double, std::milli>(chrono_duration).count() / 1000);
    cout << "Job's done in " << std::to_string(duration_seconds) << " seconds." << endl;
    return 0;
}