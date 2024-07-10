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
//quadTreeForest forest; //лес деревьев


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




    auto chrono_duration = std::chrono::steady_clock::now() - chrono_start;
    int duration_seconds = (int)round(std::chrono::duration<double, std::milli>(chrono_duration).count() / 1000);
    cout << "Job's done in " << std::to_string(duration_seconds) << " seconds." << endl;
    return 0;
}