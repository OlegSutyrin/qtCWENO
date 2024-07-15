#ifndef qtCWENO_output_H //include guard
#define qtCWENO_output_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <fstream>
#include <iostream>

#include "ProblemConfig.h"

struct dataExtrema //массивы мин и макс значений для Tecplot'а
{
    double minima[TECPLOT_FIELDS_NUMBER];
    double maxima[TECPLOT_FIELDS_NUMBER];
};

double NaNcleared(double f); //если f = NAN, возвращает -999
void ExportForest(); //вывод всего леса (кроме ghost'ов)




#endif


