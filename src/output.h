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

struct dataExtrema //������� ��� � ���� �������� ��� Tecplot'�
{
    double minima[TECPLOT_FIELDS_NUMBER];
    double maxima[TECPLOT_FIELDS_NUMBER];
};

double NaNcleared(double f); //���� f = NAN, ���������� -999
void ExportForest(); //����� ����� ���� (����� ghost'��)




#endif


