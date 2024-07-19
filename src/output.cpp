#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>

#include "main.h"
#include "CellBox.h"
#include "TreeNode.h"
#include "QuadTreeForest.h"
//#include "globalFuncs.h"

double NaNcleared(double f)  //если f = NAN, возвращает -999
{
    if (isnan(f))
        return -999.0;
    else
        return f;
}

/*void exportForestNeighbours(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;
    file_output << "ZONE T = \"Forest neighbours\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& tree : forest.trees)
    {
        for (auto& nodes_level : tree.nodes)
        {
            for (auto& cnode : nodes_level)
            {
                if (!cnode.isDeleted() && cnode.isLeaf())
                {
                    for (auto n : Neighbours)
                    {
                        if (cnode.hasNeighbour(n))
                            file_output << cnode.dumpNeighbourVector(n);
                    }
                }
            }
        }
    }
    file_output.close();
    return;
}
void exportForestNeighbours12(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;
    file_output << "ZONE T = \"Forest neighbours12\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;

    for (auto& tree : forest.trees)
    {
        for (auto& nodes_level : tree.nodes)
        {
            for (auto& cur_node : nodes_level)
            {
                if (!cnode.isDeleted() && cnode.isLeaf())
                {
                    for (auto n12 : Neighbours12)
                    {
                        if (cur_node.hasNeighbour12(n12))
                            file_output << cur_node.dumpNeighbour12Vector(n12);
                    }
                }
            }
        }
    }
    file_output.close();
    return;
}
void exportForestEdges(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\", \"dx\", \"dy\"" << endl;

    //вектора соседей
    file_output << "ZONE T = \"Forest edge neigbours\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;
    for (auto& edge : forest.edges)
    {
        if (!edge.is_deleted)
        {
            file_output << edge.dumpNeigboursVectors();
            //auto& nnode1 = TreeNode::getNode(edge.n1);
            //auto& nnode2 = TreeNode::getNode(edge.n2);
            //if (!nnode1.is_leaf || !nnode2.is_leaf)
                //cout << "suspicious edge! " << edge << endl;
        }
    }

    //точки квадратуры
    file_output << endl << "ZONE T = \"Edge quadrature points\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;
    for (auto& edge : forest.edges)
    {
        if (!edge.is_deleted)
        {
            file_output << edge.dumpQuadraturePoints();
        }
    }
    return;
}
void exportForestEdgeFluxes(std::string filename)
{
    std::ofstream file_output;
    file_output.open(filename);
    file_output << "VARIABLES = \"x\", \"y\"";
    for (auto eq : Equations)
        file_output << ", \"flux" << eq << "\"";
    file_output << endl;
    file_output << "ZONE T = \"Edge Fluxes\" ";
    file_output << "STRANDID = 1 SOLUTIONTIME = " << std::to_string(globals.time) << endl;
    for (auto& edge : forest.edges)
    {
        if (!edge.is_deleted)
        {
            file_output << edge.dumpFlux();
        }
    }
    return;
}*/

//binary writing functions for .plt files
static int writeBinary4BytesInt(__int32 value, FILE* file)
{
    __int32 dword;
    dword = value;
    fwrite(&dword, 4, 1, file);
    return 0;
}
static int writeBinary4BytesFloat(float value, FILE* file)
{
    __int32 dword;
    *((float*)&dword) = value;
    fwrite(&dword, 4, 1, file);
    return 0;
}
static int writeBinary8BytesDouble(double value, FILE* file)
{
    __int64 qword;
    *((double*)&qword) = value;
    fwrite(&qword, 8, 1, file);
    return 0;
}
static int writeBinaryString(std::string str, FILE* file)
{
    __int32 dword;
    int p;
    p = 0;
    do
    {
        dword = str.c_str()[p++];
        fwrite(&dword, 4, 1, file);
    } while (dword);
    return 0;
}
static void writeTecplotFileHeader(FILE* fp)
{
    fwrite("#!TDV111", 8, 1, fp); //Код версии данных Tecplot, нужно именно "1 символ = 1 байт"
    writeBinary4BytesInt(1, fp); //byte order
    writeBinary4BytesInt(0, fp); //filetype: 0 = full, 1 = grid, 2 = solution
    writeBinaryString("Frame title", fp);  //frame title (visible if "frame header" is enabled)
    writeBinary4BytesInt(TECPLOT_FIELDS_NUMBER, fp); //Число массивов, включая координаты
    writeBinaryString("x", fp); //Названия параметров
    writeBinaryString("y", fp);
    writeBinaryString("<greek>r</greek>", fp);
    writeBinaryString("U", fp);
    writeBinaryString("V", fp);
    writeBinaryString("P", fp);
    writeBinaryString("Mach", fp);
    writeBinaryString("level", fp);
    writeBinaryString("magGradRho", fp);
}
static void writeTecplotZoneHeader(FILE* fp)
{
    std::stringstream stream;
    writeBinary4BytesFloat(299.0f, fp); //zone marker
    stream << std::fixed << std::setprecision(10) << globals.time; //zone name
    writeBinaryString(stream.str(), fp);
    writeBinary4BytesInt(-1, fp); //parent zone
    writeBinary4BytesInt(-2, fp); //time strand (-2 = auto, -1 = none, 0+ = explicit value)
    writeBinary8BytesDouble(globals.time, fp); //время в формате double
    writeBinary4BytesInt(-1, fp); //not used, set to -1
    writeBinary4BytesInt(0, fp); //zone type, 0 = ordered
    writeBinary4BytesInt(0, fp); //0 - по блокам, 1 - по точкам
    writeBinary4BytesInt(0, fp); //var location, 0 = all data is at the nodes
    writeBinary4BytesInt(0, fp); //0 = no face neighbor data supplied
    writeBinary4BytesInt(0, fp); //0 = no user defined face connections
    int grid_points_x = (int)(config.Nx - 2) * (int)pow(2, (int)config.max_depth - 1); //число узлов при максимальном измельчении (без ghost-деревьев)
    int grid_points_y = (int)(config.Ny - 2) * (int)pow(2, (int)config.max_depth - 1);
    writeBinary4BytesInt(1, fp); //число узлов вдоль оси Z (Imax), fastest index
    writeBinary4BytesInt(grid_points_y, fp); //slower index (Jmax)
    writeBinary4BytesInt(grid_points_x, fp); //slowest index (Kmax)
    writeBinary4BytesInt(0, fp); //no auxiliary names/values
}
static void writeTecplotFinFileHeader(FILE* fp)
{
    writeBinary4BytesFloat(357.0f, fp); //end of header marker
}
static void writeTecplotZoneDataFloat(FILE* fp)
{
    writeBinary4BytesFloat(299.0f, fp); //zone marker
    //коды типов данных для каждой величины
    for (int field = 0; field < TECPLOT_FIELDS_NUMBER; field++)
        writeBinary4BytesInt(1, fp); //1 = float
    writeBinary4BytesInt(0, fp); //no passive variables
    writeBinary4BytesInt(0, fp); //no variable sharing
    writeBinary4BytesInt(-1, fp); //no connectivity list zone
    //минимумы и максимумы
    dataExtrema extrema = forest.getExtrema();
    for (int field = 0; field < TECPLOT_FIELDS_NUMBER; field++)
    {
        writeBinary8BytesDouble(extrema.minima[field], fp);
        writeBinary8BytesDouble(extrema.maxima[field], fp);
    }
    int Kmax = (int)(config.Nx - 2) * (int)pow(2, config.max_depth - 1) - 1; //-2 ghost-дерева
    int Lmax = (int)(config.Ny - 2) * (int)pow(2, config.max_depth - 1) - 1;
    double h = (extrema.maxima[0] - extrema.minima[0]) / (Kmax + 1);
    //массивы значений
    for (int field = 0; field < TECPLOT_FIELDS_NUMBER; field++)
    {
        for (int k = 0; k <= Kmax; k++)
        {
            for (int l = 0; l <= Lmax; l++)
            {
                Point p{ extrema.minima[0] + k * h + 0.5 * h, extrema.minima[1] + l * h + 0.5 * h };
                auto& ctree = forest.getTreeByCoords(p);
                auto& cnode = ctree.getNodeByCoords(p);
                //auto& dcenter = cnode.dataRef();
                //auto d = cnode.evalPolynomialAt(p);
                auto& d = cnode.dataRef();
                switch (field)
                {
                case 0:
                    writeBinary4BytesFloat((float)NaNcleared(p.x), fp);
                    break;
                case 1:
                    writeBinary4BytesFloat((float)NaNcleared(p.y), fp);
                    break;
                case 2:
                    writeBinary4BytesFloat((float)NaNcleared(d.rho()), fp);
                    break;
                case 3:
                    writeBinary4BytesFloat((float)NaNcleared(d.u()), fp);
                    break;
                case 4:
                    writeBinary4BytesFloat((float)NaNcleared(d.v()), fp);
                    break;
                case 5:
                    writeBinary4BytesFloat((float)NaNcleared(d.p()), fp);
                    break;
                case 6:
                    writeBinary4BytesFloat((float)NaNcleared(d.Mach()), fp);
                    break;
                case 7:
                    writeBinary4BytesFloat((float)cnode.tag().depth(), fp);
                    break;
                case 8:
                    writeBinary4BytesFloat((float)NaNcleared(cnode.magGradRho()), fp);
                    break;
                }
            }
        }
    }
}
static void exportForestPlt(std::string filename)
{
    FILE* fp;
    errno_t err;
    err = fopen_s(&fp, filename.c_str(), "wb");
    if (err == 0)
    {
        //ComputePolynomialCoeffs(0);
        writeTecplotFileHeader(fp);
        writeTecplotZoneHeader(fp);
        writeTecplotFinFileHeader(fp);
        writeTecplotZoneDataFloat(fp);
        fclose(fp);
    }
    else
        cout << "couldn't open file for binary writing: " << filename << endl;
}

void ExportForest() //вывод всего леса (кроме ghost'ов)
{
    std::stringstream filename;

    //if(globals.timestep_number > 0)
        //CalculateNumericalEntropyProduction();

    if (config.exportScatter)
    {
        filename.str("");
        filename << "output_" << std::setw(3) << std::setfill('0') << globals.export_number << ".dat";
        forest.exportScatter(filename.str());
        cout << "Exported to " << filename.str() << "." << endl;
        globals.file_output << "Exported to " << filename.str() << "." << endl;
    }
    if (config.exportNeighbours)
    {
        filename.str("");
        filename << "output_neighbours_" << std::setw(3) << std::setfill('0') << globals.export_number << ".dat";
        forest.exportNeighbours(filename.str()); //соседи
        cout << "Exported to " << filename.str() << "." << endl;
        globals.file_output << "Exported to " << filename.str() << "." << endl;
    }
    if (config.exportNeighbours12)
    {
        filename.str("");
        filename << "output_neighbours12_" << std::setw(3) << std::setfill('0') << globals.export_number << ".dat";
        forest.exportNeighbours12(filename.str()); //соседи12
        cout << "Exported to " << filename.str() << "." << endl;
        globals.file_output << "Exported to " << filename.str() << "." << endl;
    }
    /*if (config.exportEdges)
    {
        filename.str("");
        filename << "output_edges_" << std::setw(3) << std::setfill('0') << globals.export_number << ".dat";
        exportForestEdges(filename.str()); //ребра
        cout << "Exported to " << filename.str() << "." << endl;
        globals.file_output << "Exported to " << filename.str() << "." << endl;
    }
    if (config.exportEdgeFluxes)
    {
        filename.str("");
        filename << "output_edgeFluxes_" << std::setw(3) << std::setfill('0') << globals.export_number << ".dat";
        exportForestEdgeFluxes(filename.str()); //потоки
        cout << "Exported to " << filename.str() << "." << endl;
        globals.file_output << "Exported to " << filename.str() << "." << endl;
    }*/
    if (config.exportPlt)
    {
        filename.str(""); //очистка строки буффера
        filename << "output_" << std::setw(3) << std::setfill('0') << globals.export_number << ".plt";
        exportForestPlt(filename.str());
        cout << "Exported to " << filename.str() << "." << endl;
        globals.file_output << "Exported to " << filename.str() << "." << endl;
    }
    globals.export_number++;
    globals.last_exported_time = globals.time;
}
