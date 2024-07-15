#include "json.hpp" //NLohmann json lib

#include "main.h"
#include "CellBox.h"
#include "ProblemConfig.h"

ProblemConfig::ProblemConfig(std::string filename) //конструктор по json-файлу
{
    cout << "Reading config from " << filename << endl;
    nlohmann::json j;
    try
    {
        j = nlohmann::json::parse(std::ifstream(filename), nullptr, true, true); //allow exceptions and c-like comments in json
    }
    catch (nlohmann::json::exception& e)
    {
        cout << "Bad input file: " << e.what() << endl;
    }

    global_box = CellBox({ j["box"]["x min"], j["box"]["y min"] }, { j["box"]["x max"], j["box"]["y max"] });

    problem = j["flow type"];
    shock_position_x = j["flow geometry"]["shock position x"];
    if (problem == "layer")
    {
        layer_right = j["flow geometry"]["layer edges"]["right"];
        layer_bottom = j["flow geometry"]["layer edges"]["bottom"];
        layer_top = j["flow geometry"]["layer edges"]["top"];
    }
    else if (problem == "bubble")
    {
        bubble_axle_x = j["flow geometry"]["bubble axles"]["x"];
        bubble_axle_y = j["flow geometry"]["bubble axles"]["y"];
    }
    gamma = j["gamma"];
    Mach = j["Mach"];
    Atwood = j["Atwood"];

    coord_type = (j["coord type"] == "axisymmetric" ? CoordType::axisymmetric : CoordType::cartesian);
    boundary_conditions[static_cast<int>(Directions::left)] = (j["boundary conditions"]["x min"] == "wall" ? BCType::wall : BCType::ddn_zero);
    boundary_conditions[static_cast<int>(Directions::right)] = (j["boundary conditions"]["x max"] == "wall" ? BCType::wall : BCType::ddn_zero);
    boundary_conditions[static_cast<int>(Directions::down)] = (j["boundary conditions"]["y min"] == "wall" ? BCType::wall : BCType::ddn_zero);
    boundary_conditions[static_cast<int>(Directions::up)] = (j["boundary conditions"]["y max"] == "wall" ? BCType::wall : BCType::ddn_zero);

    rk_order = j["RK order"];
    Courant = j["Courant"];
    if (j["flux"] == "Riemann exact")
        flux = FluxType::Riemann_exact;
    else if (j["flux"] == "HLLC")
        flux = FluxType::HLLC;
    else
        flux = FluxType::LF; //необязательно, FLUX_LF по дефолту
    end_time = j["finish time"];
    exportPeriod = j["export"]["period"];
    exportScatter = j["export"]["scatter"];
    exportNeighbours = j["export"]["neighbours"];
    exportNeighbours12 = j["export"]["neighbours12"];
    exportEdges = j["export"]["edges"];
    exportEdgeFluxes = j["export"]["fluxes"];
    exportPlt = j["export"]["plt"];

    Nx = j["trees"]["Nx"];
    Ny = j["trees"]["Ny"]; //проверять на соответствие с global_box!
    max_depth = j["trees"]["depth"];
    meshRefinePeriod = j["refine period"];
    meshRefinelevel0 = j["refine criterion"]["level 0"];
    meshRefinefactor = j["refine criterion"]["factor"];
    meshInitialRefinePadding = j["initial refine padding"];
    meshRefineAll = j["refine all"]; //измельчить всё до упора и далее не склеивать

    return;
}

std::ostream& operator<<(std::ostream& os, const ProblemConfig& c) //output overload
{
    os << std::boolalpha; //для словесного вывода boolean'ов
    os << "Global box: " << c.global_box.bottom_left() << " to " << c.global_box.top_right() << endl;
    os << "Boundary conditions: top = " << (int)config.boundary_conditions[static_cast<int>(Directions::up)] << ", right = " << (int)config.boundary_conditions[static_cast<int>(Directions::right)] << ", bottom = " << (int)config.boundary_conditions[static_cast<int>(Directions::down)] << ", left = " << (int)config.boundary_conditions[static_cast<int>(Directions::left)] << std::endl;
    os << "Problem: " << c.problem << ", gamma = " << c.gamma << ", Mach = " << c.Mach << ", Atwood = " << c.Atwood << endl;
    if (config.problem == "layer")
        os << "  layer geometry: right = " << c.layer_right << ", bottom = " << c.layer_bottom << ", top = " << c.layer_top << endl;
    if (config.problem == "bubble")
        os << "  bubble geometry: x axle = " << c.bubble_axle_x << ", y axle = " << c.bubble_axle_y << endl;
    os << "Trees: " << c.Nx << " x " << c.Ny << ", max depth = " << c.max_depth << ", ref. period = " << c.meshRefinePeriod << ", by (level 0 = " << c.meshRefinelevel0 << ", factor = " << c.meshRefinefactor << ")";
    os << ", refine all = " << c.meshRefineAll << endl;
    os << "Finish time = " << c.end_time << ", export period = " << c.exportPeriod << endl;
    os << "RK order = " << c.rk_order << ", Courant = " << c.Courant << ", flux method = " << (int)c.flux << endl;
    return os;
}

