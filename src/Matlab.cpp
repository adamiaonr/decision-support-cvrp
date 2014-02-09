#include "Matlab.h"

const char * MATLAB_ROUTE_DISP_OPT[] = {
    "'-ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);", 
    "'-go', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);", 
    "'-bo', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);", 
    "'-ko', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);",
    "'--ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);", 
    "'--go', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);", 
    "'--bo', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);", 
    "'--ko', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.5,0.5,0.5]);"
};

//! Matlab constructor.
/*!
*/
Matlab::Matlab()
{
}

//! Matlab destructor.
/*!
*/
Matlab::~Matlab()
{
}

std::string GenerateCVRPFilename(std::string filename, std::string suffix) {

    char aux[16];

    sprintf(aux, "_%s", suffix.c_str());

    filename += aux + std::string(".m");
    filename.replace(filename.find(".vrp"), 4, "");
    std::replace(filename.begin(), filename.end(), '-', '_');

    return filename;
}

void GenerateCVRPFile(std::string filename, std::list <std::string> lines) {

    // Create the .m output file.
    ofstream cvrp_file;
    cvrp_file.open (filename.c_str());

    // Write the content of 'lines' into it.

    std::list <std::string>::iterator it;
    for (it = lines.begin(); it != lines.end(); it++) {

        cvrp_file << (*it);
    }

    cvrp_file.close();
}

void AddStation2Lines(
        CVRP::Station * station,
        std::string & station_x, 
        std::string & station_y,
        std::string & station_index,
        bool terminate) {

    // 10 chars, more than enough for the coordinates and indexes we're dealing 
    // with here. General purpose buffer for building up intermediate strings
    // for MATLAB commands.
    char station_buffer[10];

    // Clear station_buffer just to be sure.
    memset(station_buffer, 0, sizeof(station_buffer));

    // Add the x coordinates, separated by a ';' if another 
    // station follows it. If not, add the terminating sequence for this 
    // MATLAB set of x coordss, "], ".
    sprintf(
            station_buffer, 
            "%3d%s", 
            station->coords.x, 
            (!terminate ? ";" : "], "));

    station_x += station_buffer;

    memset(station_buffer, 0, sizeof(station_buffer));

    sprintf(
            station_buffer, 
            "%3d%s", 
            station->coords.y, 
            (!terminate ? ";" : "], "));

    station_y += station_buffer;

    memset(station_buffer, 0, sizeof(station_buffer));

    sprintf(
            station_buffer, 
            "%4d%s", 
            station->index - 1, 
            (!terminate ? ";" : "])"));

    station_index += station_buffer;
}

void Matlab::GenerateCVRPGraph(
        CVRP * graph, 
        Matlab::GraphOptions options, 
        std::string suffix) {

    // List of strings which hold the MATLAB code generated within this 
    // function.
    std::list <std::string> lines;

    // Just some hardcoded lines here...
    lines.push_back(std::string("hold on;\n"));

    CVRP::Station * depot = graph->GetDepot();

    // Simple legend for the graph.
    std::string legend = std::string(MATLAB_LEGEND);
    char aux_text[16];

    // If the DRAW_ROUTES option was selected, create several 'plot' commands, 
    // one for each route, with different line colours.
    if (options == Matlab::DRAW_ROUTES) {

        int i = 0, k = 0;

        CVRP::Station * curr_station;

        // For each route (each with a different colour), create 'plot' and 
        // 'text' commands separately for each route.
        std::list <CVRP::Route *>::iterator rit;

        for (rit = graph->routes.begin(); rit != graph->routes.end(); rit++) {

            // Remember, new 'plot' and 'text' commands.
            std::string station_text = std::string(MATLAB_TEXT);
            std::string station_index = std::string(MATLAB_NUM2STR);
            std::string station_x = std::string("[");
            std::string station_y = std::string("[");
            std::string station_circles = std::string(MATLAB_PLOT);

            // Add a depot at the start of each route.
            AddStation2Lines(depot, station_x, station_y, station_index, false);

            curr_station = (*rit)->start_station;

            do {

                AddStation2Lines(
                        curr_station, 
                        station_x, 
                        station_y, 
                        station_index, 
                        false);

                curr_station = curr_station->next;

            } while (curr_station->type != CVRP::DEPOT);

            // Add a depot at the end of each route.
            AddStation2Lines(depot, station_x, station_y, station_index, true);

            // Create the final lines and add them to the line list.
            station_text += 
                    station_x + station_y + station_index + MATLAB_TERM + "\n";

            station_circles +=
                    station_x + station_y 
                    + MATLAB_ROUTE_DISP_OPT[i++] 
                    + "\n";

            if (i >= MATLAB_ROUTE_DISP_OPT_SIZE)
                i = 0;

            lines.push_back(station_circles);
            lines.push_back(station_text);

            // Build up a small legend to the chart.
            memset(aux_text, 0, sizeof(aux_text));

            sprintf(
                    aux_text, 
                    "'Route %d',", 
                    ++k);

            legend += aux_text;
        }

        // Finalize the legend.
        legend += "0);\n";
        lines.push_back(legend);

    } else {

        // Otherwise, create a single 'plot' command, without any routes (i.e. 
        // lines connecting the circles representing the stations).
        std::string station_text = std::string(MATLAB_TEXT);
        std::string station_index = std::string(MATLAB_NUM2STR);
        std::string station_x = std::string("[");
        std::string station_y = std::string("[");
        std::string station_circles = std::string(MATLAB_PLOT);

        std::vector <CVRP::Station *>::iterator it;

        for (it = graph->stations.begin(); it != graph->stations.end(); it++) {

            AddStation2Lines(
                        (*it), 
                        station_x, 
                        station_y, 
                        station_index, 
                        false);
        }

        // Add a dummy depot, just for finalizing the 'plot' and 'text' lines.
        // FIXME: Horrible solution, fix this.
        AddStation2Lines(
                depot, 
                station_x, 
                station_y, 
                station_index, 
                true);

        // Finalize the station_text line to plot the station's indexes (text).
        station_text += 
                station_x + station_y + station_index + MATLAB_TERM + "\n";

        // Add the line to the list.
        lines.push_back(station_text);

        station_circles +=
                station_x + station_y + MATLAB_STATION_DISP_OPT + "\n";

        lines.push_back(station_circles);
    }

    // Create the .m (MATLAB) file, with a name in the format 
    // [DATASET].[YYMMDD].[HHMMSS].m.
    std::string filename = GenerateCVRPFilename(graph->dataset, suffix);

    // Just some hardcoded lines here...
    lines.push_back(std::string("xlabel('x');\n"));
    lines.push_back(std::string("ylabel('y');\n"));

    GenerateCVRPFile(filename, lines);
}

