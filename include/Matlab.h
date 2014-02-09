/*! \file   Matlab.h
    \brief  Matlab header file.
*/

#ifndef MATLAB_H_
#define MATLAB_H_

#include <fstream>
#include <iostream>
#include <algorithm>

#include "CWS.h"
#include "CVRP.h"

/*! \def    MATLAB_TEXT
    \brief  
*/
#define MATLAB_TEXT (const char *) "text("

/*! \def    MATLAB_LEGEND
    \brief  
*/
#define MATLAB_LEGEND (const char *) "legend("

/*! \def    MATLAB_PLOT
    \brief  
*/
#define MATLAB_PLOT (const char *) "plot("

/*! \def    MATLAB_AXIS
    \brief  
*/
#define MATLAB_AXIS (const char *) "axis("

/*! \def    MATLAB_NUM2STR
    \brief  
*/
#define MATLAB_NUM2STR (const char *) "num2str(["

/*! \def    MATLAB_TERM
    \brief  
*/
#define MATLAB_TERM (const char *) ");"

/*! \def    MATLAB_STATION_DISP_OPT
    \brief  Scatter plot with red bordered blue circles.
*/
#define MATLAB_STATION_DISP_OPT (const char *) "'ro', 'MarkerSize', 20, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [1,1,1]);"

/*! \def    MATLAB_ROUTE_DISP_OPT_SIZE
    \brief  8 options for displaying routes.
*/
#define MATLAB_ROUTE_DISP_OPT_SIZE 8

class Matlab {

    public:

        Matlab();
        ~Matlab();

        enum GraphOptions {SIMPLE, DRAW_ROUTES};

        static void GenerateCVRPGraph(
                CVRP * graph, 
                Matlab::GraphOptions options,
                std::string suffix);
};

#endif  // MATLAB_H_

