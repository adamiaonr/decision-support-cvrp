/*! \file   CWS.h
    \brief  CWS header file.
*/

#ifndef CWS_H_
#define CWS_H_

#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstring>
#include <sys/time.h>

#include "CVRP.h"
#include "Matlab.h"

/*! \def    PARALLEL_VERSION
    \brief  String recognized as an input option for running the 
            program with the parallel version.
*/
#define PARALLEL_VERSION (const char *) "p"

/*! \def    SEQUENTIAL_VERSION
    \brief  String recognized as an input option for running the 
            program with the sequential version.
*/
#define SEQUENTIAL_VERSION (const char *) "s"

using namespace std;

class CWS
{
    public:

        CWS();
        ~CWS();

        enum MergePosition {LEFT, RIGHT};

        enum Version {
                PARALLEL = 0, 
                SEQUENTIAL = 1};

        typedef struct Node { 

            int index;
            int demand;

            CVRP::euc_2d_coords coords;

        } Node;

        typedef struct Saving { 

            // Indexes (i,j) of the saving.
            int i, j;

            // Pointers to Station instances representing i and j.
            CVRP::Station * station_i, * station_j;

            // Distance between stations i and j.
            long dist_ij;

            // The savings associated with the link (i,j).
            long s_ij;

            // This option is used for the sequential version, as the savings 
            // list is cycled numerous times and only considering active 
            // savings.
            bool active;

        } Saving;

        static void Run(
                CVRP * graph,
                CWS::Version version);
};

#endif  // CWS_H_
