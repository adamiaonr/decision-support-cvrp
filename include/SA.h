/*! \file   SA.h
    \brief  SA header file.
*/

#ifndef SA_H_
#define SA_H_

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
#include "CWS.h"
#include "Matlab.h"

/*! \def    SA_DEFAULT_COOLING_RATE
    \brief  
*/
#define SA_DEFAULT_COOLING_RATE (double) 0.005

/*! \def    SA_DEFAULT_TERMINATION_CRITERION
    \brief  
*/
#define SA_DEFAULT_TERMINATION_CRITERION (uint32_t) 25000

/*! \def    SA_DEFAULT_HALTING_CRITERION
    \brief  
*/
#define SA_DEFAULT_HALTING_CRITERION (uint16_t) 300

/*! \def    OSMAN_ICHANGE_MAX_LAMBDA
    \brief  
*/
#define OSMAN_ICHANGE_MAX_LAMBDA 1

/*! \def    UINT32_MAX
    \brief  
*/
// #define UINT32_MAX (0xFFFFFFFF)

/*! \def    MAX_SA_POINT
    \brief  
*/
#define MAX_SA_POINT 10

class SA {

    public:

        SA();
        ~SA();

        enum Phase {PHASE_1, PHASE_2};

        typedef struct CoolingScheduleParams {

            double T_s;
            double T_f;

            uint32_t alpha;
            uint32_t gamma;

            int R;

        } CoolingScheduleParams;

        typedef struct SAPoint { 

            uint32_t k;
            uint32_t i;

            double T_k;

            long S_c;
            long S_n;
            long S_n_;
            long S_b;

        } SAPoint;

        static void Run(
                CVRP * graph, 
                double cooling_rate,
                double T_i,
                uint32_t L,
                uint16_t stop_n,
                SA::Phase start_phase,
                bool matlab);
};

#endif  // SA_H_
