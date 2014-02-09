/*! \file   Evolutionary.h
    \brief  Evolutionary header file.
*/

#ifndef EVOLUTIONARY_H_
#define EVOLUTIONARY_H_

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

#include "CWS.h"
#include "CVRP.h"
#include "Matlab.h"

/*! \def    MAX_DRAWS
    \brief  
*/
#define MAX_DRAWS 5

class Evolutionary {

    public:

        Evolutionary();
        ~Evolutionary();

        static void Run(
                CVRP * cws_parallel,
                CVRP * cws_sequential,
                CVRP * sa,
                uint32_t max_generations,
                long pop_size,
                double mutation_prob);
};

#endif  // EVOLUTIONARY_H_
