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

/*! \def    GA_DEFAULT_GENERATIONS
    \brief  
*/
#define GA_DEFAULT_GENERATIONS (uint32_t) 10000

/*! \def    GA_DEFAULT_POPULATION_SIZE
    \brief  
*/
#define GA_DEFAULT_POPULATION_SIZE (long) 120

/*! \def    GA_DEFAULT_MUTATION_PROB
    \brief  
*/
#define GA_DEFAULT_MUTATION_PROB (double) 0.5

/*! \def    MAX_DRAWS
    \brief  
*/
#define MAX_DRAWS 5

class Evolutionary {

    public:

        Evolutionary();
        ~Evolutionary();

        static void Run(
                CVRP * cws,
                CVRP * sa,
                uint32_t max_generations,
                long pop_size,
                double mutation_prob,
                bool matlab);

        static void Run(
                std::string dataset,
                uint32_t max_generations,
                long pop_size,
                double mutation_prob,
                bool matlab);
};

#endif  // EVOLUTIONARY_H_
