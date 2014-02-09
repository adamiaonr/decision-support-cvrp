/*! \file   DatasetParser.h
    \brief  DatasetParser header file.
*/

#ifndef DATASETPARSER_H_
#define DATASETPARSER_H_

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

/*! \def    MAX_CHARS_PER_LINE
    \brief  
*/
#define MAX_CHARS_PER_LINE 512

/*! \def    MAX_CHARS_PER_LINE
    \brief  
*/
#define MAX_TOKENS_PER_LINE 20

/*! \def    DELIMITER
    \brief  
*/
#define DELIMITER (const char *) " "

using namespace std;

class DatasetParser
{
    public:

        DatasetParser();
        ~DatasetParser();

        typedef struct Node { 

            int index;
            int demand;

            CVRP::euc_2d_coords coords;

        } Node;

        static int ParseDataset(
                const char * dataset,
                CVRP * graph);
};

#endif  // DATASETPARSER_H_
