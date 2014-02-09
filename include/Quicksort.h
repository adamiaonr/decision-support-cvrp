/*! \file   Quicksort.h
    \brief  Quicksort header file.
*/

#ifndef QUICKSORT_H_
#define QUICKSORT_H_

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
#include <cstring>

#include "CWS.h"

using namespace std;

class Quicksort
{
    public:

        Quicksort();
        ~Quicksort();

        static void Sort(CWS::Saving * input, long p, long r);
};

#endif  // QUICKSORT_H_
