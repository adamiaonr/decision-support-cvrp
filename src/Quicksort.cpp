#include "Quicksort.h"

//! Quicksort constructor.
/*!
*/
Quicksort::Quicksort()
{
}

//! Quicksort destructor.
/*!
*/
Quicksort::~Quicksort()
{
}

int Partition(CWS::Saving * input, long p, long r) {

    CWS::Saving pivot = input[r];

    while (p < r) {

        while (input[p].s_ij > pivot.s_ij)
            p++;

        while (input[r].s_ij < pivot.s_ij)
            r--;

        if (input[p].s_ij == input[r].s_ij)
            p++;

        else if (p < r) {

            CWS::Saving tmp = input[p];
            input[p] = input[r];
            input[r] = tmp;
        }
    }

    return r;
}

void Quicksort::Sort(CWS::Saving * input, long p, long r) {

    if (p < r) {

        long j = Partition(input, p, r);
        Sort(input, p, j - 1);
        Sort(input, j + 1, r);
    }
}

