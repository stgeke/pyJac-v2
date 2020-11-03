/*

A cog-templated skeleton for reading of initial conditions from a binary file

(C) Nicholas Curtis - 2018

Global declarations for Cog:
    - readgen: path to a serialized ReadgenRecord instance
    that may be loaded to generate this file
*/

#include "mechanism.hpp"
#include "vectorization.hpp"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

// size of a single input buffer
// total buffer size
#define BUFF_SIZE ((NN + 1))

//for sanity, the input data is expected to be in C-order

void read_initial_conditions(const char* filename, unsigned int NUM,
                             double* P_arr, double* phi,
                             const char order) {
    FILE *fp = fopen (filename, "rb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file: %s\n", filename);
        exit(-1);
    }

    double buffer[BUFF_SIZE];
    // load temperature, pressure and concentrations for all (cells)
    for (int i = 0; i < NUM; ++i)
    {
        // read line from data file
        int count = fread(buffer, sizeof(double), BUFF_SIZE, fp);
        if (count != (BUFF_SIZE))
        {
            fprintf(stderr, "File (%s) is incorrectly formatted, %d "
                "doubles were expected but only %d were read.\n",
                filename, BUFF_SIZE, count);
            exit(-1);
        }

        //fill the parameter array
        P_arr[i] = buffer[1];

        // phi fill depends on order
        if (order == 'C')
        {
            //fill in temperature
            phi[i * NN] = buffer[0];
            //fill in species moles
            for (int j = 0; j < NS; j++)
            {
                phi[i * NN + (j + 1)] = buffer[j + 2];
            }
        }
        else
        {
            //fill in temperature
            phi[i] = buffer[0];
            //fill in species moles
            for (int j = 0; j < NS; j++)
            {
                phi[(j + 1) * NUM + i] = buffer[j + 2];
            }
        }

    }
}
