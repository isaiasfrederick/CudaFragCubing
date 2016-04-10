#ifndef COPIA_PARCIAL_H
#define COPIA_PARCIAL_H

#include <iostream>
#include <map>
#include <vector>
#include "cuboid.h"

using namespace std;

typedef struct
{
	map<CubDim, map<AttVal, int*, CompAttVal>, CompCubDim> offsets;
	int total_tids;
} CopiaParcial;

int* obter_offset(CopiaParcial* copia, CubDim cub_dim, AttVal attVal)
{
	map<AttVal, int*, CompAttVal> value = copia->offsets[cub_dim];
	return value[attVal];
}

void inserir_offset(CopiaParcial* copia, CubDim cub_dim, AttVal attVal, int* offset)
{
	copia->offsets[cub_dim][attVal] = offset;
}

#endif