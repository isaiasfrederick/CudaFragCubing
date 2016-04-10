#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <sys\timeb.h>
#include <thrust\set_operations.h>
#include "gpu.h"
#include<cuda.h>

using namespace std;

typedef struct
{
	struct timeb start;
	struct timeb end;
	int elapsed_time;
} ContadoresTempo;

void reset_contadores_tempo(ContadoresTempo* contadores)
{
	contadores->elapsed_time = 0;
}

void marcar_tempo_inicial(ContadoresTempo* contadores)
{
	ftime(&contadores->start);
}

void marcar_tempo_final(ContadoresTempo* contadores)
{
	ftime(&contadores->end);
}

void contabilizar_tempo_gasto(ContadoresTempo* contadores)
{
	contadores->elapsed_time += (int)(1000.0 * (contadores->end.time - contadores->start.time) + (contadores->end.millitm - contadores->start.millitm));
}

vector<string> split(string str, char delimiter) {
	vector<string> internal;
	stringstream ss(str);
	string tok;

	while (getline(ss, tok, delimiter))
		internal.push_back(tok);

	return internal;
}

vector<int> set_intersection(vector<int>* va, vector<int>* vb)
{
	int ai = 0, bi = 0, ci = 0;
	vector<int> vc;

	int tam_a = va->size();
	int tam_b = vb->size();

	while (ai < tam_a && bi < tam_b)
	{
		if ((*va)[ai] == (*vb)[bi])
		{
			if (vc.empty( ) || (*va)[ai] != vc[ci - 1])
			{
				if ((*va)[ai] != 0)
				{
					vc.push_back((*va)[ai]);
					ci++;
				}
			}

			ai++; bi++;
		}
		else if ((*va)[ai] > (*vb)[bi])
		{
			bi++;
		}
		else if ((*va)[ai] < (*vb)[bi]) {
			ai++;
		}

	}

	return vc;

}

#endif