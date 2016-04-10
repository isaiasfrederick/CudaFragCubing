#ifndef CUBOID_H
#define CUBOID_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "constantes.h"

using namespace std;

typedef struct
{
	vector<int> dimensions;

} CubDim;

typedef struct
{
	vector<int> values;

} AttVal;

string label_attribute_values(AttVal attribute);

typedef struct
{
	bool operator() (CubDim cd1, CubDim cd2) const {
		for (register int  i = 0; i < cd1.dimensions.size() && i < cd2.dimensions.size(); i++)
			if (cd1.dimensions[i] > cd2.dimensions[i])
				return true;
			else if (cd1.dimensions[i] < cd2.dimensions[i])
				return false;

		return cd1.dimensions.size() > cd2.dimensions.size();

	}
} CompCubDim;

typedef struct
{
	bool operator() (AttVal av1, AttVal av2) const {
		for (register int i = 0; i < av1.values.size() && i < av2.values.size(); i++)
			if (av1.values[i] > av2.values[i])
				return true;
			else if (av1.values[i] < av2.values[i])
				return false;

		return av1.values.size() > av2.values.size();

	}
} CompAttVal;

typedef struct
{
	map<AttVal, vector<int>, CompAttVal> values;

} Cuboid;

typedef struct
{
	int tamanho_fragmento;
	int total_tuplas;
	vector<int> all_ids;
	map<CubDim, Cuboid, CompCubDim> cuboids;

	set<int> cardinalidades_ordenadas;
	vector<int> cardinalidades_por_coluna;
	map<int, vector<int>> map_colunas;
	vector<int> colunas_ordenadas;

} Fragmentos;

inline Cuboid* get_cuboid(Fragmentos* fragmentos, CubDim cub_dim)
{
	return &(fragmentos->cuboids[cub_dim]);
}

inline map < AttVal, vector<int>>::iterator get_attval_iterator(Cuboid* cuboid)
{
	return cuboid->values.begin();
}

inline vector<int> get_tids(Cuboid* cuboid, AttVal attVal)
{
	return cuboid->values[attVal];
}

inline bool attribute_values_equals(AttVal av1, AttVal av2)
{
	set<AttVal, CompAttVal> set;
	set.insert(av1);

	return set.find(av2) != set.end();
}

inline CubDim create_cuboid_dimensions(vector<int> dimensions)
{
	CubDim cd;

	sort(dimensions.begin(), dimensions.end());
	cd.dimensions = dimensions;

	return cd;
}

inline CubDim create_cuboid_dimensions(CubDim cub_dim, int dimensao)
{
	CubDim cd = cub_dim;

	cd.dimensions.push_back(dimensao);
	sort(cd.dimensions.begin( ), cd.dimensions.end( ));

	return cd;
}

inline CubDim create_cuboid_dimensions(int dimension)
{
	CubDim cd;

	vector<int> dimensions = { dimension };

	sort(dimensions.begin(), dimensions.end());
	cd.dimensions = dimensions;

	return cd;
}

inline CubDim create_cuboid_dimensions(CubDim cd1, CubDim cd2)
{
	CubDim cd = cd1;

	cd.dimensions.insert(cd.dimensions.end( ), cd2.dimensions.begin( ), cd2.dimensions.end( ));
	sort(cd.dimensions.begin( ), cd.dimensions.end( ));

	return cd;
}

inline AttVal create_attribute_value(vector<int> values)
{
	AttVal av;

	sort(values.begin(), values.end());
	av.values = values;

	return av;
}

inline AttVal create_attribute_value(int value)
{
	AttVal av;
	vector<int>values = { value };

	sort(values.begin( ), values.end( ));
	av.values = values;

	return av;
}

inline AttVal concatenate_attribute_value(AttVal av1, AttVal av2)
{
	AttVal av_retorno = av1;

	av_retorno.values.insert(av_retorno.values.end( ), av2.values.begin( ), av2.values.end( ));
	sort(av_retorno.values.begin( ), av_retorno.values.end( ));

	return av_retorno;
}

inline AttVal concatenate_attribute_value_ws(AttVal av1, AttVal av2)
{
	AttVal av_retorno = av1;

	av_retorno.values.insert(av_retorno.values.end( ), av2.values.begin( ), av2.values.end( ));

	return av_retorno;
}

inline AttVal concatenate_attribute_value(AttVal av1, AttVal av2, CubDim cd1, CubDim cd2)
{
	AttVal av_retorno;
	int pos_cd1 = 0;
	int pos_cd2 = 0;

	while (pos_cd1 < cd1.dimensions.size() && pos_cd2 < cd2.dimensions.size())
	{
		if (cd1.dimensions[pos_cd1] < cd2.dimensions[pos_cd2])
		{
			av_retorno.values.push_back(av1.values[pos_cd1]);
			pos_cd1++;
		}
		else
		{
			av_retorno.values.push_back(av2.values[pos_cd2]);
			pos_cd2++;
		}
	}
	
	while (pos_cd1 < cd1.dimensions.size())
	{
		av_retorno.values.push_back(av1.values[pos_cd1]);
		pos_cd1++;
	}

	while (pos_cd2 < cd2.dimensions.size())
	{
		av_retorno.values.push_back(av2.values[pos_cd2]);
		pos_cd2++;
	}

	return av_retorno;
}

inline AttVal concatenate_attribute_value(AttVal av1, int v)
{
	AttVal av_retorno = av1;

	av_retorno.values.push_back(v);
	sort(av_retorno.values.begin(), av_retorno.values.end());

	return av_retorno;
}

inline CubDim concatenate_cuboid_dimensions(CubDim cd1, CubDim cd2)
{
	CubDim cd_retorno = cd1;

	cd_retorno.dimensions.insert(cd_retorno.dimensions.end( ), cd2.dimensions.begin( ), cd2.dimensions.end( ));
	sort(cd_retorno.dimensions.begin( ), cd_retorno.dimensions.end( ));

	return cd_retorno;
}

inline CubDim concatenate_cuboid_dimensions(CubDim cd1, int d)
{
	CubDim cd_retorno = cd1;

	cd_retorno.dimensions.push_back(d);
	sort(cd_retorno.dimensions.begin( ), cd_retorno.dimensions.end( ));

	return cd_retorno;
}

template <typename T>
inline string label_vector(vector<T> v)
{
	string label = "";

	register int d;
	for (d = 0; d < v.size( ); d++)
		label += to_string(v[d]) + " ";

	return label;
}

inline string label_cuboid_dimensions(CubDim cuboid)
{
	string label = "";

	register int d;
	for (d = 0; d < cuboid.dimensions.size( ) - 1; d++)
		label += to_string(cuboid.dimensions[d]) + " ";

	label += to_string(cuboid.dimensions[d]);
	return label;
}

inline string label_attribute_values(AttVal attribute)
{
	string label = "";

	for (int v = 0; v < attribute.values.size()-1; v++)
		if (attribute.values[v] == ALL)
			label += string({AGGREGATED_CARACTER}) + " ";
		else
			label += to_string(attribute.values[v]) + " ";

	int index = attribute.values.size() - 1;
	if (attribute.values[index] == ALL)
		label += string({ AGGREGATED_CARACTER });
	else
		label += to_string(attribute.values[index]);

	return label;
}

inline string label_tuple_id_list(vector<int>* tuple_id_list)
{
	string label = "";

	int i = 0;
	for (i = 0; i < tuple_id_list->size()-1; i++)
		label += to_string((*tuple_id_list)[i]) + ", ";

	label += to_string((*tuple_id_list)[i]);
	return label;
}

#endif