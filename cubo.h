#ifndef CUBO_H
#define CUBO_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <sstream>
#include "cuboid.h"
#include "intersection.h"

using namespace std;

typedef struct
{

	map<int, map<int, vector<int>>> values;
	map<int, vector<int>> map_colunas;
	set<int> cardinalidades_ordenadas;
	vector <int> colunas_ordenadas;
	vector<int> cardinalidades_por_coluna;
	int total_tuplas;

} Cubo;

vector<int> get_all_tuple_ids(int total_tuplas)
{
	vector<int> all_tuples_ids = vector<int>();
	for (register int i = 0; i < total_tuplas; i++)
		all_tuples_ids.push_back(i + 1);

	return all_tuples_ids;
}

inline void adicionar_dimension_value(Cubo* cubo, int value, int tuple_id, int dimension_id)
{
	if (cubo->values.find(dimension_id) == cubo->values.end())
	{
		map<int, vector<int>> map_indices_invertido;
		cubo->values[dimension_id] = map_indices_invertido;
	}

	if (cubo->values[dimension_id].find(value) == cubo->values[dimension_id].end())
	{
		vector<int> lista;
		cubo->values[dimension_id][value] = lista;
	}

	cubo->values[dimension_id][value].push_back(tuple_id);

}

inline void deletar_cubo_base(Cubo** cubo_base)
{
	if (*cubo_base)
	{
		delete *cubo_base;
		*cubo_base = NULL;
	}
}

// leitor da base em formato de texto
Cubo* carregar_cubo_base(string path)
{
/*	Cubo* cubo;

	cubo = new Cubo;

	vector<string> internal;
	ifstream file(path);
	string str;

	int tuple_id = 0;
	int value;

	// Skip na primeira linha, porque essa linha é dispensável
	getline(file, str);

	while (getline(file, str))
	{
		tuple_id++;

		internal = split(str, SEPARADOR_COLUNAS_ARQUIVO);

		for (int dimension_id = 0; dimension_id < internal.size(); dimension_id++)
		{
			value = stoi(internal[dimension_id]);
			adicionar_dimension_value(cubo, value, tuple_id, dimension_id + 1);
		}

	}

	cubo->total_tuplas = tuple_id;

	return cubo;*/
	return NULL;
}

inline Cubo* carregar_cubo_base(char* path, int n_dimensions)
{
	FILE* arquivo = fopen(path, "rb");
	int n_tuplas;
	int value;
	Cubo* cubo;
	int* values = (int*) malloc(sizeof(int) * n_dimensions);

	cubo = new Cubo;

	// lendo a qtde de tuplas e jogando "fora"
	fread(&n_tuplas, sizeof(int), 1, arquivo);

	for (int dimension_id = 0; dimension_id < n_dimensions; dimension_id++)
	{
		fread(&value, sizeof(int), 1, arquivo);
		cubo->cardinalidades_ordenadas.insert(value);
		cubo->map_colunas[value].push_back(dimension_id + 1);
		cubo->colunas_ordenadas.push_back(value);
		cubo->cardinalidades_por_coluna.push_back(value);
	}

	set<int>::iterator it = cubo->cardinalidades_ordenadas.begin();

	while (it != cubo->cardinalidades_ordenadas.end())
	{
		while (!cubo->map_colunas[*it].empty())
		{
			cubo->colunas_ordenadas.push_back(cubo->map_colunas[*it][0]);
			cubo->map_colunas[*it].erase(cubo->map_colunas[*it].begin());
		}

		it++;
	}

	int tuple_id = 1;

	while (tuple_id <= n_tuplas)
	{
		fread(values, sizeof(int), n_dimensions, arquivo);

		for (int dimension_id = 0; dimension_id < n_dimensions; dimension_id++)
			adicionar_dimension_value(cubo, values[dimension_id], tuple_id, dimension_id + 1);

		tuple_id++;

	}

	cubo->total_tuplas = tuple_id - 1;

	free(values);
	fclose(arquivo);

	return cubo;

}

inline void exibir_cubo(Cubo* cubo)
{
	cout << "\n\n\nREPRESENTACAO CUBO EM INDICES INVERTIDOS\n";

	for (map<int, map<int, vector<int>>>::iterator iterator = cubo->values.begin( ); iterator != cubo->values.end(); iterator++)
	{
		cout << "\nDIMENSAO " << iterator->first;

		map<int, vector<int>> mapa = cubo->values[iterator->first];

		for (map<int, vector<int>>::iterator iterator_interno = mapa.begin(); iterator_interno != mapa.end(); iterator_interno++)
		{
			cout << "\n    Value: " << iterator_interno->first << " -> ";
			vector<int> lista = mapa[iterator_interno->first];

			for (vector<int>::iterator iterador_lista = lista.begin(); iterador_lista != lista.end(); iterador_lista++)
				cout << *iterador_lista << ", ";
		}

	}

	cout << "\n\n";
}

inline void exibir_fragmentos(Fragmentos* fragmentos)
{
	map<CubDim, Cuboid, CompCubDim>::iterator it = fragmentos->cuboids.begin();

	cout << "\n\n\nEXIBINDO " << fragmentos->cuboids.size() << " FRAGMENTOS NO TOTAL:\n";

	while (it != fragmentos->cuboids.end())
	{
		cout << "\n\nCUBOID " + label_cuboid_dimensions(it->first);

		Cuboid cuboid = fragmentos->cuboids[it->first];

		map<AttVal, vector<int>, CompAttVal>::iterator it_attribute_value = fragmentos->cuboids[it->first].values.begin();

		cout << "\nEste cuboid possui " << fragmentos->cuboids[it->first].values.size() << " AttVals!";
		cout << "\n";

		while (it_attribute_value != fragmentos->cuboids[it->first].values.end())
		{
			vector<int>* tuple_id_list = &(fragmentos->cuboids[it->first].values[it_attribute_value->first]);
			cout << "\n     " + label_attribute_values(it_attribute_value->first) + " => " + label_tuple_id_list(tuple_id_list);

			it_attribute_value++;
		}

		it++;

	}

	cout << "\n_______________________________________________________";

}

inline Cuboid calcular_cuboid(Cubo* cubo, vector<int> dimensoes)
{
	Cuboid cuboid;
	map<AttVal, vector<int>, CompAttVal>* intersecao_final;

	intersecao_final = new map < AttVal, vector<int>, CompAttVal > ;

	for (int d = 0; d < dimensoes.size(); d++)
	{
		if (d == 0)
		{
			map<int, vector<int>>::iterator it_attribute_values;

			it_attribute_values = cubo->values[dimensoes[d]].begin();

			// Adicionando o ALL para cada dimensão
			AttVal attribute_value_all = create_attribute_value({ ALL });
			(*intersecao_final).insert({attribute_value_all, get_all_tuple_ids(cubo->total_tuplas)});

			while (it_attribute_values != cubo->values[dimensoes[d]].end())
			{	
				AttVal attribute = create_attribute_value({ it_attribute_values->first });

				(*intersecao_final).insert({ attribute, cubo->values[dimensoes[d]][it_attribute_values->first] });
				(*intersecao_final).insert({ attribute, cubo->values[dimensoes[d]][it_attribute_values->first] });

				it_attribute_values++;
			}

		}
		// Com fragmento = 1 o código nunca entra aqui
		else
		{
			map<AttVal, vector<int>, CompAttVal>::iterator it_attribute_values = intersecao_final->begin();
			map<AttVal, vector<int>, CompAttVal>* intersecao_tmp;

			intersecao_tmp = new map < AttVal, vector<int>, CompAttVal >;

			while (it_attribute_values != intersecao_final->end())
			{
				map<int, vector<int>>::iterator it_value = cubo->values[dimensoes[d]].begin();

				while (it_value != cubo->values[dimensoes[d]].end())
				{
					vector<int>* tuple_ids_1 = &((*intersecao_final)[it_attribute_values->first]);
					vector<int>* tuple_ids_2 = &(cubo->values[dimensoes[d]][it_value->first]);

					vector<int> resultado = set_intersection(tuple_ids_1, tuple_ids_2);

					if (!resultado.empty())
					{
						AttVal attribute =
							concatenate_attribute_value(it_attribute_values->first, it_value->first);
						(*intersecao_tmp).insert({ attribute, resultado });

					}

					it_value++;

				}

				it_attribute_values++;

			}

			if (!intersecao_tmp->empty( ))
			{
				map<AttVal, vector<int>, CompAttVal>* intersecao_deletavel;

				intersecao_deletavel = intersecao_final;
				intersecao_final = intersecao_tmp;

				intersecao_tmp = NULL;
				delete intersecao_deletavel;

			}
			else
			{
				delete intersecao_tmp;
			}

			intersecao_tmp = NULL;
			
		}
	}

	cuboid.values = *intersecao_final;

	delete intersecao_final;
	intersecao_final = NULL;

	return cuboid;

}

inline void deletar_fragmentos(Fragmentos** fragmentos)
{

}

vector<vector<int>> gerar_permutacoes(vector<int> dimensoes)
{
	vector<vector<int>> resultado;
	map < int, vector<vector<int>>> por_dim;
	map<int, vector<int>> offsets;

	por_dim.insert({ 1, vector<vector<int>>() });
	offsets.insert({ 1, vector<int>() });

	for (register int i = 0; i < dimensoes.size(); i++)
	{
		vector<int> v = vector<int>({ dimensoes[i] });

		por_dim[1].push_back(v);
		offsets[1].push_back(i + 1);
	}

	for (register int i = 1; i < dimensoes.size(); i++)
	{
		por_dim.insert({ i + 1, vector<vector<int>>() });
		offsets.insert({ i + 1, vector<int>() });
	}

	// tamanho cuboid
	for (register int i = 1; i < por_dim.size(); i++)
	{
		vector<vector<int>> vectors = por_dim[i];
		vector<int> offsets_vector = offsets[i];

		// elemento inicio
		for (register int j = 0; j < vectors.size(); j++)
		{
			vector<int> v = vectors[j];
			int sup = offsets_vector[j];

			for (register int k = sup; k < dimensoes.size(); k++)
			{
				vector<int> v_tmp = v;
				v_tmp.push_back(dimensoes[k]);
				por_dim[i + 1].push_back(v_tmp);
				offsets[i + 1].push_back(k + 1);
			}

		}
	}

	register map < int, vector<vector<int>>>::iterator it = por_dim.begin();
	int k = 0;

	while (it != por_dim.end())
	{
		vector<vector<int>> listas = por_dim[it->first];

		for (register int j = 0; j < listas.size(); j++, k++)
			resultado.push_back(listas[j]);

		it++;
	}

	return resultado;
}

inline map<AttVal, vector<int>, CompAttVal> obterAtributos(Fragmentos* fragmentos, CubDim dimensao)
{
	return fragmentos->cuboids[dimensao].values;
}

inline Fragmentos* fragmentar_cubo(Cubo* cubo)
{
	Fragmentos* fragmentos = new Fragmentos;
	int fragmento_real_computado = 0;
	int tamanho_fragmento = 1;

	map<int, map<int, vector<int>>>::iterator iterador_dimensoes = cubo->values.begin();

	fragmentos->map_colunas = cubo->map_colunas;
	fragmentos->cardinalidades_ordenadas = cubo->cardinalidades_ordenadas;
	fragmentos->colunas_ordenadas = cubo->colunas_ordenadas;
	fragmentos->cardinalidades_por_coluna = cubo->cardinalidades_por_coluna;

	while (iterador_dimensoes != cubo->values.end())
	{
		vector<int> dimensoes;

		for (int d = 0; d < tamanho_fragmento; d++)
		{
			dimensoes.push_back(iterador_dimensoes->first);
			iterador_dimensoes++;

			if (iterador_dimensoes == cubo->values.end())
			{

				if (!fragmento_real_computado)
				{
					fragmento_real_computado = 1;
					fragmentos->tamanho_fragmento = d + 1;
				}

				d = tamanho_fragmento;

			}

		}

		Cuboid cuboid_retorno = calcular_cuboid(cubo, dimensoes);

		if (!cuboid_retorno.values.empty())
		{
			CubDim cuboid_dimensions = create_cuboid_dimensions(dimensoes);
			fragmentos->cuboids.insert({ cuboid_dimensions, cuboid_retorno });
		}

	}

	fragmentos->total_tuplas = cubo->total_tuplas;
	fragmentos->all_ids = get_all_tuple_ids(cubo->total_tuplas);

	return fragmentos;

}

#endif