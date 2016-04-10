#ifndef QUERY_H
#define QUERY_H

#include <vector>
#include <map>
#include <iostream>
#include <list>
#include <map>
#include <fstream>
#include <string>
#include <thrust\copy.h>
#include <Windows.h>
#include <sstream>
#include <thread>
#include <chrono>
#include "cuboid.h"
#include "intersection.h"
#include "copiaparcial.h"
#include "constantes.h"
#include "cubo.h"
#include "gpu.h"
#include "thrust\set_operations.h"

using namespace std;

#define DEBUG_QUERY_H 0

typedef map<AttVal, vector<int>, CompAttVal> AttValsIndexados;
typedef map<CubDim, AttValsIndexados, CompCubDim> ResultadoInquired;
typedef map<CubDim, AttValsIndexados, CompCubDim> AttValsIndexadosPorCuboid;

typedef struct
{

	string comando_str;
	vector<string> comando;

	vector<int> dimensoes_instanciadas;
	vector<int> dimensoes_inquired;
	// Ordenadas por cardinalidade
	set<int> set_dimensoes_inquired;
	vector<int> dimensoes_aggregated;
	vector<int> valores_instanciados;
	vector<int> mapeamento_colunas;

} Query;

inline void indexar_tuplas(Query query, map<string, long>* detector_colisoes, AttValsIndexados resultado);

inline string construir_attribute_label(Query query, AttVal attribute_value);

inline void processar_consulta(gpu::GpuVars* gpu_vars, Fragmentos* fragmentos, Query query, map<string, long>* detector_colisoes);

inline AttValsIndexados* processar_dimensoes_instanciadas(gpu::GpuVars* gpu_vars, Fragmentos* fragmentos, Query query);

inline AttValsIndexados* processar_dimensoes_instanciadas(Fragmentos* fragmentos, Query query);

inline ResultadoInquired* processar_dimensoes_inquired(Fragmentos* fragmentos, Query query, AttValsIndexados* dimensoes_instanciadas);

inline ResultadoInquired* processar_dimensoes_inquired(gpu::GpuVars* gpu_vars, Fragmentos* fragmentos, Query query, AttValsIndexados* dimensoes_instanciadas);

inline ResultadoInquired* merge_cuboid_inquired(ResultadoInquired* dimensoes_inquired);

inline AttValsIndexados* processar_dimensoes_agregadas(Fragmentos* fragmentos, Query query);

inline AttValsIndexadosPorCuboid *coletar_todas_dimensoes_inquired(Fragmentos* fragmentos, Query query);

inline void mensagem_executando_query(string query)
{
	std::cout << "\n    Executing query: " << query;
}

inline void mensagem_resultado_query(int total_rows, int ms)
{
	std::cout << "\n\n    " << total_rows << " rows returned. Query executed in " << ms << " ms.";
}

inline void processar_consulta(gpu::GpuVars* gpu_vars, Fragmentos* fragmentos, Query query, map<string, long>* detector_colisoes)
{
	ContadoresTempo contadores_tempo;

	contadores_tempo.elapsed_time = 0;
	marcar_tempo_inicial(&contadores_tempo);

	AttValsIndexados* dimensoes_instanciadas = NULL;
	AttValsIndexados::iterator attribute_it;

	ResultadoInquired* dimensoes_inquired = NULL;
	ResultadoInquired* dimensoes_inquired_merged = NULL;

	AttValsIndexados* dimensoes_agregadas = NULL;

	if (!query.dimensoes_instanciadas.empty())
	{
		#if EXECUTAR_EM_CPU
		dimensoes_instanciadas = processar_dimensoes_instanciadas(fragmentos, query);
		#else
		dimensoes_instanciadas = processar_dimensoes_instanciadas(gpu_vars, fragmentos, query);
		#endif

	}
	else
	{
		dimensoes_instanciadas = new AttValsIndexados;
	}

	if ((!dimensoes_instanciadas->empty() || !query.dimensoes_aggregated.empty()) && !query.dimensoes_inquired.empty())
	{
		if (query.dimensoes_instanciadas.empty() || !dimensoes_instanciadas->empty())
		{
			// processando dimensões inquired
			dimensoes_inquired = processar_dimensoes_inquired(gpu_vars, fragmentos, query, dimensoes_instanciadas);
			//dimensoes_inquired = processar_dimensoes_inquired(fragmentos, query, dimensoes_instanciadas);
			dimensoes_inquired_merged = merge_cuboid_inquired(dimensoes_inquired);
		}
		else
		{
			dimensoes_inquired = new AttValsIndexadosPorCuboid;
			dimensoes_inquired_merged = new AttValsIndexadosPorCuboid;
		}
	}
	else
	{
		dimensoes_inquired = new AttValsIndexadosPorCuboid;
		dimensoes_inquired_merged = new AttValsIndexadosPorCuboid;
	}

	// cubo completamente materializado
	if (query.dimensoes_inquired.size( ) == query.comando.size( ))
	{
		// processando dimensões inquired
		dimensoes_inquired = coletar_todas_dimensoes_inquired(fragmentos, query);
		dimensoes_inquired_merged = merge_cuboid_inquired(dimensoes_inquired);
	}

	if (query.dimensoes_instanciadas.empty() && query.dimensoes_inquired.empty() && !query.dimensoes_aggregated.empty())
	{
		dimensoes_agregadas = processar_dimensoes_agregadas(fragmentos, query);
	}
	else
	{
		dimensoes_agregadas = new AttValsIndexados;
	}

	if (!dimensoes_agregadas->empty())
	{
		int total_rows = 0;

	#if EXIBIR_TUPLAS
		CubDim cuboid_dimensions = create_cuboid_dimensions(query.dimensoes_aggregated);
		AttVal attribute_value = dimensoes_agregadas->begin()->first;

		mensagem_executando_query(query.comando_str);
		string label = construir_attribute_label(query, attribute_value);
		std::cout << "\n    " + label << " : " << (*dimensoes_agregadas)[attribute_value].size();

	#if EXECUTAR_HAN
		(*detector_colisoes)[label] = (*dimensoes_agregadas)[attribute_value].size();
	#endif

	#endif

		total_rows++;

		marcar_tempo_final(&contadores_tempo);
		contabilizar_tempo_gasto(&contadores_tempo);

		mensagem_resultado_query(total_rows, contadores_tempo.elapsed_time);

	}
	// Instanciadas
	else if (!query.dimensoes_instanciadas.empty() && query.dimensoes_inquired.empty())
	{
		int total_rows = 0;

		AttValsIndexados::iterator it_attribute_value;
		CubDim cuboid_dimensions = create_cuboid_dimensions(query.dimensoes_instanciadas);
		it_attribute_value = dimensoes_instanciadas->begin();

		mensagem_executando_query(query.comando_str);

		if (!dimensoes_instanciadas->empty())
		{
			total_rows++;

#if EXIBIR_TUPLAS
			vector<int>* tuple_ids = &((*dimensoes_instanciadas)[it_attribute_value->first]);
			string label = construir_attribute_label(query, it_attribute_value->first);
			vector<int> tmp = dimensoes_instanciadas->begin()->first.values;
			std::cout << "\n    " << label << " : " << tuple_ids->size();

#if EXECUTAR_HAN
			(*detector_colisoes)[label] = tuple_ids->size();
#endif
#endif

		}

		marcar_tempo_final(&contadores_tempo);
		contabilizar_tempo_gasto(&contadores_tempo);

		mensagem_resultado_query(total_rows, contadores_tempo.elapsed_time);

	}
	// Inquired
	else if (!dimensoes_inquired_merged->empty() && query.dimensoes_instanciadas.empty())
	{
		int total_rows;
		ResultadoInquired::iterator it_cuboid_dimensions;
		it_cuboid_dimensions = dimensoes_inquired_merged->begin();

		AttValsIndexados::iterator it_attribute_value;
		it_attribute_value = ((*dimensoes_inquired_merged)[it_cuboid_dimensions->first]).begin( );

		total_rows = ((*dimensoes_inquired_merged)[it_cuboid_dimensions->first]).size( );

		mensagem_executando_query(query.comando_str);

#if EXIBIR_TUPLAS
		while (it_attribute_value != ((*dimensoes_inquired_merged)[it_cuboid_dimensions->first]).end())
		{
			vector<int>* tuple_ids = &((*dimensoes_inquired_merged)[it_cuboid_dimensions->first][it_attribute_value->first]);

			string label = construir_attribute_label(query, it_attribute_value->first);
			std::cout << "\n    " << label << " : " << tuple_ids->size();

#if EXECUTAR_HAN
			(*detector_colisoes)[label] = tuple_ids->size();
#endif
			it_attribute_value++;
		}
#endif

		marcar_tempo_final(&contadores_tempo);
		contabilizar_tempo_gasto(&contadores_tempo);

		mensagem_resultado_query(total_rows, contadores_tempo.elapsed_time);

	}
	// Consulta
	else if (query.dimensoes_instanciadas.empty() && (query.dimensoes_inquired.empty() || (!query.dimensoes_inquired.empty() && query.dimensoes_aggregated.empty())))
	{
		std::cout << "\nInvalid command.\n";
	}
	// Inquired + Instanciadas
	else if (!dimensoes_inquired_merged->empty() && !query.dimensoes_instanciadas.empty())
	{
		int total_rows = 0;

		mensagem_executando_query(query.comando_str);

		#if EXECUTAR_EM_CPU
		ResultadoInquired::iterator it_cuboid_dimensions;

		it_cuboid_dimensions = dimensoes_inquired_merged->begin();

		AttValsIndexados::iterator it_attribute_value;

		it_attribute_value = ((*dimensoes_inquired_merged)[it_cuboid_dimensions->first]).begin();
		total_rows = ((*dimensoes_inquired_merged)[it_cuboid_dimensions->first]).size();

		#if EXIBIR_TUPLAS
		while (it_attribute_value != ((*dimensoes_inquired_merged)[it_cuboid_dimensions->first]).end())
		{

			vector<int>* tuple_ids = &((*dimensoes_inquired_merged)[it_cuboid_dimensions->first][it_attribute_value->first]);
			string label = construir_attribute_label(query, it_attribute_value->first);
			std::cout << "\n    " << label << " : " << tuple_ids->size();

			#if EXECUTAR_HAN
			(*detector_colisoes)[label] = tuple_ids->size();
			#endif

			it_attribute_value++;
		}
		#endif
		#else
		// Executando em GPU



		#endif

		marcar_tempo_final(&contadores_tempo);
		contabilizar_tempo_gasto(&contadores_tempo);

		mensagem_resultado_query(total_rows, contadores_tempo.elapsed_time);

	}

	if (dimensoes_instanciadas != NULL)
		delete dimensoes_instanciadas;

	if (dimensoes_inquired != NULL)
		delete dimensoes_inquired;

	if (dimensoes_inquired_merged != NULL)
		delete dimensoes_inquired_merged;

	if (dimensoes_agregadas != NULL)
		delete dimensoes_agregadas;
}

inline bool create_query(Fragmentos* fragmentos, string query_str, Query* query_retorno)
{
	vector<string> parametros = split(query_str, ':');
	string inquired_symbol = string({ INQUIRED_CARACTER });
	string aggregated_symbol = string({ AGGREGATED_CARACTER });
	set<int>::iterator it = query_retorno->set_dimensoes_inquired.begin();

	query_retorno->comando_str = query_str;

	for (int i = 0; i < parametros.size( ); i++)
	{
		// se é dimensao inquired
		if (parametros[i].compare(inquired_symbol) == 0)
		{
			query_retorno->dimensoes_inquired.push_back(i + 1);
			query_retorno->set_dimensoes_inquired.insert(it, i + 1);
			query_retorno->mapeamento_colunas.push_back(DIMENSAO_INQUIRED);
		}
		// se é dimensao aggregated
		else if (parametros[i].compare(aggregated_symbol) == 0)
		{
			query_retorno->dimensoes_aggregated.push_back(i + 1);
			query_retorno->mapeamento_colunas.push_back(DIMENSAO_AGREGADA);
		}
		// se não é dimensao instanciada
		else if (parametros[i].length( ) > 0)
		{
			query_retorno->dimensoes_instanciadas.push_back(i + 1);
			query_retorno->mapeamento_colunas.push_back(DIMENSAO_INSTANCIADA);
			query_retorno->valores_instanciados.push_back(stoi(parametros[i]));
		}
		// faça nada ainda
		else
		{
			return false;
		}
	}

	query_retorno->comando = parametros;

	return true;

}

inline AttValsIndexados* processar_dimensoes_agregadas(Fragmentos* fragmentos, Query query)
{
	vector<int> attribute_values_vector;
	AttValsIndexados* resultado = new AttValsIndexados;

	for (register int d = 0; d < query.dimensoes_aggregated.size(); d++)
		attribute_values_vector.push_back(ALL);

	AttVal attribute_value = create_attribute_value(attribute_values_vector);

	(*resultado)[attribute_value] = fragmentos->all_ids;

	return resultado;
}

inline AttValsIndexados* processar_dimensoes_instanciadas(Fragmentos* fragmentos, Query query)
{
	AttValsIndexados* resultado;
	AttVal last_attribute_value;
	CubDim last_cuboid_dimensions;

	resultado = new AttValsIndexados;

	for (int d = 0; d < query.dimensoes_instanciadas.size(); d++)
	{
		CubDim cuboid_dimensions = create_cuboid_dimensions( query.dimensoes_instanciadas[d] );
		AttVal attribute_value = create_attribute_value( query.valores_instanciados[d] );

		if (d == 0)
		{
			if (fragmentos->cuboids[cuboid_dimensions].values.find(attribute_value) != fragmentos->cuboids[cuboid_dimensions].values.end())
			{
				last_attribute_value = attribute_value;
				last_cuboid_dimensions = cuboid_dimensions;
				(*resultado)[attribute_value] = fragmentos->cuboids[cuboid_dimensions].values[attribute_value];
			}
			else
				return resultado;
		}
		else
		{
			if (fragmentos->cuboids[cuboid_dimensions].values.find(attribute_value) != fragmentos->cuboids[cuboid_dimensions].values.end())
			{
				AttVal attribute_value_all = create_attribute_value( ALL );
				vector<int>* tuple_ids_1 = &(fragmentos->cuboids[cuboid_dimensions].values[attribute_value]);
				vector<int>* tuple_ids_2 = &(*resultado)[last_attribute_value];
				vector<int> retorno;

				if (attribute_values_equals(attribute_value_all, attribute_value))
				{
					retorno = *tuple_ids_2;
				}
				else if (attribute_values_equals(attribute_value_all, last_attribute_value))
				{
					retorno = *tuple_ids_1;
				}
				else
				{
					retorno = set_intersection(tuple_ids_1, tuple_ids_2);
				}

				if (!retorno.empty())
				{

					AttVal novo_attribute_value = concatenate_attribute_value
					(
						last_attribute_value,
						attribute_value,
						last_cuboid_dimensions,
						cuboid_dimensions
					);

					last_cuboid_dimensions = concatenate_cuboid_dimensions(last_cuboid_dimensions, cuboid_dimensions);
					last_attribute_value = novo_attribute_value;

					resultado->clear();

					(*resultado)[novo_attribute_value] = retorno;

				}
				else
				{
					resultado->clear();
					return resultado;
				}

			}
			else
			{
				resultado->clear();
				return resultado;
			}
		}
	}

	return resultado;

}

inline ResultadoInquired *coletar_todas_dimensoes_inquired(Fragmentos* fragmentos, Query query)
{
	ResultadoInquired* resultado;
	map<CubDim, Cuboid, CompCubDim>::iterator fragmentos_iterator;

	resultado = new map < CubDim, AttValsIndexados, CompCubDim>;
	fragmentos_iterator = fragmentos->cuboids.begin();

	while (fragmentos_iterator != fragmentos->cuboids.end())
	{
		(*resultado)[fragmentos_iterator->first] = fragmentos->cuboids[fragmentos_iterator->first].values;
		fragmentos_iterator++;
	}

	return resultado;

}


inline AttValsIndexados* processar_dimensoes_instanciadas(gpu::GpuVars* gpu_vars, Fragmentos* fragmentos, Query query)
{
	AttValsIndexados* resultado;
	AttVal last_attribute_value;
	CubDim last_cuboid_dimensions;

	int grid_size;
	int block_size;

	resultado = new AttValsIndexados;

	if (query.dimensoes_instanciadas.size() == 1)
	{
		CubDim cub_dim = create_cuboid_dimensions(query.dimensoes_instanciadas[0]);
		AttVal att_val = create_attribute_value(query.valores_instanciados[0]);

		(*resultado)[att_val] = fragmentos->cuboids[cub_dim].values[att_val];
		vector<int>* fonte = &(fragmentos->cuboids[cub_dim].values[att_val]);

		gpu::copiar_cpu_para_gpu(gpu_vars->inicio_memoria, &((*fonte)[0]), fonte->size( ) * sizeof(int));
		gpu_vars->total_numeros_na_gpu = fonte->size();

		return resultado;
	}

	#if KERNELS
	gpu_vars->total_kernels_simultaneos = KERNELS;
	#else
	gpu_vars->total_kernels_simultaneos = 1;
	#endif

	int threshold = query.dimensoes_instanciadas.size() / gpu_vars->total_kernels_simultaneos;

	#if DEBUG_STREAM
	cudaEvent_t inicio_event;
	cudaEvent_t fim_event;
	float ms;
	#endif

	gpu::KernelVars* kernel_vars = gpu_vars->kernel_vars;

	vector<int> inicio_chunks;
	vector<int> fim_chunks;

	int pos_inicio_chunk = 0;
	int s = 0;

	gpu::particionar_memoria(gpu_vars);

	for (register int d = 0; d < query.dimensoes_instanciadas.size(); d++)
	{
		if (s + 1 == gpu_vars->total_kernels_simultaneos)
			threshold += query.dimensoes_instanciadas.size() % gpu_vars->total_kernels_simultaneos;

		CubDim cub_dim = create_cuboid_dimensions(query.dimensoes_instanciadas[d]);
		AttVal att_val = create_attribute_value(query.valores_instanciados[d]);

		vector<int>* tid_list = &fragmentos->cuboids[cub_dim].values[att_val];
		pos_inicio_chunk += tid_list->size();

		if (!gpu_vars->kernel_vars[s].total_elementos_vetores_menores)
		{
			gpu_vars->kernel_vars[s].inicio_espaco_busca = gpu_vars->kernel_vars[s].inicio_vetores_menores + tid_list->size();
			gpu_vars->kernel_vars[s].inicio_vetores_maiores = gpu_vars->kernel_vars[s].inicio_vetores_menores + tid_list->size();		
			gpu_vars->kernel_vars[s].total_elementos_vetores_menores = tid_list->size();

			gpu::add_buffer(gpu_vars, &(*tid_list)[0], tid_list->size());
		}
		else
		{
			gpu_vars->kernel_vars[s].total_intersecoes++;
			gpu_vars->kernel_vars[s].total_elementos_vetores_maiores += tid_list->size();

			gpu::add_buffer(gpu_vars, &(*tid_list)[0], tid_list->size());

			inicio_chunks.push_back(pos_inicio_chunk - tid_list->size());
			fim_chunks.push_back(pos_inicio_chunk);

			if (gpu_vars->kernel_vars[s].total_intersecoes == threshold || d + 1 == query.dimensoes_instanciadas.size())
			{
				gpu::add_buffer(gpu_vars, &(inicio_chunks)[0], inicio_chunks.size());
				gpu::add_buffer(gpu_vars, &(fim_chunks)[0], fim_chunks.size());

				gpu_vars->kernel_vars[s].resultado = gpu_vars->kernel_vars[s].inicio_espaco_busca +
					gpu_vars->kernel_vars[s].total_elementos_vetores_maiores + (gpu_vars->kernel_vars[s].total_intersecoes * 2);

				#if DEBUG_STREAM
				cudaEventCreate(&inicio_event);
				cudaEventCreate(&fim_event);
				cudaEventRecord(inicio_event);
				#endif

				cudaMemcpyAsync(gpu_vars->kernel_vars[s].inicio_vetores_menores, gpu_vars->buffer.vetor,
					gpu_vars->buffer.pos_buffer * sizeof(int), cudaMemcpyHostToDevice, gpu_vars->kernel_vars[s].stream);

				#if DEBUG_STREAM
				cudaEventRecord(fim_event);
				cudaEventSynchronize(fim_event);
				cudaEventElapsedTime(&ms, inicio_event, fim_event);

				std::cout << "\nStream " << s << " - Tempo gasto com a Copia H2D: " << ms << " ms";

				cudaEventDestroy(inicio_event);
				cudaEventDestroy(fim_event);
				#endif

				gpu_vars->buffer.pos_buffer = 0;
				s++;

				if (s == gpu_vars->total_kernels_simultaneos)
				{
					s = 0;
				}
				else
				{
					gpu_vars->kernel_vars[s].total_intersecoes = 0;
					gpu_vars->kernel_vars[s].total_elementos_vetores_maiores = 0;
					gpu_vars->kernel_vars[s].total_elementos_vetores_menores = 0;

					inicio_chunks.clear();
					fim_chunks.clear();

					pos_inicio_chunk = 0;
				}

			}

		}

	}

	for (int s = 0; s < gpu_vars->total_kernels_simultaneos; s++)
	{
		kernel_vars[s].resultado =
			kernel_vars[s].inicio_vetores_menores + (kernel_vars[s].total_intersecoes * 2) +
			kernel_vars[s].total_elementos_vetores_menores + kernel_vars[s].total_elementos_vetores_maiores;

		cudaOccupancyMaxPotentialBlockSize(&grid_size, &block_size, gpu::intersection_with_chunks);
		grid_size = (kernel_vars[s].total_elementos_vetores_menores/block_size) + 1;

		// std::cout << "\n\nTotal de threads: " << kernel_vars[s].total_elementos_vetores_menores << "\n\n";

		#if DEBUG_STREAM
		cudaEventCreate(&inicio_event);
		cudaEventCreate(&fim_event);
		cudaEventRecord(inicio_event);
		#endif

		gpu::intersection_with_chunks << <grid_size, block_size, 0, kernel_vars[s].stream >> >(kernel_vars[s]);

		#if DEBUG_STREAM
		cudaEventRecord(fim_event);
		cudaEventSynchronize(fim_event);
		cudaEventElapsedTime(&ms, inicio_event, fim_event);

		std::cout << "\nStream " << s << " - Tempo gasto com a intersecao: " << ms << " ms";

		cudaEventDestroy(inicio_event);
		cudaEventDestroy(fim_event);
		#endif

	}
	
	gpu::KernelVars kernel_vars_merge = kernel_vars[0];

	int tamanho_resposta;
	int total_elementos_copiados = 0;
	int conjunto_vazio = 0;

	kernel_vars_merge.total_intersecoes = 0;
	kernel_vars_merge.total_elementos_vetores_maiores = 0;
	kernel_vars_merge.total_elementos_vetores_menores = 0;

	pos_inicio_chunk = 0;

	inicio_chunks.clear();
	fim_chunks.clear();

	// copy D2H
	for (int s = 0; s < gpu_vars->total_kernels_simultaneos; s++)
	{
		if (gpu_vars->total_kernels_simultaneos == 1)
		{
			gpu::gpu_stream_compaction
			(
				kernel_vars[s].resultado,
				kernel_vars_merge.inicio_vetores_menores + total_elementos_copiados,
				kernel_vars[s].total_elementos_vetores_menores,
				&tamanho_resposta
			);
		}
		else
		{
			
		}

		if (tamanho_resposta)
		{
			total_elementos_copiados += tamanho_resposta;

			if (!s)
			{
				kernel_vars_merge.inicio_espaco_busca = kernel_vars_merge.inicio_vetores_menores + tamanho_resposta;
				kernel_vars_merge.inicio_vetores_maiores = kernel_vars_merge.inicio_vetores_menores + tamanho_resposta;
				kernel_vars_merge.total_elementos_vetores_menores = tamanho_resposta;
			}
			else
			{
				inicio_chunks.push_back(pos_inicio_chunk);
				fim_chunks.push_back(pos_inicio_chunk + tamanho_resposta);

				kernel_vars_merge.total_intersecoes++;
				kernel_vars_merge.total_elementos_vetores_maiores += tamanho_resposta;
			}

			pos_inicio_chunk += tamanho_resposta;
			cudaStreamDestroy(kernel_vars[s].stream);

		}
		else
		{
			conjunto_vazio = 1;

			while (s < gpu_vars->total_kernels_simultaneos)
			{
				cudaStreamDestroy(kernel_vars[s].stream);
				s++;
			}

		}

	}

	if (gpu_vars->total_kernels_simultaneos > 1)
	{
		if (!conjunto_vazio)
		{

			gpu::add_buffer(gpu_vars, &(inicio_chunks[0]), inicio_chunks.size());
			gpu::add_buffer(gpu_vars, &(fim_chunks[0]), fim_chunks.size());

			cudaMemcpy(kernel_vars_merge.inicio_vetores_menores + total_elementos_copiados,
				gpu_vars->buffer.vetor, gpu_vars->buffer.pos_buffer * 4, cudaMemcpyHostToDevice);

			gpu_vars->buffer.pos_buffer = 0;

			kernel_vars_merge.resultado = kernel_vars_merge.inicio_vetores_menores + total_elementos_copiados + (kernel_vars_merge.total_intersecoes * 2);

			cudaOccupancyMaxPotentialBlockSize(&grid_size, &block_size, gpu::intersection_with_chunks);
			grid_size = (kernel_vars_merge.total_elementos_vetores_menores / block_size) + 1;

			gpu::intersection_with_chunks << < grid_size, block_size >> > (kernel_vars_merge);

			gpu::gpu_stream_compaction(kernel_vars_merge.resultado, kernel_vars_merge.inicio_vetores_menores,
				kernel_vars_merge.total_elementos_vetores_menores, &tamanho_resposta);

			cudaMemcpy(gpu_vars->buffer.vetor, kernel_vars_merge.inicio_vetores_menores, tamanho_resposta * sizeof(int), cudaMemcpyDeviceToHost);

			vector<int> tid_resposta(gpu_vars->buffer.vetor, gpu_vars->buffer.vetor + tamanho_resposta);
			AttVal novo_av = create_attribute_value(query.valores_instanciados);

			(*resultado)[novo_av] = tid_resposta;

		}
	}
	else
	{
		if (!conjunto_vazio)
		{
			cudaMemcpy(gpu_vars->buffer.vetor, kernel_vars[0].inicio_vetores_menores, tamanho_resposta * sizeof(int), cudaMemcpyDeviceToHost);

			vector<int> tid_resposta(gpu_vars->buffer.vetor, gpu_vars->buffer.vetor + tamanho_resposta);
			AttVal novo_av = create_attribute_value(query.valores_instanciados);

			(*resultado)[novo_av] = tid_resposta;
		}
	}

	// adicionar contador numeros da GPU
	int freq_atributos_dim_instanciadas = (*resultado)[resultado->begin()->first].size();
	gpu::set_numeros_na_gpu(gpu_vars, freq_atributos_dim_instanciadas);

	return resultado;

}

ResultadoInquired* processar_dimensoes_inquired(Fragmentos* fragmentos, Query query, AttValsIndexados* dimensoes_instanciadas)
{
	ResultadoInquired* resultado;
	vector<int>* tuple_id_dimensoes_instanciadas;

	if (!query.dimensoes_instanciadas.empty())
		tuple_id_dimensoes_instanciadas = &((*dimensoes_instanciadas)[dimensoes_instanciadas->begin()->first]);

	resultado = new ResultadoInquired;

	for (int d = 0; d < query.dimensoes_inquired.size(); d++)
	{
		CubDim cuboid_dimensions = create_cuboid_dimensions( query.dimensoes_inquired[d] );
		AttValsIndexados::iterator it_attribute_value;

		it_attribute_value = fragmentos->cuboids[cuboid_dimensions].values.begin( );

		while (it_attribute_value != fragmentos->cuboids[cuboid_dimensions].values.end())
		{
			if (!query.dimensoes_instanciadas.empty())
			{
				AttVal attribute_value_all = create_attribute_value( ALL );
				vector<int>* tuple_id_list = &(fragmentos->cuboids[cuboid_dimensions].values[it_attribute_value->first]);
				vector<int> retorno;

				if (attribute_values_equals(attribute_value_all, it_attribute_value->first))
				{
					retorno = *tuple_id_dimensoes_instanciadas;
				}
				else if (attribute_values_equals(attribute_value_all, dimensoes_instanciadas->begin()->first))
				{
					retorno = *tuple_id_list;
				}
				else
				{
					retorno = set_intersection(tuple_id_dimensoes_instanciadas, tuple_id_list);
				}

				if (!retorno.empty())
				{
					if (resultado->find(cuboid_dimensions) == resultado->end())
						(*resultado)[cuboid_dimensions] = AttValsIndexados();

					(*resultado)[cuboid_dimensions][it_attribute_value->first] = retorno;
				}

			}
			else
			{
				if (resultado->find(cuboid_dimensions) == resultado->end())
					(*resultado)[cuboid_dimensions] = AttValsIndexados();

				(*resultado)[cuboid_dimensions][it_attribute_value->first] = fragmentos->cuboids[cuboid_dimensions].values[it_attribute_value->first];

			}

			it_attribute_value++;

		}

	}

	return resultado;

}

ResultadoInquired*
merge_cuboid_inquired(ResultadoInquired* dimensoes_inquired)
{
	map<CubDim, AttValsIndexados>::iterator it_cuboid_dimensions = dimensoes_inquired->begin();
	ResultadoInquired* resultado;

	resultado = new AttValsIndexadosPorCuboid;

	// Iterando dimensoes de dimensions_inquired
	while (it_cuboid_dimensions != dimensoes_inquired->end())
	{
		if (it_cuboid_dimensions == dimensoes_inquired->begin())
		{
			(*resultado)[it_cuboid_dimensions->first] = (*dimensoes_inquired)[it_cuboid_dimensions->first];
		}
		else
		{
			AttValsIndexados::iterator iterator_attribute_value;
			iterator_attribute_value = (*dimensoes_inquired)[it_cuboid_dimensions->first].begin();

			ResultadoInquired* resultado_tmp;
			resultado_tmp = new AttValsIndexadosPorCuboid;

			// iterando attribute values de dimensoes inquired;
			while (iterator_attribute_value != (*dimensoes_inquired)[it_cuboid_dimensions->first].end())
			{
				ResultadoInquired::iterator cuboid_resultado;
				cuboid_resultado = resultado->begin();

				AttValsIndexados::iterator iterator_attribute_value_resultado;
				iterator_attribute_value_resultado = (*resultado)[cuboid_resultado->first].begin();

				vector<int>* tuple_ids_1 = &(*dimensoes_inquired)[it_cuboid_dimensions->first][iterator_attribute_value->first];

				// Iterando attribute values de resultado
				while (iterator_attribute_value_resultado != (*resultado)[cuboid_resultado->first].end())
				{
					AttVal attribute_value_all = create_attribute_value(ALL);
					vector<int>* tuple_ids_2 = &(*resultado)[cuboid_resultado->first][iterator_attribute_value_resultado->first];
					vector<int> resultado;

					if (attribute_values_equals(attribute_value_all, iterator_attribute_value->first))
					{
						resultado = *tuple_ids_2;
					}
					else if (attribute_values_equals(attribute_value_all, iterator_attribute_value_resultado->first))
					{
						resultado = *tuple_ids_1;
					}
					else
					{
						resultado = set_intersection(tuple_ids_1, tuple_ids_2);
					}

					if (!resultado.empty())
					{
						CubDim novo_cuboid_dimensions = concatenate_cuboid_dimensions(cuboid_resultado->first, it_cuboid_dimensions->first);
						AttVal novo_attribute_values = concatenate_attribute_value
						(
							iterator_attribute_value_resultado->first,
							iterator_attribute_value->first,
							cuboid_resultado->first,
							it_cuboid_dimensions->first
						);

						(*resultado_tmp)[novo_cuboid_dimensions].insert({ novo_attribute_values, resultado });

					}

					iterator_attribute_value_resultado++;
				}

				iterator_attribute_value++;

			}

			if (!resultado_tmp->empty())
			{
				delete resultado;
				resultado = resultado_tmp;
				resultado_tmp = NULL;
			}
			else
			{
				delete resultado;
			}

		}

		it_cuboid_dimensions++;
	}

	return resultado;
}

void limpar_query(Query* query)
{
	query->comando.clear();
	query->dimensoes_inquired.clear();
	query->dimensoes_instanciadas.clear();
	query->mapeamento_colunas.clear();
	query->valores_instanciados.clear();
	query->dimensoes_aggregated.clear();
}

void inline indexar_tuplas(Query query, map<string, long>* detector_colisoes, AttValsIndexados resultado)
{
	AttValsIndexados::iterator iterador;

	iterador = resultado.begin();

	while (iterador != resultado.end())
	{
		(*detector_colisoes)[construir_attribute_label(query, iterador->first)] = resultado[iterador->first].size();
		iterador++;
	}
}

string inline construir_attribute_label(Query query, AttVal attribute_value)
{
	register int comando_size = query.comando.size();
	register string label = "";
	register int ia = 0;

	for (register int d = 0; d < comando_size; d++)
	{
		if (query.comando[d].compare("?") == 0)
		{
			if (attribute_value.values[ia] != ALL)
				label += to_string(attribute_value.values[ia]);
			else
				label += "*";
			ia++;
		}
		else if (query.comando[d].compare("*") == 0)
		{
			label += "*";
		}
		else
			label += query.comando[d];

		if (d + 1 < comando_size)
			label += " ";
	}
	
	return label;
}

void thread_gerar_cuboids(vector<int> valores, vector<vector<int>>* labels, bool* terminado)
{
	*labels = gerar_permutacoes(valores);
	*terminado = true;
}

inline int descobrir_cuboid_celula(vector<pair<int, int>>* deslocamentos_cuboids_vector_pairs, int fim)
{
	register int i = 0;

	while (i < deslocamentos_cuboids_vector_pairs->size( ))
	{
		// cout << "\n\nSecond: " << (*deslocamentos_cuboids_vector_pairs)[i].second;
		// cout << "\nFim: " << fim;

		if (fim == (*deslocamentos_cuboids_vector_pairs)[i].second)
			return i;
		else if (fim > (*deslocamentos_cuboids_vector_pairs)[i].second)
			return i + 1;

		i++;
	}

	return i;
}

inline ResultadoInquired*
processar_dimensoes_inquired(gpu::GpuVars* gpu_vars, Fragmentos* fragmentos, Query query, AttValsIndexados* dimensoes_instanciadas)
{
	ResultadoInquired* resultado;

	map < int, vector<int>> inq_dim_por_card;
	vector<int> cardinalidades_inquired_ordenadas;
	string label_dimensoes_instanciadas = "INST";
	vector<int> vector_dimensoes_instanciadas = query.dimensoes_instanciadas;

	vector<CubDim> deslocamentos_cuboids_vector;
	vector <pair<int, int>> deslocamentos_cuboids_vector_pairs;
	set<CubDim, CompCubDim> deslocamentos_cuboids_set;

	resultado = new ResultadoInquired;

	for (register int i = 0; i < query.dimensoes_inquired.size( ); i++)
	{
		int inq_dim = query.dimensoes_inquired[i];
		int i_card_inq_dim = fragmentos->cardinalidades_por_coluna[inq_dim - 1];

		if (inq_dim_por_card.find(i_card_inq_dim) == inq_dim_por_card.end( ))
		{
			vector<int> vetor;
			inq_dim_por_card[i_card_inq_dim] = vetor;
			cardinalidades_inquired_ordenadas.push_back(i_card_inq_dim);
		}

		inq_dim_por_card[i_card_inq_dim].push_back(inq_dim);

	}

	std::sort(cardinalidades_inquired_ordenadas.begin(), cardinalidades_inquired_ordenadas.end( ));
	std::sort(vector_dimensoes_instanciadas.begin( ), vector_dimensoes_instanciadas.end( ));

	std::cout << "\nCardinalidades INQUIRED ordenadas: " << label_tuple_id_list(&cardinalidades_inquired_ordenadas);
	std::cout << "\nDimensoes instanciadas: " << label_tuple_id_list(&vector_dimensoes_instanciadas);
	getchar( );

	vector<vector<int>> permutacoes;
	bool thread_terminada = false;

	std::thread thread_cuboids(thread_gerar_cuboids, cardinalidades_inquired_ordenadas, &permutacoes, &thread_terminada);

	vector<int> tids_instanciadas = (*dimensoes_instanciadas)[dimensoes_instanciadas->begin( )->first];

	vector<int> inicio;
	vector<int> fim;

	int freq_atributos_dim_instanciadas = tids_instanciadas.size( );
	map<AttVal, vector<int>>::iterator it_attval_ref;
	AttVal attval_all = create_attribute_value(ALL);
	int freq_valor_atributo_inq = 0;

	// total tids adicionadas a gpu
	int total_tids_attval_inquired = 0;
	int replicas;

	vector<string> celulas;
	CubDim cuboid_dimensoes_instanciadas = create_cuboid_dimensions(vector_dimensoes_instanciadas);

	celulas.push_back(label_dimensoes_instanciadas);

	pair<int, int> pair_deslocamento = pair<int, int>(0, freq_atributos_dim_instanciadas);

	deslocamentos_cuboids_set.insert(cuboid_dimensoes_instanciadas);
	deslocamentos_cuboids_vector_pairs.push_back(pair_deslocamento);
	deslocamentos_cuboids_vector.push_back(cuboid_dimensoes_instanciadas);

	vector<int> deslocamento_celulas_inicio;
	vector<int> deslocamento_celulas_fim;

	// poolling na variavel
	while (!thread_terminada);

	thread_cuboids.join( );

	int pair_deslocamento_first;
	int pair_deslocamento_last;

	// iterando cardinalidades inquired já ordenadas
	for (int i_card = 0; i_card < cardinalidades_inquired_ordenadas.size( ); i_card++)
	{
		int cardinalidade_dim_inquired = cardinalidades_inquired_ordenadas[i_card];

		// obtendo dimensões inquired com a cardinalidade em questão
		vector<int>* dimensoes_inq = &(inq_dim_por_card[cardinalidade_dim_inquired]);

		int i_dim_inq = 0;

		CubDim cuboid_dimensao_processada = create_cuboid_dimensions(cuboid_dimensoes_instanciadas, (*dimensoes_inq)[i_dim_inq]);
		pair_deslocamento_first = freq_atributos_dim_instanciadas;
		pair_deslocamento_last;

		// Processando primeira dimensão inquired
		if (!i_card && !i_dim_inq)
		{
			CubDim cub_dim = create_cuboid_dimensions((*dimensoes_inq)[i_dim_inq]);
			Cuboid cuboid = *get_cuboid(fragmentos, cub_dim);
			map<AttVal, vector<int>>::iterator it_attval_ref = get_attval_iterator(&cuboid);

			// se trocou de dimensao
			total_tids_attval_inquired = 0;

			#if DEBUG
			std::cout << "\n\nDIMENSAO: " << (*dimensoes_inq)[i_dim_inq];
			std::cout << "\nCardinalidade " << cardinalidade_dim_inquired;
			#endif

			AttVal attval;
			map<AttVal, vector<int>>::iterator it_attval = get_attval_iterator(&cuboid);

			#if DEBUG
			std::cout << "\n\nAtributos da dimensao inquired " << (*dimensoes_inq)[i_dim_inq];
			#endif

			int max_replicas_gpu = gpu::calcular_replicacao(gpu_vars, &cuboid, &it_attval_ref, freq_atributos_dim_instanciadas, TAMANHO_FIXO, freq_atributos_dim_instanciadas);
			replicas = gpu::replicar_tids(gpu_vars, max_replicas_gpu, freq_atributos_dim_instanciadas);

			#if DEBUG
			gpu::exibir_memoria(gpu_vars->inicio_memoria, freq_atributos_dim_instanciadas * max_replicas_gpu);
			std::cout << "\n\nTotal de replicas feitas: " << replicas << "\n";
			#endif

			int cont_attval_inq = 0;

			// iterando atributos
			while (it_attval != cuboid.values.end())
			{
				vector<int> tids_attval_dim_inq = get_tids(&cuboid, it_attval->first);

				// se nao e ALL
				if (attribute_values_equals(it_attval->first, attval_all) == false)
				{
					#if DEBUG
					std::cout << "\n\n[" << label_attribute_values(it_attval->first) << "] => " << label_tuple_id_list(&tids_attval_dim_inq) << "\nCOUNT: " << tids_attval_dim_inq.size();
					#endif

					freq_valor_atributo_inq = tids_attval_dim_inq.size( );

					#if DEBUG
					std::cout << "\nTotal bytes ocupados: " << gpu_vars->total_bytes_ocupados;
					std::cout << "\nTotal numeros na GPU: " << gpu_vars->total_numeros_na_gpu;
					#endif

					bool gpu_cheia = (it_attval == it_attval_ref);

					if (!gpu_cheia && cont_attval_inq < cardinalidade_dim_inquired)
					{
						#if DEBUG
						std::cout << "\n\n\nA memoria nao encheu. A GPU tem " << ((gpu::gpu_bytes_livres(gpu_vars) / 1024.0) / 1024.0) << " MB livres...";
						#endif
						// adicionando mais um atributo ao buffer
						total_tids_attval_inquired += freq_valor_atributo_inq;

						inicio.push_back(total_tids_attval_inquired - freq_valor_atributo_inq);
						fim.push_back(total_tids_attval_inquired);

						// vendo se o buffer irá saturar
						if (gpu::buffer_cheio(&gpu_vars->buffer, freq_valor_atributo_inq) == true)
						{
							#if DEBUG
							std::cout << "\nBuffer encheu! Hora de copiar para a GPU!";
							#endif
							// flush buffer
							gpu::esvaziar_buffer(gpu_vars);
						}
					
						gpu::add_buffer(gpu_vars, &(tids_attval_dim_inq[0]), tids_attval_dim_inq.size( ));

						#if DEBUG
						std::cout << "\n\nAdicionando TIDs de " << label_attribute_values(it_attval->first) << " ao buffer...";
						std::cout << "\nTotal elementos no buffer: " << gpu_vars->buffer.pos_buffer;
						#endif

						// incrementando contador de valor de atributos percorridos
						cont_attval_inq++;

						// ir pro proximo valor de atributo
						it_attval++;

					}

					if (gpu_cheia || cont_attval_inq == cardinalidade_dim_inquired)
					{
						if (gpu_cheia)
						{
							// ainda nao implementa esta parte
							printf("\nA GPU ENCHEU! ABORTANDO...");
							exit(0);
						}

						gpu::add_buffer(gpu_vars, &(inicio[0]), inicio.size( ));
						gpu::add_buffer(gpu_vars, &(fim[0]), fim.size( ));

						gpu::esvaziar_buffer(gpu_vars);

						int tamanho_bloco;
						int tamanho_grid;

						int total_threads = replicas * freq_atributos_dim_instanciadas;
						gpu::balancear_carga(gpu_vars, total_threads, &tamanho_bloco, &tamanho_grid);

						inicio.clear( );
						fim.clear( );

						gpu::KernelVarsSubcube kernel_vars;

						kernel_vars.replicas = replicas;
						kernel_vars.cardinalidade = cardinalidade_dim_inquired;
						kernel_vars.inicio_memoria = gpu_vars->inicio_memoria;
						kernel_vars.freq_atributos_dim_instanciadas = freq_atributos_dim_instanciadas;
						kernel_vars.inicio_espaco_busca = kernel_vars.inicio_memoria + (kernel_vars.replicas * kernel_vars.freq_atributos_dim_instanciadas);
						kernel_vars.offset_inicio = kernel_vars.inicio_espaco_busca + total_tids_attval_inquired;
						kernel_vars.offset_fim = kernel_vars.offset_inicio + kernel_vars.cardinalidade;
						kernel_vars.resposta = kernel_vars.offset_fim + kernel_vars.cardinalidade;

						int total_numeros_resposta = kernel_vars.freq_atributos_dim_instanciadas * kernel_vars.cardinalidade;

						gpu::subcubo_com_replicacao << <tamanho_grid, tamanho_bloco >> >(kernel_vars);
						cudaDeviceSynchronize( );

						// fazendo a compactação das células
						int contador_celulas = 0;
						int* inicio_celula_agregada = kernel_vars.inicio_memoria + freq_atributos_dim_instanciadas;
						int offset_celula_agregada = freq_atributos_dim_instanciadas;

						// adicionando offsets da celula base composta so de valores instanciados
						deslocamento_celulas_inicio.push_back(offset_celula_agregada - freq_atributos_dim_instanciadas);
						deslocamento_celulas_fim.push_back(offset_celula_agregada);

						gpu_vars->total_numeros_na_gpu = freq_atributos_dim_instanciadas;

						while (contador_celulas < cont_attval_inq)
						{

							if (offset_celula_agregada + freq_atributos_dim_instanciadas > freq_atributos_dim_instanciadas * replicas)
							{

							}

							int tamanho_array_compacto;
							gpu::gpu_stream_compaction(kernel_vars.resposta + (contador_celulas * freq_atributos_dim_instanciadas), inicio_celula_agregada, freq_atributos_dim_instanciadas, &tamanho_array_compacto);

							if (tamanho_array_compacto)
							{
								inicio_celula_agregada = inicio_celula_agregada + tamanho_array_compacto;
								gpu_vars->total_numeros_na_gpu += tamanho_array_compacto;

								// ----
								string valor = "." + to_string(cardinalidades_inquired_ordenadas[i_card] - contador_celulas);
								string label = label_dimensoes_instanciadas + ", " + to_string((*dimensoes_inq)[i_dim_inq]) + valor;
								celulas.push_back(label);
								// ----

								offset_celula_agregada += tamanho_array_compacto;

								// adicionando mais uma celula do sub-cubo
								deslocamento_celulas_inicio.push_back(offset_celula_agregada - tamanho_array_compacto);
								deslocamento_celulas_fim.push_back(offset_celula_agregada);

							}

							contador_celulas++;

						}

						pair_deslocamento_last = gpu_vars->total_numeros_na_gpu;

					}

				}
				else
				{
					// se não é ALL
					it_attval++;
				}

			}

			// se acabou esta dimensao inquired, vamos para a proxima
			i_dim_inq++;

			if (deslocamentos_cuboids_set.find(cuboid_dimensao_processada) == deslocamentos_cuboids_set.end())
			{
				pair_deslocamento_last = gpu_vars->total_numeros_na_gpu;
				pair<int, int> pair_cuboid(pair_deslocamento_first, pair_deslocamento_last);
				deslocamentos_cuboids_set.insert(cuboid_dimensao_processada);
				deslocamentos_cuboids_vector_pairs.push_back(pair_cuboid);
				deslocamentos_cuboids_vector.push_back(cuboid_dimensao_processada);
			}

		}

		// PROCESSANDO DIMENSOES INQUIRED DIFERENTES DA PRIMEIRA
		// somatorio do count de tids de todas as celulas obtidas na iteracao anterior
		int freq_tids_celulas_agregadas = gpu_vars->total_numeros_na_gpu;

		// total de celulas obtidas da outra agregacacao
		int total_celulas_derivadas = deslocamento_celulas_inicio.size( );

		while (i_dim_inq < dimensoes_inq->size( ))
		{
			int dimensao_inquired = (*dimensoes_inq)[i_dim_inq];
			CubDim cub_dim = create_cuboid_dimensions(dimensao_inquired);
			Cuboid cuboid = *get_cuboid(fragmentos, cub_dim);
			map<AttVal, vector<int>>::iterator it_attval_ref = get_attval_iterator(&cuboid);

			// se trocou de dimensao
			total_tids_attval_inquired = 0;

			AttVal attval;
			map<AttVal, vector<int>>::iterator it_attval = get_attval_iterator(&cuboid);
			
			int max_replicas_gpu = gpu::calcular_replicacao(gpu_vars, &cuboid, &it_attval_ref, freq_tids_celulas_agregadas, TAMANHO_VARIAVEL, total_celulas_derivadas);
			
			replicas = gpu::replicar_tids(gpu_vars, max_replicas_gpu, freq_tids_celulas_agregadas);

			int i_attval_inq = 0;

			// iterando atributos
			while (it_attval != cuboid.values.end( ))
			{
				vector<int> tids_attval_dim_inq = get_tids(&cuboid, it_attval->first);

				if (attribute_values_equals(it_attval->first, attval_all) == false)
				{
					#if DEBUG
					std::cout << "\n\t[" << label_attribute_values(it_attval->first) << "] => "
						<< label_tuple_id_list(&tids_attval_dim_inq) << "\nCOUNT: " << tids_attval_dim_inq.size();
					#endif

					freq_valor_atributo_inq = tids_attval_dim_inq.size( );

					bool gpu_cheia = (it_attval == it_attval_ref);

					if (!gpu_cheia && i_attval_inq < cardinalidade_dim_inquired)
					{
						// adicionando mais um atributo ao Buffer
						total_tids_attval_inquired += freq_valor_atributo_inq;

						// deslocamentos do espaco de busca dos valores de atributo
						inicio.push_back(total_tids_attval_inquired - freq_valor_atributo_inq);
						fim.push_back(total_tids_attval_inquired);

						// vendo se o buffer irá saturar
						if (gpu::buffer_cheio(&gpu_vars->buffer, freq_valor_atributo_inq) == true)
						{
							std::cout << "\nBuffer encheu! Hora de copiar para a GPU!";
							gpu::esvaziar_buffer(gpu_vars);
						}

						gpu::add_buffer(gpu_vars, &(tids_attval_dim_inq[0]), tids_attval_dim_inq.size( ));

						i_attval_inq++;
						it_attval++;

					}

					if (gpu_cheia || i_attval_inq == cardinalidade_dim_inquired)
					{

						if (gpu_cheia)
						{
							exit(0);
						}

						gpu::add_buffer(gpu_vars, &(inicio[0]), inicio.size( ));
						gpu::add_buffer(gpu_vars, &(fim[0]), fim.size());
						gpu::add_buffer(gpu_vars, &(deslocamento_celulas_inicio[0]), deslocamento_celulas_inicio.size());
						gpu::add_buffer(gpu_vars, &(deslocamento_celulas_fim[0]), deslocamento_celulas_fim.size());

						gpu::esvaziar_buffer(gpu_vars);

						int tamanho_bloco;
						int tamanho_grid;

						int total_threads = replicas * freq_tids_celulas_agregadas;;
						gpu::balancear_carga(gpu_vars, total_threads, &tamanho_bloco, &tamanho_grid);

						gpu::KernelVarsSubcube kernel_vars;

						kernel_vars.replicas = replicas;
						kernel_vars.cardinalidade = cardinalidade_dim_inquired;
						kernel_vars.total_celulas_derivadas = total_celulas_derivadas;
						kernel_vars.inicio_memoria = gpu_vars->inicio_memoria;
						kernel_vars.freq_atributos_celulas_agregadas = freq_tids_celulas_agregadas;
						kernel_vars.inicio_espaco_busca = kernel_vars.inicio_memoria + (kernel_vars.replicas * kernel_vars.freq_atributos_celulas_agregadas);

						kernel_vars.offset_inicio = kernel_vars.inicio_espaco_busca + total_tids_attval_inquired;
						kernel_vars.offset_fim = kernel_vars.offset_inicio + kernel_vars.cardinalidade;

						kernel_vars.deslocamento_celulas_derivadas_inicio = kernel_vars.offset_fim + kernel_vars.cardinalidade;
						kernel_vars.deslocamento_celulas_derivadas_fim = kernel_vars.deslocamento_celulas_derivadas_inicio + total_celulas_derivadas;

						kernel_vars.resposta = kernel_vars.deslocamento_celulas_derivadas_fim + total_celulas_derivadas;

						int total_numeros_resposta = kernel_vars.freq_atributos_celulas_agregadas * kernel_vars.cardinalidade;

						gpu::subcubo_com_replicacao_tamanho_variavel << <tamanho_grid, tamanho_bloco >> >(kernel_vars);
						cudaDeviceSynchronize( );

						// fazendo a compactação das células
						int* inicio_celula_agregada = kernel_vars.inicio_memoria + freq_tids_celulas_agregadas;
						int offset_celula_agregada = freq_tids_celulas_agregadas;

						gpu_vars->total_numeros_na_gpu = freq_tids_celulas_agregadas;

						int total_tids_celulas_resultado = 0;
						int total_valores_atributo_inquired = inicio.size( );

						int total_celulas_inq = inicio.size( );

						int i_cel_derivadas = 0;
						int num_celula_inq = 0;

						while (i_cel_derivadas < total_celulas_derivadas)
						{
							int i_cel_inq = 0;
							int index = 0;

							while (i_cel_inq < total_celulas_inq)
							{
								// celulas da dimensao inquired processada
								int tam_celula_resultante = deslocamento_celulas_fim[index] - deslocamento_celulas_inicio[index];
								int tamanho_array_compacto;

								gpu::gpu_stream_compaction(kernel_vars.resposta + total_tids_celulas_resultado, inicio_celula_agregada,
									tam_celula_resultante, &tamanho_array_compacto);

								total_tids_celulas_resultado += tam_celula_resultante;

								if (tamanho_array_compacto > 0)
								{
									cout << "\n\n\nCuboid processado: " << dimensao_inquired;
									cout << "\nInicio da celula agregada de base: " << deslocamento_celulas_inicio[index];
									cout << "\nFim da celula agregada de base: " << deslocamento_celulas_fim[index];

									int indice_cuboid = descobrir_cuboid_celula(&deslocamentos_cuboids_vector_pairs, deslocamento_celulas_fim[index]);
									cout << "\nIndice do cuboid: " << indice_cuboid;
									cout << "\nCuboid antigo: " << label_cuboid_dimensions(deslocamentos_cuboids_vector[indice_cuboid]);

									CubDim novo_cuboid = create_cuboid_dimensions(deslocamentos_cuboids_vector[indice_cuboid], dimensao_inquired);
									cout << "\nCuboid criado: " << label_cuboid_dimensions(novo_cuboid);

									if (deslocamentos_cuboids_set.find(novo_cuboid) == deslocamentos_cuboids_set.end( ))
									{
										deslocamentos_cuboids_set.insert(novo_cuboid);
										pair_deslocamento_first = freq_tids_celulas_agregadas;
										pair_deslocamento_last = -1;
										cout << "\n>>> Cuboid novo: " << label_cuboid_dimensions(novo_cuboid);
										cout << "\nO inicio do cuboid novo e: " << pair_deslocamento_first;
									}
									else
									{
										if (pair_deslocamento_last == -1)
										{

										}
									}


									inicio_celula_agregada = inicio_celula_agregada + tamanho_array_compacto;
									gpu_vars->total_numeros_na_gpu += tamanho_array_compacto;

									// --------
									//int fator = cardinalidades_inquired_ordenadas[i_card] - (num_celula_inq / total_celulas_derivadas);
									int fator_2 = num_celula_inq % total_celulas_derivadas;

									string valor = "." + to_string(fator_2);
									string label = celulas[fator_2] + ", " + to_string((*dimensoes_inq)[i_dim_inq]) + valor;
									celulas.push_back(label);
									// --------

									offset_celula_agregada += tamanho_array_compacto;
									freq_tids_celulas_agregadas += tamanho_array_compacto;

									deslocamento_celulas_inicio.push_back(offset_celula_agregada - tamanho_array_compacto);
									deslocamento_celulas_fim.push_back(offset_celula_agregada);

								}

								num_celula_inq++;
								i_cel_inq++;

								index++;

							}
					
							i_cel_derivadas++;

						}

						inicio.clear( );
						fim.clear( );

#if DEBUG
						printf("\n\n");
						gpu::exibir_memoria(gpu_vars->inicio_memoria, gpu_vars->total_numeros_na_gpu);

						std::cout << "\nDeslocamentos derivadas inicio: " << label_tuple_id_list(&deslocamento_celulas_inicio);
						std::cout << "\nDeslocamentos derivadas fim:     " << label_tuple_id_list(&deslocamento_celulas_fim);
						std::cout << "\nTotal celulas agregadas: " << deslocamento_celulas_inicio.size( );
#endif
						total_tids_attval_inquired = freq_atributos_dim_instanciadas + offset_celula_agregada;

					}

				}
				else
				{
					// Se não é ALL
					it_attval++;
				}

			}

			std::cout << "\n\n\n";

			total_celulas_derivadas = deslocamento_celulas_inicio.size();
#if DEBUG
			std::cout << "\nAcabei a dimensao " << (*dimensoes_inq)[i_dim_inq];
#endif

			// se acabou esta dimensao inquired, vamos para a proxima
			i_dim_inq++;

		}

	}

	gpu::copiar_gpu_para_cpu(gpu_vars->buffer.vetor, gpu_vars->inicio_memoria, gpu_vars->total_numeros_na_gpu * sizeof(int));
	gpu_vars->buffer.pos_buffer = gpu_vars->total_numeros_na_gpu;
	gpu_vars->total_numeros_na_gpu = 0;

	#if EXIBIR_RESULTADO
	for (int j = 0; j < celulas.size( ); j++)
	{
		int tam = deslocamento_celulas_fim[j] - deslocamento_celulas_inicio[j];
		std::cout << "\n(" + celulas[j] + ")  = COUNT: " << tam;
	}
	#endif

	std::cout << "\n\n\nTOTAL CELULAS: " << deslocamento_celulas_inicio.size();
	std::cout << "\nTotal threads: " << gpu_vars->num_cpu_cores;

	gpu::exibir_buffer(gpu_vars->buffer);

	exit(0);

	return NULL;
}

//esta funcao preenche o buffer com os atributos seguintes que não
// serão processdos nesse ciclo de processamento
int preencher_resto_buffer( )
{
	return 0;
}

#endif