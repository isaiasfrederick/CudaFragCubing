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
	vector<int> dimensoes_com_intervalo;

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

inline ResultadoInquired* processar_dimensoes_inquired(gpu::GpuVars* gpu_vars,
	Fragmentos* fragmentos, Query query, AttValsIndexados* dimensoes_instanciadas, vector<string>* celulas, vector<int>* medidas);

inline ResultadoInquired* merge_cuboid_inquired(ResultadoInquired* dimensoes_inquired);

inline AttValsIndexados* processar_dimensoes_agregadas(Fragmentos* fragmentos, Query query);

inline AttValsIndexadosPorCuboid *coletar_todas_dimensoes_inquired(Fragmentos* fragmentos, Query query);

inline void mensagem_executando_query(string query)
{
	std::cout << "\n    Executing query: " << query;
}

inline void mensagem_resultado_query(int total_linhas, int ms)
{
	std::cout << "\n\n    " << total_linhas << " rows returned. Query executed in " << ms << " ms.";
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

	vector<string> celulas;
	vector<int> medidas;
	int flag_gpu = 0;

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
			#if EXECUTAR_EM_CPU
			dimensoes_inquired = processar_dimensoes_inquired(fragmentos, query, dimensoes_instanciadas);
			dimensoes_inquired_merged = merge_cuboid_inquired(dimensoes_inquired);
			#else
			dimensoes_inquired = processar_dimensoes_inquired(gpu_vars, fragmentos, query, dimensoes_instanciadas, &celulas, &medidas);
			flag_gpu = 1;
			#endif
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
	if (query.dimensoes_inquired.size() == query.comando.size())
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

	// consulta a sub-cubo em GPU
	if (flag_gpu)
	{
		int total_linhas = 0;

		mensagem_executando_query(query.comando_str);

		#if EXIBIR_TUPLAS
		cout << "\n";
		register int total_celulas = celulas.size( );
		for (register int i = 0; i < total_celulas; i++)
			cout << "\n    (" << celulas[i] << "): " << medidas[i];
		#endif

		marcar_tempo_final(&contadores_tempo);

		total_linhas = celulas.size( );
		contabilizar_tempo_gasto(&contadores_tempo);
		mensagem_resultado_query(total_linhas, contadores_tempo.elapsed_time);

	}
	else
	{
		if (!dimensoes_agregadas->empty( ))
		{
			int total_linhas = 0;

			#if EXIBIR_TUPLAS
			CubDim dimensoes_cuboid = create_cuboid_dimensions(query.dimensoes_aggregated);
			AttVal attribute_value = dimensoes_agregadas->begin()->first;

			mensagem_executando_query(query.comando_str);
			string label = construir_attribute_label(query, attribute_value);
			std::cout << "\n    " + label << " : " << (*dimensoes_agregadas)[attribute_value].size();

			#if EXECUTAR_HAN
			(*detector_colisoes)[label] = (*dimensoes_agregadas)[attribute_value].size();
#endif
#endif

			total_linhas++;

			marcar_tempo_final(&contadores_tempo);
			contabilizar_tempo_gasto(&contadores_tempo);

			mensagem_resultado_query(total_linhas, contadores_tempo.elapsed_time);

		}
		// instanciadas
		else if (!query.dimensoes_instanciadas.empty( ) && query.dimensoes_inquired.empty( ))
		{
			int total_linhas = 0;

			AttValsIndexados::iterator it_attribute_value;
			CubDim dimensoes_cuboid = create_cuboid_dimensions(query.dimensoes_instanciadas);
			it_attribute_value = dimensoes_instanciadas->begin();

			mensagem_executando_query(query.comando_str);

			if (!dimensoes_instanciadas->empty())
			{
				total_linhas++;

				#if EXIBIR_TUPLAS
				vector<int>* tuple_ids = &((*dimensoes_instanciadas)[it_attribute_value->first]);
				string label = construir_attribute_label(query, it_attribute_value->first);
				vector<int> tmp = dimensoes_instanciadas->begin()->first.values;
				std::cout << "\n    " << label << " : " << tuple_ids->size();

				#if EXECUTAR_HAN
				(*detector_colisoes)[label] = tuple_ids->size( );
				#endif
				#endif

			}

			marcar_tempo_final(&contadores_tempo);
			contabilizar_tempo_gasto(&contadores_tempo);

			mensagem_resultado_query(total_linhas, contadores_tempo.elapsed_time);

		}
		// inquired
		else if (!dimensoes_inquired_merged->empty() && query.dimensoes_instanciadas.empty())
		{
			int total_linhas;
			ResultadoInquired::iterator it_dimensoes_cuboid;
			it_dimensoes_cuboid = dimensoes_inquired_merged->begin();

			AttValsIndexados::iterator it_attribute_value;
			it_attribute_value = ((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first]).begin();

			total_linhas = ((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first]).size();

			mensagem_executando_query(query.comando_str);

#if EXIBIR_TUPLAS
			while (it_attribute_value != ((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first]).end())
			{
				vector<int>* tuple_ids = &((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first][it_attribute_value->first]);

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

			mensagem_resultado_query(total_linhas, contadores_tempo.elapsed_time);

		}
		// Consulta
		else if (query.dimensoes_instanciadas.empty() && (query.dimensoes_inquired.empty() || (!query.dimensoes_inquired.empty() && query.dimensoes_aggregated.empty())))
		{
			std::cout << "\nInvalid command.\n";
		}
		// Inquired + Instanciadas
		else if (!dimensoes_inquired_merged->empty() && !query.dimensoes_instanciadas.empty())
		{
			int total_linhas = 0;

			mensagem_executando_query(query.comando_str);

#if EXECUTAR_EM_CPU
			ResultadoInquired::iterator it_dimensoes_cuboid;

			it_dimensoes_cuboid = dimensoes_inquired_merged->begin();

			AttValsIndexados::iterator it_attribute_value;

			it_attribute_value = ((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first]).begin();
			total_linhas = ((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first]).size();

#if EXIBIR_TUPLAS
			while (it_attribute_value != ((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first]).end())
			{

				vector<int>* tuple_ids = &((*dimensoes_inquired_merged)[it_dimensoes_cuboid->first][it_attribute_value->first]);
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
#endif

			mensagem_resultado_query(total_linhas, contadores_tempo.elapsed_time);

		}
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

	for (register int i = 0; i < parametros.size(); i++)
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
		else if (parametros[i].length() > 0)
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
	CubDim last_dimensoes_cuboid;

	resultado = new AttValsIndexados;

	for (int d = 0; d < query.dimensoes_instanciadas.size(); d++)
	{
		CubDim dimensoes_cuboid = create_cuboid_dimensions(query.dimensoes_instanciadas[d]);
		AttVal attribute_value = create_attribute_value(query.valores_instanciados[d]);

		if (d == 0)
		{
			if (fragmentos->cuboids[dimensoes_cuboid].values.find(attribute_value) != fragmentos->cuboids[dimensoes_cuboid].values.end())
			{
				last_attribute_value = attribute_value;
				last_dimensoes_cuboid = dimensoes_cuboid;
				(*resultado)[attribute_value] = fragmentos->cuboids[dimensoes_cuboid].values[attribute_value];
			}
			else
				return resultado;
		}
		else
		{
			if (fragmentos->cuboids[dimensoes_cuboid].values.find(attribute_value) != fragmentos->cuboids[dimensoes_cuboid].values.end())
			{
				AttVal attribute_value_all = create_attribute_value(ALL);
				vector<int>* tuple_ids_1 = &(fragmentos->cuboids[dimensoes_cuboid].values[attribute_value]);
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
						last_dimensoes_cuboid,
						dimensoes_cuboid
					);

					last_dimensoes_cuboid = concatenate_cuboid_dimensions(last_dimensoes_cuboid, dimensoes_cuboid);
					last_attribute_value = novo_attribute_value;

					resultado->clear();

					(*resultado)[novo_attribute_value] = retorno;

				}
				else
				{
					resultado->clear( );
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
	CubDim last_dimensoes_cuboid;

	int grid_size;
	int block_size;

	resultado = new AttValsIndexados;

	if (query.dimensoes_instanciadas.size( ) == 1)
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
		grid_size = (kernel_vars[s].total_elementos_vetores_menores / block_size) + 1;

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
		CubDim dimensoes_cuboid = create_cuboid_dimensions(query.dimensoes_inquired[d]);
		AttValsIndexados::iterator it_attribute_value;

		it_attribute_value = fragmentos->cuboids[dimensoes_cuboid].values.begin();

		while (it_attribute_value != fragmentos->cuboids[dimensoes_cuboid].values.end())
		{
			if (!query.dimensoes_instanciadas.empty())
			{
				AttVal attribute_value_all = create_attribute_value(ALL);
				vector<int>* tuple_id_list = &(fragmentos->cuboids[dimensoes_cuboid].values[it_attribute_value->first]);
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
					if (resultado->find(dimensoes_cuboid) == resultado->end())
						(*resultado)[dimensoes_cuboid] = AttValsIndexados();

					(*resultado)[dimensoes_cuboid][it_attribute_value->first] = retorno;
				}

			}
			else
			{
				if (resultado->find(dimensoes_cuboid) == resultado->end())
					(*resultado)[dimensoes_cuboid] = AttValsIndexados();

				(*resultado)[dimensoes_cuboid][it_attribute_value->first] = fragmentos->cuboids[dimensoes_cuboid].values[it_attribute_value->first];

			}

			it_attribute_value++;

		}

	}

	return resultado;

}

ResultadoInquired*
merge_cuboid_inquired(ResultadoInquired* dimensoes_inquired)
{
	map<CubDim, AttValsIndexados>::iterator it_dimensoes_cuboid = dimensoes_inquired->begin();
	ResultadoInquired* resultado;

	resultado = new AttValsIndexadosPorCuboid;

	// Iterando dimensoes de dimensions_inquired
	while (it_dimensoes_cuboid != dimensoes_inquired->end())
	{
		if (it_dimensoes_cuboid == dimensoes_inquired->begin())
		{
			(*resultado)[it_dimensoes_cuboid->first] = (*dimensoes_inquired)[it_dimensoes_cuboid->first];
		}
		else
		{
			AttValsIndexados::iterator iterator_attribute_value;
			iterator_attribute_value = (*dimensoes_inquired)[it_dimensoes_cuboid->first].begin();

			ResultadoInquired* resultado_tmp;
			resultado_tmp = new AttValsIndexadosPorCuboid;

			// iterando attribute values de dimensoes inquired;
			while (iterator_attribute_value != (*dimensoes_inquired)[it_dimensoes_cuboid->first].end())
			{
				ResultadoInquired::iterator cuboid_resultado;
				cuboid_resultado = resultado->begin();

				AttValsIndexados::iterator iterator_attribute_value_resultado;
				iterator_attribute_value_resultado = (*resultado)[cuboid_resultado->first].begin();

				vector<int>* tuple_ids_1 = &(*dimensoes_inquired)[it_dimensoes_cuboid->first][iterator_attribute_value->first];

				// iterando attribute values de resultado
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
						CubDim novo_dimensoes_cuboid = concatenate_cuboid_dimensions(cuboid_resultado->first, it_dimensoes_cuboid->first);
						AttVal novo_attribute_values = concatenate_attribute_value
							(
							iterator_attribute_value_resultado->first,
							iterator_attribute_value->first,
							cuboid_resultado->first,
							it_dimensoes_cuboid->first
							);

						(*resultado_tmp)[novo_dimensoes_cuboid].insert({ novo_attribute_values, resultado });

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

		it_dimensoes_cuboid++;
	}

	return resultado;
}

void limpar_query(Query* query)
{
	query->comando.clear();
	query->dimensoes_inquired.clear();
	query->dimensoes_instanciadas.clear();
	query->dimensoes_com_intervalo.clear();
	query->mapeamento_colunas.clear();
	query->valores_instanciados.clear();
	query->dimensoes_aggregated.clear();
}

// indexar tuplas no detector de colisoes
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

	while (i < deslocamentos_cuboids_vector_pairs->size())
	{
		if (fim == (*deslocamentos_cuboids_vector_pairs)[i].second)
			return i;
		else if (fim >(*deslocamentos_cuboids_vector_pairs)[i].second)
			return i + 1;

		i++;
	}

	return i;
}

inline ResultadoInquired* processar_dimensoes_inquired(gpu::GpuVars* gpu_vars,
	Fragmentos* fragmentos, Query query, AttValsIndexados* dimensoes_instanciadas, vector<string>* celulas, vector<int>* medidas)
{
	vector<int> vector_dimensoes_instanciadas = query.dimensoes_instanciadas;
	vector<int> cardinalidades_inquired_ordenadas;
	map < int, vector<int>> map_inq_dim_por_card;
	string label_dimensoes_instanciadas = "INST";
	vector<int> tids_instanciadas;

	for (register int i = 0; i < query.dimensoes_inquired.size(); i++)
	{
		int inq_dim = query.dimensoes_inquired[i];
		int i_card_inq_dim = fragmentos->cardinalidades_por_coluna[inq_dim - 1];

		if (map_inq_dim_por_card.find(i_card_inq_dim) == map_inq_dim_por_card.end())
		{
			vector<int> vetor;
			map_inq_dim_por_card[i_card_inq_dim] = vetor;
			cardinalidades_inquired_ordenadas.push_back(i_card_inq_dim);
		}

		map_inq_dim_por_card[i_card_inq_dim].push_back(inq_dim);

	}

	std::sort(cardinalidades_inquired_ordenadas.begin(), cardinalidades_inquired_ordenadas.end());
	std::sort(vector_dimensoes_instanciadas.begin(), vector_dimensoes_instanciadas.end());

	tids_instanciadas = (*dimensoes_instanciadas)[dimensoes_instanciadas->begin( )->first];

	// inicio espaco busca
	vector<int> inicio_espaco_busca;

	// fim espaco busca
	vector<int> fim_espaco_busca;

	int freq_atributos_dim_instanciadas = tids_instanciadas.size( );
	map<AttVal, vector<int>>::iterator it_attval_ref;
	AttVal attval_all = create_attribute_value(ALL);
	int freq_valor_atributo_inq = 0;

	// total tids adicionadas a gpu
	int total_tids_attval_inquired = 0;
	int replicas;

	CubDim cuboid_dimensoes_instanciadas = create_cuboid_dimensions(vector_dimensoes_instanciadas);

	celulas->push_back(label_dimensoes_instanciadas);
	medidas->push_back(freq_atributos_dim_instanciadas);

	vector<int> deslocamento_celulas_inicio;
	vector<int> deslocamento_celulas_fim;

	int dimensao_inquired;

	// iterando cardinalidades inquired já ordenadas
	for (int i_card = 0; i_card < cardinalidades_inquired_ordenadas.size( ); i_card++)
	{
		int cardinalidade_dim_inquired = cardinalidades_inquired_ordenadas[i_card];

		// obtendo dimensões inquired com a cardinalidade em questão
		vector<int>* dimensoes_inq = &(map_inq_dim_por_card[cardinalidade_dim_inquired]);

		int i_dim_inq = 0;

		CubDim cuboid_dimensao_processada = create_cuboid_dimensions(cuboid_dimensoes_instanciadas, (*dimensoes_inq)[i_dim_inq]);

		// processando primeira dimensão inquired
		if (!i_card && !i_dim_inq)
		{
			dimensao_inquired = (*dimensoes_inq)[i_dim_inq];
			CubDim cub_dim = create_cuboid_dimensions(dimensao_inquired);
			Cuboid cuboid = *get_cuboid(fragmentos, cub_dim);
			map<AttVal, vector<int>>::iterator it_attval_ref = get_attval_iterator(&cuboid);

			// se trocou de dimensao
			total_tids_attval_inquired = 0;

			AttVal attval;
			map<AttVal, vector<int>>::iterator it_attval = get_attval_iterator(&cuboid);

			// calculando o maximo de replicas possiveis
			int max_replicas_gpu = gpu::calcular_replicacao
			(
				gpu_vars, &cuboid,
				&it_attval_ref,
				freq_atributos_dim_instanciadas,
				TAMANHO_FIXO,
				freq_atributos_dim_instanciadas
			);

			// fazendo replicas
			replicas = gpu::replicar_tids(gpu_vars, max_replicas_gpu, freq_atributos_dim_instanciadas);

			// contador valor de atributo inquired
			int cont_attval_inq = 0;

			// iterando valores de atributo
			while (it_attval != cuboid.values.end())
			{
				// TIDs do valor de atributo atual
				vector<int> tids_attval_dim_inq = get_tids(&cuboid, it_attval->first);

				// se nao e all
				if (attribute_values_equals(it_attval->first, attval_all) == false)
				{
					freq_valor_atributo_inq = tids_attval_dim_inq.size();

					bool gpu_cheia = (it_attval == it_attval_ref);

					if (!gpu_cheia && cont_attval_inq < cardinalidade_dim_inquired)
					{
						// adicionando mais um atributo ao buffer
						total_tids_attval_inquired += freq_valor_atributo_inq;

						inicio_espaco_busca.push_back(total_tids_attval_inquired - freq_valor_atributo_inq);
						fim_espaco_busca.push_back(total_tids_attval_inquired);

						// vendo se o buffer irá saturar
						if (gpu::buffer_cheio(&gpu_vars->buffer, freq_valor_atributo_inq) == true)
						{
							// flush buffer
							gpu::esvaziar_buffer(gpu_vars);
						}

						if (gpu::add_buffer(gpu_vars, &(tids_attval_dim_inq[0]), tids_attval_dim_inq.size( )) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						// incrementando contador de valor de atributos percorridos
						cont_attval_inq++;

						// ir pro proximo valor de atributo
						it_attval++;

					}

					// se a gpu encheu ou percorrel dos valores de atributo para a dimensao em questao
					if (gpu_cheia || cont_attval_inq == cardinalidade_dim_inquired)
					{
						if (gpu_cheia)
						{
							// ainda nao implementa esta parte
							printf("\nA GPU ENCHEU! ABORTANDO...");
							exit(0);
						}

						if (gpu::add_buffer(gpu_vars, &(inicio_espaco_busca[0]), inicio_espaco_busca.size( )) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						if (gpu::add_buffer(gpu_vars, &(fim_espaco_busca[0]), fim_espaco_busca.size( )) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						if (gpu_vars->buffer.pos_buffer)
							gpu::esvaziar_buffer(gpu_vars);

						int tamanho_bloco;
						int tamanho_grid;

						int total_threads = replicas * freq_atributos_dim_instanciadas;
						gpu::balancear_carga(gpu_vars, total_threads, &tamanho_bloco, &tamanho_grid);

						inicio_espaco_busca.clear( );
						fim_espaco_busca.clear( );

						gpu::KernelVarsSubcube kernel_vars;

						kernel_vars.replicas = replicas;
						kernel_vars.cardinalidade = cardinalidade_dim_inquired;
						kernel_vars.inicio_memoria = gpu_vars->inicio_memoria;
						kernel_vars.freq_atributos_dim_instanciadas = freq_atributos_dim_instanciadas;
						kernel_vars.inicio_espaco_busca = kernel_vars.inicio_memoria + (kernel_vars.replicas * kernel_vars.freq_atributos_dim_instanciadas);
						kernel_vars.deslocamento_celulas_inquired_inicio = kernel_vars.inicio_espaco_busca + total_tids_attval_inquired;
						kernel_vars.deslocamento_celulas_inquired_fim = kernel_vars.deslocamento_celulas_inquired_inicio + kernel_vars.cardinalidade;
						kernel_vars.resposta = kernel_vars.deslocamento_celulas_inquired_fim + kernel_vars.cardinalidade;

						gpu::subcubo_com_replicacao << <tamanho_grid, tamanho_bloco >> >(kernel_vars);
						cudaDeviceSynchronize( );

						// fazendo a compactação das células
						int contador_celulas = 0;
						int* inicio_celula_derivada = kernel_vars.inicio_memoria + freq_atributos_dim_instanciadas;
						int offset_celulas_derivadas = freq_atributos_dim_instanciadas;

						// adicionando offsets da celula base composta so de valores instanciados
						deslocamento_celulas_inicio.push_back(offset_celulas_derivadas - freq_atributos_dim_instanciadas);
						deslocamento_celulas_fim.push_back(offset_celulas_derivadas);

						gpu_vars->total_numeros_na_gpu = freq_atributos_dim_instanciadas;
						
						#if COMPACTACAO
						// total celulas resultantes
						int total_celulas_novas = cont_attval_inq;
						int total_celulas_derivadas = 1;

						// utilizando minha propria compactacao
						// tamanho celulas resultantes
						int* tam_celulas_gpu;
						int* tam_celulas_cpu;

						cudaMalloc(&tam_celulas_gpu, total_celulas_novas * sizeof(int));

						tam_celulas_cpu = remover_esparsidade
						(
							gpu_vars,
							kernel_vars.resposta,
							total_celulas_novas,
							total_celulas_derivadas,
							tam_celulas_gpu,
							freq_atributos_dim_instanciadas,
							inicio_celula_derivada
						);

						cudaFree(tam_celulas_gpu);
						tam_celulas_gpu = NULL;

						while (contador_celulas < cont_attval_inq)
						{
							int tamanho_resultado_compacto = tam_celulas_cpu[contador_celulas];

							if (tamanho_resultado_compacto)
							{
								inicio_celula_derivada = inicio_celula_derivada + tamanho_resultado_compacto;
								gpu_vars->total_numeros_na_gpu += tamanho_resultado_compacto;

								// ----
								register string valor = "." + to_string(cardinalidades_inquired_ordenadas[i_card] - contador_celulas);
								register string label = label_dimensoes_instanciadas + ", " + to_string((*dimensoes_inq)[i_dim_inq]) + valor;

								celulas->push_back(label);
								medidas->push_back(tamanho_resultado_compacto);
								// ----

								offset_celulas_derivadas += tamanho_resultado_compacto;

								// adicionando mais uma celula do sub-cubo
								deslocamento_celulas_inicio.push_back(offset_celulas_derivadas - tamanho_resultado_compacto);
								deslocamento_celulas_fim.push_back(offset_celulas_derivadas);

							}

							contador_celulas++;

						}

						free(tam_celulas_cpu);
						tam_celulas_cpu = NULL;

						#else
						while (contador_celulas < cont_attval_inq)
						{

							int tamanho_resultado_compacto;
							gpu::gpu_stream_compaction
							(
								kernel_vars.resposta + (contador_celulas * freq_atributos_dim_instanciadas),
								inicio_celula_derivada,
								freq_atributos_dim_instanciadas,
								&tamanho_resultado_compacto
							);

							if (tamanho_resultado_compacto)
							{
								inicio_celula_derivada = inicio_celula_derivada + tamanho_resultado_compacto;
								gpu_vars->total_numeros_na_gpu += tamanho_resultado_compacto;

								// ----
								register string valor = "." + to_string(cardinalidades_inquired_ordenadas[i_card] - contador_celulas);
								register string label = label_dimensoes_instanciadas + ", " + to_string((*dimensoes_inq)[i_dim_inq]) + valor;

								celulas->push_back(label);
								medidas->push_back(tamanho_resultado_compacto);
								// ----

								offset_celulas_derivadas += tamanho_resultado_compacto;

								// adicionando mais uma celula do sub-cubo
								deslocamento_celulas_inicio.push_back(offset_celulas_derivadas - tamanho_resultado_compacto);
								deslocamento_celulas_fim.push_back(offset_celulas_derivadas);

							}

							contador_celulas++;

						}
						#endif

					}

				}
				else
				{
					// se é all, pule
					it_attval++;
				}

			}

			// se acabou esta dimensao inquired, vamos para a proxima
			i_dim_inq++;

		}

		// PROCESSANDO DIMENSOES INQUIRED DIFERENTES DA PRIMEIRA
		// somatorio do count de tids de todas as celulas obtidas na iteracao anterior
		int freq_tids_celulas_derivadas = gpu_vars->total_numeros_na_gpu;

		// total de celulas obtidas da outra agregacacao
		int total_celulas_derivadas = deslocamento_celulas_inicio.size( );

		while (i_dim_inq < dimensoes_inq->size( ))
		{
			dimensao_inquired = (*dimensoes_inq)[i_dim_inq];
			CubDim cub_dim = create_cuboid_dimensions(dimensao_inquired);
			Cuboid cuboid = *get_cuboid(fragmentos, cub_dim);
			map<AttVal, vector<int>>::iterator it_attval_ref = get_attval_iterator(&cuboid);

			// se trocou de dimensao
			total_tids_attval_inquired = 0;

			AttVal attval;
			map<AttVal, vector<int>>::iterator it_attval = get_attval_iterator(&cuboid);

			int max_replicas_gpu = gpu::calcular_replicacao(gpu_vars, &cuboid, &it_attval_ref, freq_tids_celulas_derivadas,
				TAMANHO_VARIAVEL, total_celulas_derivadas);

			replicas = gpu::replicar_tids(gpu_vars, max_replicas_gpu, freq_tids_celulas_derivadas);

			int i_attval_inq = 0;

			// iterando atributos
			while (it_attval != cuboid.values.end( ))
			{
				vector<int> tids_attval_dim_inq = get_tids(&cuboid, it_attval->first);

				if (attribute_values_equals(it_attval->first, attval_all) == false)
				{
					freq_valor_atributo_inq = tids_attval_dim_inq.size();

					bool gpu_cheia = (it_attval == it_attval_ref);

					if (!gpu_cheia && i_attval_inq < cardinalidade_dim_inquired)
					{
						// adicionando mais um atributo ao Buffer
						total_tids_attval_inquired += freq_valor_atributo_inq;

						// deslocamentos do espaco de busca dos valores de atributo
						inicio_espaco_busca.push_back(total_tids_attval_inquired - freq_valor_atributo_inq);
						fim_espaco_busca.push_back(total_tids_attval_inquired);

						gpu::add_buffer(gpu_vars, &(tids_attval_dim_inq[0]), tids_attval_dim_inq.size( ));

						i_attval_inq++;
						it_attval++;

					}

					if (gpu_cheia || i_attval_inq == cardinalidade_dim_inquired)
					{

						if (gpu_cheia)
						{
							cout << "\nA GPU ENCHEU!";
							exit(0);
						}

						// deslocamento celulas inquired
						if (gpu::add_buffer(gpu_vars, &(inicio_espaco_busca[0]), inicio_espaco_busca.size()) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						if (gpu::add_buffer(gpu_vars, &(fim_espaco_busca[0]), fim_espaco_busca.size( )) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						// deslocamento celulas agregadas
						if (gpu::add_buffer(gpu_vars, &(deslocamento_celulas_inicio[0]), deslocamento_celulas_inicio.size( )) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						if (gpu::add_buffer(gpu_vars, &(deslocamento_celulas_fim[0]), deslocamento_celulas_fim.size( )) != SUCESSO_COPIA)
						{
							cout << "\nA GPU encheu! Abortando...";
							exit(0);
						}

						// se buffer possui elemento
						if (gpu_vars->buffer.pos_buffer)
							gpu::esvaziar_buffer(gpu_vars);

						int tamanho_bloco;
						int tamanho_grid;

						int total_threads = replicas * freq_tids_celulas_derivadas;;
						gpu::balancear_carga(gpu_vars, total_threads, &tamanho_bloco, &tamanho_grid);

						gpu::KernelVarsSubcube kernel_vars;

						kernel_vars.replicas = replicas;
						kernel_vars.cardinalidade = cardinalidade_dim_inquired;
						kernel_vars.total_celulas_derivadas = total_celulas_derivadas;
						kernel_vars.inicio_memoria = gpu_vars->inicio_memoria;
						kernel_vars.freq_atributos_celulas_derivadas = freq_tids_celulas_derivadas; // total de TIDs com todas celulas derivadas
						kernel_vars.inicio_espaco_busca = kernel_vars.inicio_memoria + (kernel_vars.replicas * kernel_vars.freq_atributos_celulas_derivadas);

						kernel_vars.deslocamento_celulas_inquired_inicio = kernel_vars.inicio_espaco_busca + total_tids_attval_inquired;
						kernel_vars.deslocamento_celulas_inquired_fim = kernel_vars.deslocamento_celulas_inquired_inicio + kernel_vars.cardinalidade;

						kernel_vars.deslocamento_celulas_derivadas_inicio = kernel_vars.deslocamento_celulas_inquired_fim + kernel_vars.cardinalidade;
						kernel_vars.deslocamento_celulas_derivadas_fim = kernel_vars.deslocamento_celulas_derivadas_inicio + total_celulas_derivadas;

						kernel_vars.resposta = kernel_vars.deslocamento_celulas_derivadas_fim + total_celulas_derivadas;

						gpu::subcubo_com_replicacao_tamanho_variavel << <tamanho_grid, tamanho_bloco >> >(kernel_vars);
						cudaDeviceSynchronize( );

						// fazendo a compactação das células
						int* inicio_celula_derivada = kernel_vars.inicio_memoria + freq_tids_celulas_derivadas;
						int offset_celulas_derivadas = freq_tids_celulas_derivadas;

						gpu_vars->total_numeros_na_gpu = freq_tids_celulas_derivadas;

						int total_valores_atributo_inquired = inicio_espaco_busca.size( );
						int total_tids_celulas_resultado = 0;
						int i_celula_derivada = 0;

						int total_celulas_inq = inicio_espaco_busca.size( );

						int i_cel_derivada = 0;
						int num_celula_inq = 0;

						#if COMPACTACAO
						// total celulas resultantes
						int total_celulas_novas = total_celulas_derivadas * total_celulas_inq;

						// utilizando minha propria compactacao
						// tamanho celulas resultantes
						int* tam_celulas_gpu;
						int* tam_celulas_cpu;

						cudaMalloc(&tam_celulas_gpu, total_celulas_novas * sizeof(int));

						tam_celulas_cpu = remover_esparsidade
						(
							gpu_vars,
							kernel_vars.resposta,
							total_celulas_novas,
							total_celulas_derivadas,
							tam_celulas_gpu,
							kernel_vars.deslocamento_celulas_derivadas_inicio,
							kernel_vars.deslocamento_celulas_derivadas_fim,
							inicio_celula_derivada
						);

						cudaFree(tam_celulas_gpu);
						tam_celulas_gpu = NULL;

						while (i_cel_derivada < total_celulas_derivadas)
						{
							int i_cel_inq = 0;

							while (i_cel_inq < total_celulas_inq)
							{
								// celulas da dimensao inquired processada
								int tamanho_resultado_compacto = tam_celulas_cpu[num_celula_inq];

								if (tamanho_resultado_compacto > 0)
								{
									inicio_celula_derivada = inicio_celula_derivada + tamanho_resultado_compacto;
									gpu_vars->total_numeros_na_gpu += tamanho_resultado_compacto;

									int cel_agre = num_celula_inq % total_celulas_derivadas;
									int cel_inq = num_celula_inq / total_celulas_derivadas;

									register string valor = "." + to_string(cardinalidades_inquired_ordenadas[i_card] - cel_inq);
									register string label = (*celulas)[cel_agre] + ", " + to_string((*dimensoes_inq)[i_dim_inq]) + valor;

									celulas->push_back(label);
									medidas->push_back(tamanho_resultado_compacto);
									// ---

									offset_celulas_derivadas += tamanho_resultado_compacto;
									freq_tids_celulas_derivadas += tamanho_resultado_compacto;

									deslocamento_celulas_inicio.push_back(offset_celulas_derivadas - tamanho_resultado_compacto);
									deslocamento_celulas_fim.push_back(offset_celulas_derivadas);

								}

								num_celula_inq++;
								i_cel_inq++;

								if (i_celula_derivada + 1 == total_celulas_derivadas)
									i_celula_derivada = 0;
								else
									i_celula_derivada++;

							}
							
							i_cel_derivada++;

						}

						inicio_espaco_busca.clear( );
						fim_espaco_busca.clear( );

						total_tids_attval_inquired = freq_atributos_dim_instanciadas + offset_celulas_derivadas;
						#else

						// usando inumeras operacoes de compactacao
						i_cel_derivada = 0;
						num_celula_inq = 0;

						while (i_cel_derivada < total_celulas_derivadas)
						{
							int i_cel_inq = 0;

							while (i_cel_inq < total_celulas_inq)
							{
								// celulas da dimensao inquired processada
								int tam_celula_resultante = deslocamento_celulas_fim[i_celula_derivada] - deslocamento_celulas_inicio[i_celula_derivada];
								int tamanho_resultado_compacto;

								gpu::gpu_stream_compaction(kernel_vars.resposta + total_tids_celulas_resultado, inicio_celula_derivada,
									tam_celula_resultante, &tamanho_resultado_compacto);

								total_tids_celulas_resultado += tam_celula_resultante;

								if (tamanho_resultado_compacto > 0)
								{
									inicio_celula_derivada = inicio_celula_derivada + tamanho_resultado_compacto;
									gpu_vars->total_numeros_na_gpu += tamanho_resultado_compacto;

									int cel_agre = num_celula_inq % total_celulas_derivadas;
									int cel_inq = num_celula_inq / total_celulas_derivadas;

									register string valor = "." + to_string(cardinalidades_inquired_ordenadas[i_card] - cel_inq);
									register string label = (*celulas)[cel_agre] + ", " + to_string((*dimensoes_inq)[i_dim_inq]) + valor;

									celulas->push_back(label);
									medidas->push_back(tamanho_resultado_compacto);
									// ---

									offset_celulas_derivadas += tamanho_resultado_compacto;
									freq_tids_celulas_derivadas += tamanho_resultado_compacto;

									deslocamento_celulas_inicio.push_back(offset_celulas_derivadas - tamanho_resultado_compacto);
									deslocamento_celulas_fim.push_back(offset_celulas_derivadas);

								}

								num_celula_inq++;
								i_cel_inq++;

								if (i_celula_derivada + 1 == total_celulas_derivadas)
									i_celula_derivada = 0;
								else
									i_celula_derivada++;

							}

							i_cel_derivada++;

						}

						inicio_espaco_busca.clear();
						fim_espaco_busca.clear();

						total_tids_attval_inquired = freq_atributos_dim_instanciadas + offset_celulas_derivadas;
						#endif

					}

				}
				else
				{
					// se é all, pule
					it_attval++;
				}

			}

			total_celulas_derivadas = deslocamento_celulas_inicio.size();

			// se acabou esta dimensao inquired, vamos para a proxima
			i_dim_inq++;

		}

	}

	gpu::copiar_gpu_para_cpu(gpu_vars->buffer.vetor, gpu_vars->inicio_memoria, gpu_vars->total_numeros_na_gpu * sizeof(int));
	gpu_vars->buffer.pos_buffer = gpu_vars->total_numeros_na_gpu;
	gpu_vars->total_numeros_na_gpu = 0;

	return NULL;

}

#endif