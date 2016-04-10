#ifndef GPU_H
#define GPU_H

#include <iostream>
#include <cstdlib>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <math.h>

#include "cuda_runtime.h"
#include "constantes.h"
#include "device_launch_parameters.h"

using namespace std;

namespace gpu
{
	typedef struct Buffer
	{
		int max_tam_buffer;
		int max_barramento;
		int pos_buffer;
		int* vetor;
	} Buffer;

	typedef struct KernelVars
	{
		int total_intersecoes;

		// tabela de offsets
		int* inicio_vetores_menores;
		int* inicio_vetores_maiores;

		// contadores
		int total_elementos_vetores_menores;
		int total_elementos_vetores_maiores;

		int* inicio_espaco_busca;
		int* fim_espaco_busca;

		int* deslocamento_threads;
		int* resultado;

		int replicas;
		int cardinalidade_dim_inq;

		size_t max_bytes_alocados;
		size_t bytes_utilizados;

		cudaStream_t stream;

	} KernelVars;

	typedef struct KernelVarsSubcube
	{
		int replicas;
		int* inicio_memoria;
		int freq_atributos_dim_instanciadas;
		int freq_atributos_celulas_derivadas;
		int cardinalidade;

		int* deslocamento_celulas_derivadas_inicio;
		int* deslocamento_celulas_derivadas_fim;
		int total_celulas_dimensao_inquired;
		int* tamanho_novas_celulas;
		int total_celulas_derivadas;

		int* inicio_espaco_busca;
		int* deslocamento_celulas_inquired_inicio;
		int* deslocamento_celulas_inquired_fim;
		int* resposta;

	} KernelVarsSubcube;

	typedef struct GpuVars
	{
		int* inicio_memoria;
		long total_bytes_alocados;
		long total_bytes_ocupados;
		int total_numeros_na_gpu;
		int max_kernels_simultaneos;
		int total_kernels_simultaneos;
		int block_size;
		int max_blocks_por_multiprocessador;
		int num_cpu_cores;

		gpu::Buffer buffer;
		gpu::KernelVars* kernel_vars;

		cudaDeviceProp prop;

	} GpuVars;

	typedef struct SubcubeCicleVars
	{
		int total_threads;
		int total_atributos_copiados;
		int total_tids_copiadas;
	} SubcubeCicleVars;

	bool gpu_cheia(GpuVars* gpu_vars, int novos_bytes);
	int esvaziar_buffer(GpuVars* gpu_vars);

	void iniciar_subcube_cicle_vars(SubcubeCicleVars* sccv)
	{
		sccv->total_atributos_copiados = 0;
		sccv->total_tids_copiadas = 0;
		sccv->total_threads = 0;
	}

	struct is_not_invalid
	{
		__host__ __device__ bool operator()(const int x)
		{
			return (x != INVALID_NUMBER);
		}
	};

	GpuVars*  criar_gpu_vars( )
	{
		GpuVars* gpu_vars = (GpuVars*)malloc(sizeof(GpuVars));
		cudaGetDeviceProperties(&gpu_vars->prop, 0);

		return gpu_vars;
	}

	void desalocar_gpu(gpu::GpuVars* gpu)
	{
		cudaFree(&gpu->inicio_memoria);
		cudaFree(gpu->kernel_vars);

		gpu->kernel_vars = NULL;

		gpu->total_bytes_alocados = 0;
		gpu->total_bytes_ocupados = 0;
		gpu->inicio_memoria = NULL;
	}

	inline void balancear_carga(gpu::GpuVars* gpu, int total_threads, int* tamanho_bloco, int* tamanho_grid)
	{
		int divisao = total_threads / gpu->prop.warpSize;
		vector<int> m;

		*tamanho_bloco = gpu->prop.warpSize;

		int total_warps = total_threads % gpu->prop.warpSize ? divisao + 1 : divisao;

		if (total_warps > gpu->prop.multiProcessorCount)
		{
			int mi = 0;

			for (register int i = 0; i < gpu->prop.multiProcessorCount; i++)
				m.push_back(1);

			*tamanho_grid = gpu->prop.multiProcessorCount;
			total_warps -= gpu->prop.multiProcessorCount;

			while (total_warps)
			{
				if (((m[mi] + 1) * gpu->prop.warpSize) > gpu->prop.maxThreadsPerBlock)
				{
					int div_tmp = total_threads / gpu->prop.maxThreadsPerBlock;
					*tamanho_grid = total_threads % gpu->prop.maxThreadsPerBlock ? div_tmp + 1 : div_tmp;
					*tamanho_bloco = gpu->prop.maxThreadsPerBlock;
					return;
				}

				m[mi] = 1 + m[mi];

				if (mi + 1 == gpu->prop.multiProcessorCount)
				{
					mi = 0;

					if (*tamanho_bloco < gpu->prop.maxThreadsPerBlock)
						*tamanho_bloco = *tamanho_bloco + gpu->prop.warpSize;
				}
				else if (total_warps == 1 && mi < gpu->prop.multiProcessorCount - 1)
				{
					*tamanho_bloco = *tamanho_bloco + gpu->prop.warpSize;
				}

				mi++;
				total_warps--;

			}

		}
		else
			*tamanho_grid = total_warps;

	}

	void alocar_gpu(gpu::GpuVars* gpu)
	{
		size_t total, livre;

		cudaMemGetInfo(&livre, &total);

		cudaGetDeviceProperties(&gpu->prop, 0);

		gpu->total_bytes_alocados = (size_t)(livre * PER_CENT_MEM);
		gpu->total_bytes_ocupados = 0;
		gpu->block_size = DEFAULT_BLOCK_SIZE;
		gpu->max_blocks_por_multiprocessador = DEFAULT_BLOCKS_PER_MULTIPROCESSOR;
		gpu->buffer.max_barramento = MAX_BARRAMENTO;

		gpu->kernel_vars = (gpu::KernelVars*) malloc(sizeof(gpu::KernelVars) * gpu->max_kernels_simultaneos);

#if PINNED
		cudaHostAlloc(&gpu->inicio_memoria, gpu->total_bytes_alocados, cudaHostAllocDefault);
#else
		cudaMalloc(&gpu->inicio_memoria, gpu->total_bytes_alocados);
#endif
	}

	void inline particionar_memoria(gpu::GpuVars* gpu)
	{
		for (int s = 0; s < gpu->total_kernels_simultaneos; s++)
		{
			gpu->kernel_vars[s].total_intersecoes = 0;
			gpu->kernel_vars[s].total_elementos_vetores_maiores = 0;
			gpu->kernel_vars[s].total_elementos_vetores_menores = 0;
			gpu->kernel_vars[s].max_bytes_alocados = (gpu->total_bytes_alocados / gpu->total_kernels_simultaneos);
			gpu->kernel_vars[s].bytes_utilizados = 0;
			gpu->kernel_vars[s].inicio_vetores_menores = gpu->inicio_memoria + ((gpu->kernel_vars[s].max_bytes_alocados / sizeof(int)) * s);

			cudaStreamCreate(&(gpu->kernel_vars[s].stream));

		}

	}

	inline void configurar_kernel(int threads_necessarias, int* grid_size, int* block_size)
	{
		int block_size_tmp = 64;

		*block_size = block_size_tmp;
		*grid_size = threads_necessarias / *block_size + 1;
	}

	inline bool copiar_cpu_para_gpu(int* dest, int* source, int bytes)
	{
		return cudaMemcpy(dest, source, bytes, cudaMemcpyHostToDevice) == cudaSuccess;
	}

	inline bool copiar_gpu_para_cpu(int* dest, int* source, int bytes)
	{
		return cudaMemcpy(dest, source, bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
	}

	__device__ void contains(int* vetor, int tam_vetor, int key, int* retorno)
	{
		int e = 0;
		int d = tam_vetor - 1;
		int pivo = (e + d) / 2;

		while (e <= d)
		{
			if (key == vetor[pivo]) {
				*retorno = key;
				return;
			}
			else if (key > vetor[pivo]) {
				e = pivo + 1;
			}
			else if (key < vetor[pivo]) {
				d = pivo - 1;
			}

			pivo = (e + d) / 2;

		}

		*retorno = INVALID_NUMBER;
	}

	__global__ void subcubo_com_replicacao(KernelVarsSubcube kernel_vars)
	{
		int thread_id = (blockDim.x * blockIdx.x) + threadIdx.x;

		if (thread_id < kernel_vars.replicas * kernel_vars.freq_atributos_dim_instanciadas)
		{
			int key = kernel_vars.inicio_memoria[thread_id];
			int retorno;

			// descobrir tamanho do chunk
			int atributos_por_replica = (kernel_vars.cardinalidade / kernel_vars.replicas) + (kernel_vars.cardinalidade % kernel_vars.replicas);
			int tidlist_atr_i = (thread_id / kernel_vars.freq_atributos_dim_instanciadas) * atributos_por_replica;

			int limite_tidlist;

			limite_tidlist = tidlist_atr_i + atributos_por_replica;

			if ((thread_id / kernel_vars.freq_atributos_dim_instanciadas) + 1 == kernel_vars.replicas)
				if (kernel_vars.cardinalidade % atributos_por_replica)
					limite_tidlist = tidlist_atr_i + (kernel_vars.cardinalidade % atributos_por_replica);

			int index;

			while (tidlist_atr_i < limite_tidlist)
			{
				index = (thread_id % kernel_vars.freq_atributos_dim_instanciadas) + (tidlist_atr_i * kernel_vars.freq_atributos_dim_instanciadas);

				int* vetor = kernel_vars.inicio_espaco_busca + kernel_vars.deslocamento_celulas_inquired_inicio[tidlist_atr_i];
				int tamanho = kernel_vars.deslocamento_celulas_inquired_fim[tidlist_atr_i] - kernel_vars.deslocamento_celulas_inquired_inicio[tidlist_atr_i];

				contains(vetor, tamanho, key, &retorno);
				kernel_vars.resposta[index] = retorno;

				tidlist_atr_i++;
			}

		}

		__syncthreads();

	}

	__global__ void subcubo_com_replicacao_tamanho_variavel(KernelVarsSubcube kernel_vars)
	{
		int thread_id = (blockDim.x * blockIdx.x) + threadIdx.x;

		if (thread_id < kernel_vars.replicas * kernel_vars.freq_atributos_celulas_derivadas)
		{
			int num_replica = thread_id / kernel_vars.freq_atributos_celulas_derivadas;
			int num_cel = 0;
			int attval_min = (kernel_vars.cardinalidade / kernel_vars.replicas) * num_replica;
			int attval_max = attval_min + (kernel_vars.cardinalidade / kernel_vars.replicas);

			while (num_cel < kernel_vars.total_celulas_derivadas && (thread_id % kernel_vars.freq_atributos_celulas_derivadas) >= kernel_vars.deslocamento_celulas_derivadas_fim[num_cel])
				num_cel++;

			if (num_replica + 1 == kernel_vars.replicas)
				attval_max = attval_min + ((kernel_vars.cardinalidade % kernel_vars.replicas) + (kernel_vars.cardinalidade / kernel_vars.replicas));

			while (attval_min < attval_max)
			{
				int* inicio_espaco_busca = kernel_vars.inicio_espaco_busca + kernel_vars.deslocamento_celulas_inquired_inicio[attval_min];
				int tamanho_espaco_busca = kernel_vars.deslocamento_celulas_inquired_fim[attval_min] - kernel_vars.deslocamento_celulas_inquired_inicio[attval_min];
				int pos_resposta = (attval_min * kernel_vars.freq_atributos_celulas_derivadas) + (thread_id % (kernel_vars.freq_atributos_celulas_derivadas));

				contains(inicio_espaco_busca, tamanho_espaco_busca, kernel_vars.inicio_memoria[thread_id], kernel_vars.resposta + pos_resposta);

				attval_min++;

			}

		}

		__syncthreads();

	}

	__global__ void intersection_with_chunks(KernelVars vars)
	{
		int threadId = (blockDim.x * blockIdx.x) + threadIdx.x;

		if (threadId < vars.total_elementos_vetores_menores)
		{
			int chunk_atual = 0;

			int elemento_chave = vars.inicio_vetores_menores[threadId];

			int* ptr_inicio = vars.inicio_espaco_busca + vars.total_elementos_vetores_maiores;
			int* ptr_fim = ptr_inicio + vars.total_intersecoes;

			vars.resultado[threadId] = elemento_chave;

			while (chunk_atual < vars.total_intersecoes)
			{
				int* ptr_vetor_maior = vars.inicio_vetores_menores + ptr_inicio[chunk_atual];
				int tam_vetor_maior = ptr_fim[chunk_atual] - ptr_inicio[chunk_atual];
				int retorno;

				contains(ptr_vetor_maior, tam_vetor_maior, elemento_chave, &retorno);

				if (retorno == INVALID_NUMBER)
				{
					// skip no laço
					chunk_atual = vars.total_intersecoes;
					vars.resultado[threadId] = INVALID_NUMBER;
				}
				else
				{
					// ir fazer busca binaria no proximo chunk
					chunk_atual++;
				}

			}
		}

		__syncthreads( );
	}

	inline bool buffer_cheio(Buffer* buffer, int qtde_novos_elementos)
	{
		return !(buffer->pos_buffer + qtde_novos_elementos < buffer->max_tam_buffer);
	}

	inline bool init_buffer(Buffer* b, int tam)
	{
		b->pos_buffer = 0;
		b->vetor = (int*)malloc(sizeof(int) * tam);
		b->max_tam_buffer = tam;

		if (!b->vetor)
			return false;

		return true;
	}

	inline void terminar_buffer(Buffer* buffer)
	{
		free(buffer->vetor);
		buffer->vetor = NULL;
		buffer->pos_buffer = 0;
	}

	inline int add_buffer(GpuVars* gpu_vars, int* entrada, int qtde)
	{
		if (qtde + gpu_vars->buffer.pos_buffer > gpu_vars->buffer.max_tam_buffer)
		{
			if (gpu::gpu_cheia(gpu_vars, (gpu_vars->total_numeros_na_gpu + qtde * sizeof(int))))
				return ERRO_GPU;

			esvaziar_buffer(gpu_vars);

		}

		// copia host2host
		cudaMemcpy(&gpu_vars->buffer.vetor[gpu_vars->buffer.pos_buffer], entrada, qtde * sizeof(int), cudaMemcpyHostToHost);
		gpu_vars->buffer.pos_buffer += qtde;

		return SUCESSO_COPIA;
	}

	inline int espaco_vazio_buffer(Buffer* buffer)
	{
		return buffer->max_tam_buffer - buffer->pos_buffer;
	}

	inline int espaco_vazio_gpu_bytes(GpuVars* gpu_vars)
	{
		return gpu_vars->total_bytes_alocados - gpu_vars->total_bytes_ocupados;
	}

	inline void contabilizar_numeros_na_gpu(GpuVars* gpu_vars, int numeros)
	{
		gpu_vars->total_numeros_na_gpu += numeros;
		gpu_vars->total_bytes_ocupados += numeros * sizeof(int);
	}

	inline void set_valor_array(int* ptr, int tam, int valor)
	{
		thrust::device_ptr<int> device_ptr(ptr);
		thrust::fill(device_ptr, device_ptr + tam, valor);
	}

	inline int replicar_tids(GpuVars* gpu_vars, int qtde_replicas, int tam, int k)
	{
		if (qtde_replicas == 1)
			return 1;

		int* ptr_gpu_tids = gpu_vars->inicio_memoria;
		int* final_replica = ptr_gpu_tids + tam;
		register int r = 1;

		while (r < qtde_replicas)
		{
			/* cudaMemcpy(final_replica, ptr_gpu_tids, sizeof(int)* tam*fator, cudaMemcpyDeviceToDevice);
			final_replica = final_replica + tam;
			r += fator;

			if ((qtde_replicas - r) > fator)
			fator = r;
			else
			fator = qtde_replicas - r; */

			cudaMemcpy(final_replica, ptr_gpu_tids, sizeof(int)* tam, cudaMemcpyDeviceToDevice);
			final_replica = final_replica + tam;
			r++;

		}

		contabilizar_numeros_na_gpu(gpu_vars, (qtde_replicas - 1) * tam);
		return r;
	}

	inline int replicar_tids(GpuVars* gpu_vars, int qtde_replicas, int tam)
	{
		if (qtde_replicas == 1)
			return 1;

		int* ptr_gpu_tids = gpu_vars->inicio_memoria;
		int* final_replica = ptr_gpu_tids + tam;
		int total_copia = 1;

		register int r = 1;

		while (r < qtde_replicas)
		{
			if (r < qtde_replicas - r)
				total_copia = r;
			else
				total_copia = qtde_replicas - r;

			cudaMemcpy(final_replica, ptr_gpu_tids, sizeof(int) * tam * total_copia, cudaMemcpyDeviceToDevice);

			final_replica = final_replica + tam * total_copia;
			r += total_copia;
		}

		contabilizar_numeros_na_gpu(gpu_vars, (qtde_replicas - 1) * tam);
		return r;
	}

	void exibir_buffer(Buffer buffer)
	{
		printf("\nBUFFER: %d: ", buffer.pos_buffer);
		for (int i = 0; i < buffer.pos_buffer; i++)
			cout << (buffer.vetor)[i] << ", ";
	}

	// se a GPU encherá com os novos bytes a colocar
	inline bool gpu_cheia(GpuVars* gpu_vars, int novos_bytes)
	{
		return gpu_vars->total_bytes_alocados < (gpu_vars->total_bytes_ocupados + novos_bytes);
	}

	size_t gpu_bytes_livres(GpuVars* gpu_vars)
	{
		return gpu_vars->total_bytes_alocados - gpu_vars->total_bytes_ocupados;
	}

	void set_numeros_na_gpu(GpuVars* gpu_vars, int numeros)
	{
		gpu_vars->total_numeros_na_gpu = numeros;
		gpu_vars->total_bytes_ocupados = numeros * sizeof(int);
	}

	inline int* gpu_stream_compaction(int* gpu_esparse_array, int* gpu_compact_array, int tamanho, int* tamanho_compact_array)
	{
		thrust::device_ptr<int> esparse(gpu_esparse_array);
		thrust::device_ptr<int> compact(gpu_compact_array);

		thrust::device_ptr<int> retorno = thrust::copy_if(esparse, esparse + tamanho, compact, is_not_invalid());

		int* ptr_retorno = thrust::raw_pointer_cast(retorno);

		*tamanho_compact_array = (ptr_retorno - gpu_compact_array);

		return ptr_retorno;

	}

	__global__ void
		sequential_stream_compaction(int* gpu_esparse_array, int* gpu_compact_array, int tamanho, int* tamanho_novo)
	{
		int threadId = (blockIdx.x * blockDim.x) + threadIdx.x;
		int t = 0;

		if (!threadId)
			for (int i = 0; i < tamanho; i++)
				if (gpu_esparse_array[i] != INVALID_NUMBER)
				{
					gpu_compact_array[t] = gpu_esparse_array[i];
					t++;
				}

		*tamanho_novo = t;
	}

	// copia o buffer para a memoria da GPU de forma a concatenar os dados
	inline int esvaziar_buffer(GpuVars* gpu_vars)
	{
		Buffer* buffer = &(gpu_vars->buffer);

		if (buffer->pos_buffer == 0)
			return ERRO_BUFFER;

		int* ptr_gpu = gpu_vars->inicio_memoria + gpu_vars->total_numeros_na_gpu;

		if (copiar_cpu_para_gpu(ptr_gpu, buffer->vetor, buffer->pos_buffer * sizeof(int)))
		{
			contabilizar_numeros_na_gpu(gpu_vars, buffer->pos_buffer);
			buffer->pos_buffer = 0;

			return SUCESSO_COPIA;
		}

		return ERRO_BUFFER;

	}

	void resetar_buffer(Buffer* buffer)
	{
		buffer->pos_buffer = 0;
	}

	__global__ void kernel_exibir_memoria(int* ptr, int qtde)
	{
		printf("\n# Memoria GPU [%d]:  ", qtde);
		for (int e = 0; e < qtde; e++)
			printf("%d, ", ptr[e]);
	}

	void exibir_memoria(int* ptr_memoria_gpu, int qtde)
	{
		gpu::kernel_exibir_memoria << <1, 1 >> >(ptr_memoria_gpu, qtde);
		cudaDeviceSynchronize();
	}

	// calcula o maximo de replicacoes feitas para aumentar o paralelismo
	inline int calcular_replicacao(GpuVars* gpu_vars, Cuboid* cuboid, map<AttVal, vector<int>>::iterator* it_attval_ref, int freq_tids, int tipo, int qtde_celulas)
	{
		int limitante_threads = gpu_vars->max_blocks_por_multiprocessador * gpu_vars->prop.multiProcessorCount * gpu_vars->block_size;
		int total_valores_atributos = 0;
		int bytes_a_ocupar = 0;

		AttVal att_val_all = create_attribute_value(ALL);

		while (*it_attval_ref != cuboid->values.end())
		{
			if (attribute_values_equals((*it_attval_ref)->first, att_val_all) == false)
			{
				vector<int> tids_attval_dim_inq = get_tids(cuboid, (*it_attval_ref)->first);
				int novos_bytes;

				if (tipo == TAMANHO_FIXO)
				{
					novos_bytes = ((tids_attval_dim_inq.size() + 2) + (freq_tids * 2)) * sizeof(int);
				}
				else if (tipo == TAMANHO_VARIAVEL)
					novos_bytes = ((tids_attval_dim_inq.size( ) + 2) + (freq_tids * 2) + (qtde_celulas * 2)) * sizeof(int);

				if (gpu::gpu_cheia(gpu_vars, bytes_a_ocupar + novos_bytes) == false)
				{
					bytes_a_ocupar += novos_bytes;
					total_valores_atributos++;
				}
				else break;

			}

			(*it_attval_ref)++;

		}

		// threshold de threads
		int max_replicas_tmp;

		max_replicas_tmp = (gpu::espaco_vazio_gpu_bytes(gpu_vars) - bytes_a_ocupar) / (freq_tids * sizeof(int));
		max_replicas_tmp = (limitante_threads / freq_tids) < max_replicas_tmp ? (limitante_threads / freq_tids) : max_replicas_tmp;

		return total_valores_atributos < max_replicas_tmp ? total_valores_atributos : max_replicas_tmp;

	}

	__global__ void remover_esparsidade(int*resposta, int total_celulas_novas, int total_celulas_derivadas, int* tam_celulas_gpu, int* inicio_derivada, int* fim_derivada, int* saida)
	{
		int threadId = (blockDim.x * blockIdx.x) + threadIdx.x;

		if (threadId < total_celulas_novas)
		{
			int indice_inicial = ((threadId / total_celulas_derivadas) * fim_derivada[total_celulas_derivadas - 1]) + inicio_derivada[threadId % total_celulas_derivadas];
			int indice_final = indice_inicial + fim_derivada[threadId % total_celulas_derivadas] - inicio_derivada[threadId % total_celulas_derivadas];
			int total_tids = 0;

			while (indice_inicial < indice_final)
			{
				if (resposta[indice_inicial] == INVALID_NUMBER)
				{
					int percorre = indice_inicial + 1;

					while (percorre < indice_final)
					{
						if (resposta[percorre] != INVALID_NUMBER)
						{
							resposta[indice_inicial] = resposta[percorre];
							indice_inicial++;
							total_tids++;
						}

						percorre++;

					}

					if (percorre == indice_final)
						indice_inicial = indice_final;

				}
				else
				{
					indice_inicial++;
					total_tids++;
				}
			}

			tam_celulas_gpu[threadId] = total_tids;

		}

		__syncthreads( );

		// memcpy sequenciais
		if (threadId == 0)
		{
			int deslocamento = 0;

			for (int celula = 0; celula < total_celulas_novas; celula++)
			{
				if (tam_celulas_gpu[celula] > 0)
				{
					int indice_inicial = ((celula / total_celulas_derivadas) * fim_derivada[total_celulas_derivadas - 1]) + inicio_derivada[threadId % total_celulas_derivadas];
					memcpy(saida + deslocamento, resposta + indice_inicial, tam_celulas_gpu[celula] * sizeof(int));
					deslocamento += tam_celulas_gpu[celula];
				}
			}
		}

	}

	__global__ void remover_esparsidade(int*resposta, int total_celulas_novas, int total_celulas_derivadas, int* tam_celulas_gpu, int total_tids_instanciadas, int* saida)
	{
		int threadId = (blockDim.x * blockIdx.x) + threadIdx.x;

		if (threadId < total_celulas_novas)
		{
			int indice_inicial = threadId * total_tids_instanciadas;
			int indice_final = indice_inicial + total_tids_instanciadas;
			int total_tids = 0;

			while (indice_inicial < indice_final)
			{
				if (resposta[indice_inicial] == INVALID_NUMBER)
				{
					int percorre = indice_inicial + 1;

					while (percorre < indice_final)
					{
						if (resposta[percorre] != INVALID_NUMBER)
						{
							resposta[indice_inicial] = resposta[percorre];
							indice_inicial++;
							total_tids++;
						}

						percorre++;

					}

					if (percorre == indice_final)
						indice_inicial = indice_final;

				}
				else
				{
					indice_inicial++;
					total_tids++;
				}
			}

			tam_celulas_gpu[threadId] = total_tids;

		}

		__syncthreads( );

		// memcpy sequenciais
		if (threadId == 0)
		{
			int deslocamento = 0;

			for (int celula = 0; celula < total_celulas_novas; celula++)
			{
				if (tam_celulas_gpu[celula] > 0)
				{
					int indice_inicial = celula * total_tids_instanciadas;
					memcpy(saida + deslocamento, resposta + indice_inicial, tam_celulas_gpu[celula] * sizeof(int));
					deslocamento += tam_celulas_gpu[celula];
				}
			}
		}
	}

	inline int* remover_esparsidade(GpuVars* gpu_vars, int*resposta, int total_celulas_novas, int total_celulas_derivadas, int* tam_celulas_gpu, int* inicio_derivada, int* fim_derivada, int* saida)
	{
		int tamanho_grid, tamanho_bloco;

		gpu::balancear_carga(gpu_vars, total_celulas_novas, &tamanho_bloco, &tamanho_grid);

		int* tamanho_celulas_resultantes = (int*)malloc(sizeof(int) * total_celulas_novas);

		remover_esparsidade << <tamanho_grid, tamanho_bloco >> >
			(resposta, total_celulas_novas, total_celulas_derivadas, tam_celulas_gpu, inicio_derivada, fim_derivada, saida);

		cudaMemcpy(tamanho_celulas_resultantes, tam_celulas_gpu, sizeof(int) * total_celulas_novas, cudaMemcpyDeviceToHost);

		return tamanho_celulas_resultantes;
	}

	inline int* remover_esparsidade(GpuVars* gpu_vars, int*resposta, int total_celulas_novas, int total_celulas_derivadas, int* tam_celulas_gpu, int total_tids_instanciadas, int* saida)
	{
		int tamanho_grid, tamanho_bloco;

		gpu::balancear_carga(gpu_vars, total_celulas_novas, &tamanho_bloco, &tamanho_grid);

		int* tamanho_celulas_resultantes = (int*) malloc(sizeof(int) * total_celulas_novas);

		remover_esparsidade << <tamanho_grid, tamanho_bloco >> >
			(resposta, total_celulas_novas, total_celulas_derivadas, tam_celulas_gpu, total_tids_instanciadas, saida);

		cudaDeviceSynchronize( );
		cudaMemcpy(tamanho_celulas_resultantes, tam_celulas_gpu, sizeof(int) * total_celulas_novas, cudaMemcpyDeviceToHost);

		return tamanho_celulas_resultantes;
	}
}

#endif