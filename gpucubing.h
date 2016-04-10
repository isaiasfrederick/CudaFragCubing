#ifndef GPU_QUERY_CUBING_H
#define GPU_QUERY_CUBING_H

#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "query.h"
#include "cubo.h"
#include "query.h"
#include "frag.h"
#include "gpu.h"
#include "constantes.h"

using namespace std;

namespace gpucubing
{
	inline map<string, long> execute_query(string query_str, long* tempo_gasto, Fragmentos* fragmentos)
	{
		gpu::GpuVars gpu_vars;
		Query query;
		ContadoresTempo contadores;
		map<string, long> map_resultado;

		gpu_vars.max_kernels_simultaneos = MAX_KERNELS_SIMULTANEOS;

		gpu::alocar_gpu(&gpu_vars);
		gpu::init_buffer(&gpu_vars.buffer, MAX_TAM_BUFFER);

		create_query(fragmentos, query_str, &query);

		contadores.elapsed_time = 0;
		marcar_tempo_inicial(&contadores);

		processar_consulta(&gpu_vars, fragmentos, query, &map_resultado);

		marcar_tempo_final(&contadores);
		contabilizar_tempo_gasto(&contadores);

		limpar_query(&query);

		*tempo_gasto = contadores.elapsed_time;

		gpu::terminar_buffer(&gpu_vars.buffer);
		gpu::desalocar_gpu(&gpu_vars);

		return map_resultado;
	}

	inline map<string, long> execute_query(gpu::GpuVars* gpu_vars, string query_str, long* tempo_gasto, Fragmentos* fragmentos)
	{
		Query query;
		ContadoresTempo contadores;
		map<string, long> map_resultado;

		create_query(fragmentos, query_str, &query);

		contadores.elapsed_time = 0;
		marcar_tempo_inicial(&contadores);

		processar_consulta(gpu_vars, fragmentos, query, &map_resultado);

		marcar_tempo_final(&contadores);
		contabilizar_tempo_gasto(&contadores);

		limpar_query(&query);

		gpu_vars->buffer.pos_buffer = 0;
		gpu_vars->total_numeros_na_gpu = 0;

		*tempo_gasto = contadores.elapsed_time;

		return map_resultado;
	}
}


#endif