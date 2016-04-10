#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gpucubing.h"

const char* path_log_amostras = "ultima_saida.txt";

void exibir_parametrizacao_correta()
{
	cout << "\n\nParametrizacao errada! As corretas sao:";
	cout << "\n\n1. Processanto lote:";
	cout << "\n<base> <numero-colunas> <path-queries>\n";
	cout << "\n\n2. Shell iterativo:";
	cout << "\n<base> <numero-colunas> -s\n";
	cout << "\n\n3. Repetir consulta:";
	cout << "\n<base> <numero-colunas> -b <comando> <repeticoes>\n";
}

void salvar_amostras(char* base, char* comando, vector<long> amostras)
{
	FILE* arquivo = fopen(path_log_amostras, "w");

	cout << "\n\n\nGerando output...";

	if (arquivo)
	{
		cout << "\nSalvando saida...";
		fprintf(arquivo, "Base: %s\n\nQuery:\n\n%s\n\nTotal amostras: %d", base, comando, amostras.size( ));
		fprintf(arquivo, "\n\n");

		for (int i = 0; i < amostras.size( ); i++)
			fprintf(arquivo, "%d, ", amostras[i]);

		fclose(arquivo);

		cout << "\nSaida salva!";
	}
	else
	{
		cout << "\nErro! Nao foi possivel criar o arquivo para salvar as amostras...";
	}
}

bool comparar_resultados(map<string, long>* mapa_isaias, map<string, long>* mapa_han)
{
	map<string, long>::iterator iterator = mapa_isaias->begin( );

	cout << "\n\n\n\n";

	while (iterator != mapa_isaias->end())
	{
		cout << "\n    " << iterator->first;

		if (mapa_han->find(iterator->first) != mapa_han->end())
		{
			if ((*mapa_han)[iterator->first] != (*mapa_isaias)[iterator->first])
			{
				cout << "\n    ERRADO!   Isaias - Measure: " << (*mapa_isaias)[iterator->first]
					<< " | Han - Measure: " << (*mapa_han)[iterator->first];
				cout << "\n\n_______________________________________________________________________\n\n";

				return false;
			}
			else
			{
				//	cout << "\n    Isaias: " << (*mapa_isaias)[iterator->first]
				//	<< " | Han: " << (*mapa_han)[iterator->first];
			}
		}
		else
		{
			cout << "\n    ERRADO! Value não encontrado em ambos hashs!";
			cout << "\n\n_______________________________________________________________________\n\n";

			return false;
		}

		iterator++;

	}

	return true;
	
}

int _main_(int argc, char** argv)
{
	int* ptr_destino;
	int* ptr_source;
	int* inicio;
	int* fim;
	int total_celulas = 8;
	int* tamanho_celulas;

	cudaMalloc(&ptr_destino, 40 * sizeof(int));
	cudaMalloc(&ptr_source, 40 * sizeof(int));
	cudaMalloc(&tamanho_celulas, total_celulas * sizeof(int));

	gpu::exibir_memoria(ptr_destino, 40);

	vector<int> dados = { 0, 1, -1, 2, 3, -1, -1, 4, 5, -1, -1, 6, 7, -1, -1, 8, 9, -1, -1, 10, 11, -1, -1, 12, -1, -1, 1, -1, -1, -1, -1, -1 };
	vector<int> inicio_cpu = { 0, 4, 8, 12, 16, 20, 24, 28 };
	vector<int> fim_cpu = { 4, 8, 12, 16, 20, 24, 28, 32 };

	cudaMalloc(&inicio, inicio_cpu.size() * sizeof(int));
	cudaMalloc(&fim, fim_cpu.size() * sizeof(int));

	cudaMemcpy(ptr_source, &dados[0], dados.size() * 4, cudaMemcpyHostToDevice);

	cudaMemcpy(inicio, &inicio_cpu[0], inicio_cpu.size() * 4, cudaMemcpyHostToDevice);
	cudaMemcpy(fim, &fim_cpu[0], fim_cpu.size() * 4, cudaMemcpyHostToDevice);

	printf("\nTotal copias: %d", dados.size() / 4);

	/*int* d = gpu::remover_esparsidade(1, 32, ptr_source, tamanho_celulas, inicio, fim, total_celulas, total_celulas, ptr_destino);

	int sum = 0;
	printf("\n\n");
	for (int i = 0; i < total_celulas; i++)
	{
	sum += d[i];
	printf("%d, ", d[i]);
	}

	printf("\n\nMemoria GPU: ");
	gpu::exibir_memoria(ptr_destino, sum);
	getchar( );



	free(d);
	d = NULL;
	*/
	return 0;
}

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		exibir_parametrizacao_correta();
		return 0;
	}

	vector<string> array_queries;
	int num_cpu_cores;
	string path_file_queries;
	ifstream myfile;
	string line;

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	num_cpu_cores = sysinfo.dwNumberOfProcessors;

	if (argc == 5)
	{
		if (!strcmp(argv[3], "-c"))
		{
			string comando = string(argv[4]);
			array_queries.push_back(comando);
		}
	}
	else if (strcmp(argv[3], "-s"))
	{
		path_file_queries = string(argv[4]);
		myfile = ifstream(path_file_queries);

		while (getline(myfile, line))
			array_queries.push_back(line);

		myfile.close( );
	}

	Cubo* cubo;
	Fragmentos* fragmentos;
	int numero_colunas = atoi(argv[2]);

	cubo = carregar_cubo_base(argv[1], numero_colunas);

	cout << "\nFragmentando cubo...";
	fragmentos = fragmentar_cubo(cubo);
	cout << "\nCubo fragmentado.";

	deletar_cubo_base(&cubo);

	string path_file_isaias = string(argv[1]);
	map<string, long> map_isaias;

	#if EXECUTAR_HAN
	map<string, long> map_han;
	char* path_base_han = argv[1];

	// iniciando variaveis do frag-cube
	frag_init_vars(path_base_han, 1, 1, numero_colunas);

	long frag_cubing_time = 0;

	#endif

	long gpu_query_cubing_time = 0;

	// Obtendo amostras
	if (argc == 6 && !strcmp(argv[3], "-b"))
	{
		string comando = string(argv[4]);
		int repeticoes = atoi(argv[5]);
		vector<long> amostras;

		gpu::GpuVars gpu_vars;

		gpu_vars.max_kernels_simultaneos = MAX_KERNELS_SIMULTANEOS;
		gpu_vars.num_cpu_cores = num_cpu_cores;

		gpu::alocar_gpu(&gpu_vars);
		gpu::init_buffer(&gpu_vars.buffer, MAX_TAM_BUFFER);

		for (int r = 0; r < repeticoes; r++)
		{
			cout << "\nPression <enter> para realizar a execucao da outra query...";
			fflush(stdin);
			getchar( );

			long gpu_query_cubing_time_tmp = 0;
			cout << "\nExecutando algoritmo GPU Cubing:" << endl;
			map_isaias = gpucubing::execute_query(&gpu_vars, comando, &gpu_query_cubing_time_tmp, fragmentos);
			gpu_query_cubing_time += gpu_query_cubing_time_tmp;
			amostras.push_back(gpu_query_cubing_time_tmp);
		}

		gpu::terminar_buffer(&gpu_vars.buffer);
		gpu::desalocar_gpu(&gpu_vars);

		cout << "\n\nQuery: " << comando;
		cout << "\n\nBase: " << string(argv[1]);
		cout << "\n\nAmostras:\n";

		for (int a = 0; a < amostras.size(); a++)
			cout << amostras[a] << ", ";

		salvar_amostras(argv[1], argv[4], amostras);

		cout << "\n\n";

		cudaDeviceReset( );

		return 0;

	}
	else if (argc == 4 && !strcmp(argv[3], "-s"))		// shell interativo
	{
		string comando = "";
		gpu::GpuVars gpu_vars;

		gpu_vars.max_kernels_simultaneos = MAX_KERNELS_SIMULTANEOS;
		gpu_vars.num_cpu_cores = num_cpu_cores;

		gpu::alocar_gpu(&gpu_vars);
		gpu::init_buffer(&gpu_vars.buffer, MAX_TAM_BUFFER);

		while (comando.compare("exit") != 0)
		{
			cout << "\ngpucubing >> ";
			getline(cin, comando);

			comando = replace_char(comando, ' ', ':');
			
			if (comando.compare("exit") && comando.compare(""))
			{
				long gpu_query_cubing_time_tmp = 0;

				cout << "\nExecutando algoritmo GPU Cubing:" << endl;
				map_isaias = gpucubing::execute_query(&gpu_vars, comando, &gpu_query_cubing_time_tmp, fragmentos);
				gpu_query_cubing_time += gpu_query_cubing_time_tmp;
			}

		}

		gpu::terminar_buffer(&gpu_vars.buffer);
		gpu::desalocar_gpu(&gpu_vars);

	}
	else
	{

		for (int q = 0; q < array_queries.size( ); q++)
		{
			long gpu_query_cubing_time_tmp = 0;
			cout << "\nExecutando algoritmo Isaias:" << endl;
			map_isaias = gpucubing::execute_query(array_queries[q], &gpu_query_cubing_time_tmp, fragmentos);
			gpu_query_cubing_time += gpu_query_cubing_time_tmp;

			#if EXECUTAR_HAN
			long frag_cubing_time_tmp = 0;
			cout << "\nExecutando algoritmo Frag-Cubing:" << endl;
			map_han = execute_query(array_queries[q], &frag_cubing_time_tmp);
			frag_cubing_time += frag_cubing_time_tmp;

			cout << "    Tempo total Isaias: " << gpu_query_cubing_time << " ms." << endl;
			getchar();

			if (map_han.size() == map_isaias.size())
				if (comparar_resultados(&map_isaias, &map_han))
				{
					cout << "\n\n    A consulta\n\n" << array_queries[q] << "\n\n    foi executada corretamente!\n\n" << endl;
					cout << "    Tempo total Isaias: " << gpu_query_cubing_time << " ms." << endl;
					cout << "    Tempo total Han: " << frag_cubing_time << " ms." << endl;
					cout << "\n\n";
				}
				else
					cout << "\n\n    A consulta " << array_queries[q] << " possui resultados diferentes!" << endl;
			else
				cout << "\n\nERRO! O tamanho do retorno de ambos algoritmos retornam valores distintos!";
			#endif

		} // for
	}

	#if EXECUTAR_HAN
	frag_delete_vars( );
	#else
	cout << "\n    Tempo total Isaias: " << gpu_query_cubing_time << " ms." << endl;
	#endif

	cudaDeviceReset( );

	return 0;

}
