#ifndef FRAG_H
#define FRAG_H
#include <iostream> 
#include <sstream> 
#include <fstream> 
#include <string>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <windows.h>

#define DEBUG false
#define ONLINE_DEBUG false
#define SHRINK 25
#define MAX_LINE_LENGTH 4096
#define LOOK 10

#define NUM_QUERIES 1000
#define NUM_INQ 4
#define NUM_INS 3
#define NUM_DIMS 20
#define CARD 10

#define tipoint64 long long

using namespace std;

struct cell {
	tipoint64 num_trans;
	tipoint64 capacity;
	tipoint64 *trans;
};

struct node {
	tipoint64 name;
	node **child;
	tipoint64 num_cells;
	cell **cells;
};

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone2
{
	tipoint64  tz_minuteswest; /* minutes W of Greenwich */
	tipoint64  tz_dsttime;     /* type of dst correction */
};

	// prototypes
	void read_datafile(char*);
	void fragment(node*, cell*, vector< tipoint64 >, tipoint64*, tipoint64);
	void online_fragment(node*, node*, cell*, vector< tipoint64 >, tipoint64*, tipoint64, tipoint64*, tipoint64*, bool, map<string,long>*);
	void intersect(cell*, cell*, cell*);
	void io_intersect(cell*, cell*, cell*);
	void elapsed_time(timeval*);
	void mark_time(timeval*);
	void interfaceM(string, map<string, long>*, long*);
	void execute_query(string, bool, bool, map<string, long>*);
	tipoint64 stoiM(const string &s);
	void multi_intersect(cell*, cell **, tipoint64);
	cell* cell_of(int*, tipoint64);
	cell* online_cell_of(node*, tipoint64*, tipoint64, tipoint64*);
	void free_tree(node*, tipoint64*, tipoint64);
	string itos(tipoint64);
	void exp_time_taken(timeval*);
	tipoint64 gettimeofdayM(struct timeval *tv, struct timezone2 *tz);

	// funcoes inseridas pelo Isaias
	void frag_delete_vars();
	long final_time_taken(timeval *st);
	string replace_char(string, char, char);
	int frag_init_vars(char* path_file, int min_sup_arg, int frag_size_arg, int numero_coluna_arg);

	// global vars
	tipoint64 minsup;
	tipoint64 n_dimensions;
	tipoint64 n_tuples;
	tipoint64 *cardinality;
	node *root;
	timeval *prev_time;
	tipoint64 n_frags;
	tipoint64 frag_size;
	tipoint64 *frag;
	tipoint64 *frag_index;
	tipoint64 num_rows;
	tipoint64 _IOCOUNT;
	tipoint64 _IOMAX;
	tipoint64 _CUBEIO;

	tipoint64 numero_coluna;

	string replace_char(string str, char ch1, char ch2) {
		for (int i = 0; i < str.length(); ++i) {
			if (str[i] == ch1)
				str[i] = ch2;
		}

		return str;
	}


	// output file
	ofstream outfile("output.full");

	void frag_delete_vars( )
	{
		if (prev_time != NULL)
			delete prev_time;
		if (cardinality != NULL)
			delete cardinality;

	}

	int frag_init_vars(char* path_file, int min_sup_arg, int frag_size_arg, int numero_coluna_arg)
	{
		tipoint64 *cell_name;
		vector<tipoint64 > cuboid_name;

		// ISAIAS if (argc != 7) {
		// ISAIAS 	cout << argv[0] << " <datafile> <minsup> <fragment size> <numero_coluna> <query>" << endl;
		// ISAIAS 	exit(1);
		// ISAIAS }

		// ISAIAS cout << "Initializing";

		// mark start time
		prev_time = new timeval();
		gettimeofdayM(prev_time, NULL);

		// init vars
		n_dimensions = 0;
		minsup = min_sup_arg;
		frag_size = frag_size_arg;
		numero_coluna = numero_coluna_arg;

		// read data file
		read_datafile(path_file);

		cell_name = new tipoint64[n_dimensions];
		memset(cell_name, -1, sizeof(tipoint64) * n_dimensions);

		// ISAIAS cout << "Computing shell fragments";

		// for each fragment
		for (tipoint64 i = 0; i < root->num_cells; i++) {

			// for each node in the fragment
			for (tipoint64 j = 0; j < root->child[i]->num_cells; j++) {
				cuboid_name.push_back(root->child[i]->child[j]->name);
				fragment(root->child[i]->child[j], NULL, cuboid_name,
					cell_name, root->child[i]->num_cells);
				cuboid_name.pop_back();
			}
		}

		// tempo de fragmentação
		elapsed_time(prev_time);

		// ISAIAS cout << endl << n_dimensions << " dimensional data loaded." << endl;
		// ISAIAS cout << n_tuples << " tuples read." << endl;
		// ISAIAS cout << n_frags << " shell fragments of size " << frag_size << 
		// ISAIAS	" constructed." << endl << endl;

		// clean up a bit
		delete[] cell_name;

		/*long final_time = 0;

		for (int q = 0; q < array_queries.size(); q++)
		{
			// online component
			string comando_novo = "q " + replace_char(array_queries[q], ':', ' ');
			interfaceM(comando_novo, mapa, &final_time);
		}

		cout << "\n    Tempo total necessario: " << final_time << " ms.";*/

		// ISAIAS han_delete_vars();

		return 0;
	}


	// funcao criada por isaias
	map<string, long> execute_query(string query_str, long* final_time)
	{

		map<string, long> mapa_retorno;
		string comando_novo = "q " + replace_char(query_str, ':', ' ');

		interfaceM(comando_novo, &mapa_retorno, final_time);

		return mapa_retorno;

	}



	void read_datafile(char *filename)
	{
		FILE *f;
		tipoint64 temp, index, size;
		// ISAIAS tipoint64 word;
		//char line[MAX_LINE_LENGTH], *word;
		cell *c;

		// open data file
		if ((f = fopen(filename, "rb")) == NULL) {
			printf("Error: cannot open file %s.\n", filename);
			exit(-1);
		}

		// read first line of data file 
		//fgets(line, MAX_LINE_LENGTH, f);

		// read in the number of tuples
		/*word = strtok(line, " ");
		n_tuples = atoi(word);

		n_dimensions = 0;
		word = strtok(NULL, " ");*/
		fread(&n_tuples, sizeof(int), 1, f);

		// count the number of dimensions (numero de colunas)
		/*while (word != NULL && (0 != strcmp(word, "\n"))) {
		n_dimensions++;
		word = strtok(NULL, " ");
		}*/
		n_dimensions = numero_coluna;

		// rewind to beginning of line and skip to next word
		/*rewind(f);
		fgets(line, MAX_LINE_LENGTH, f);
		word = strtok(line, " ");*/

		// store the cardinality and offset of each dimension
		cardinality = new tipoint64[n_dimensions];

		// read cardinalities
		int temp2;
		for (tipoint64 i = 0; i < n_dimensions; i++) {
			//word = strtok(NULL, " ");
			//cardinality[i] = atoi(word) + 1;
			fread(&temp2, sizeof(int), 1, f);
			cardinality[i] = temp2 + 1;
		}

		if (DEBUG) {
			// ISAIAS cout << endl << "cardinality: ";
			for (tipoint64 i = 0; i < n_dimensions; i++)
				cout << cardinality[i] << " ";
			cout << endl;
		}

		// establish the shell fragments
		n_frags = (tipoint64)ceil((double)n_dimensions / frag_size);

		// create the root of the entire forest
		root = new node;
		root->name = -1000;
		root->child = new node*[n_frags];
		root->num_cells = n_frags;
		root->cells = NULL;

		index = 0;

		// allocate each cuboid group tree
		for (tipoint64 i = 0; i < n_frags; i++) {

			size = frag_size;

			// last shell fragment.  this fragment needs special care
			// because its size can be smaller than frag_size
			if (i == n_frags - 1 && n_dimensions % frag_size != 0) {
				size = n_dimensions % frag_size;
			}

			// create the root of this particular fragment
			root->child[i] = new node;
			root->child[i]->name = -100;
			root->child[i]->num_cells = size;
			root->child[i]->cells = NULL;
			root->child[i]->child = new node*[size];

			// now create all the nodes in this fragment
			for (tipoint64 j = 0; j < size; j++) {
				root->child[i]->child[j] = new node;
				root->child[i]->child[j]->name = index;
				root->child[i]->child[j]->num_cells = cardinality[index];
				root->child[i]->child[j]->cells = new cell*[cardinality[index]];
				root->child[i]->child[j]->cells[0] = NULL;
				root->child[i]->child[j]->child = NULL;

				// note that this starts at 1 because the generator program
				// starts at 1.
				for (tipoint64 k = 1; k < cardinality[index]; k++) {
					root->child[i]->child[j]->cells[k] = new cell;
					root->child[i]->child[j]->cells[k]->num_trans = 0;

					// allocate space assuming that the data is uniform
					// distribution
					root->child[i]->child[j]->cells[k]->capacity = n_tuples /
						(cardinality[index] - 1);
					root->child[i]->child[j]->cells[k]->trans = new
						tipoint64[root->child[i]->child[j]->cells[k]->capacity];
				}

				index++;
			}
		}

		index = -1;

		// keeps track of which fragment each dimension belongs to
		frag = new tipoint64[n_dimensions];

		// keeps track of where in each fragment the dimension is
		frag_index = new tipoint64[n_dimensions];

		// create the fragment indices of each dimension
		for (tipoint64 i = 0; i < n_dimensions; i++) {
			if (i % frag_size == 0)
				index++;

			frag[i] = index;
			frag_index[i] = i % frag_size;
		}

		elapsed_time(prev_time);
		// ISAIAS cout << "Reading input file";

		// read in all data into the one-dimensional cells
		for (tipoint64 i = 0; i < n_tuples; i++) {

			//word = strtok(line, " ");

			for (tipoint64 j = 0; j < n_dimensions; j++) {
				//temp = atoi(word);
				fread(&temp2, sizeof(int), 1, f);
				temp = temp2;

				// find the cell for this entry
				c = root->child[frag[j]]->child[frag_index[j]]->cells[temp];

				// no transaction list yet, make one
				if (c->capacity == 0) {
					c->capacity = 4;
					c->trans = new tipoint64[4];
					c->num_trans = 1;
					c->trans[0] = i;
				}
				else {

					// need to expand the transaction list
					if (c->num_trans >= c->capacity) {

						tipoint64 *old = c->trans;
						tipoint64 old_size = c->capacity;

						// double capacity unless it's too big
						c->capacity *= 2;
						if (c->capacity > n_tuples) c->capacity = n_tuples;

						// copy over the older stuff
						c->trans = new tipoint64[c->capacity];
						memcpy(c->trans, old, sizeof(tipoint64) * old_size);

						// delete the older stuff
						delete[] old;
					}

					c->trans[c->num_trans] = i;
					c->num_trans++;
				}

				//word = strtok(NULL, " ");
				//fread(&temp2, sizeof(int), 1, f);
			}
		}

		// prune transaction lists
		for (tipoint64 i = 0; i < root->num_cells; i++) {

			for (tipoint64 j = 0; j < root->child[i]->num_cells; j++) {

				for (tipoint64 k = 0; k < root->child[i]->child[j]->num_cells; k++) {

					if (root->child[i]->child[j]->cells[k] != NULL) {

						c = root->child[i]->child[j]->cells[k];

						// iceberg pruning
						if (c->num_trans > 0 && c->num_trans < minsup) {
							delete[] c->trans;
							c->capacity = 0;
							c->num_trans = 0;
						}

						// in case the transaction list capacity is too big,
						// shrink it down
						if (c->capacity - c->num_trans > SHRINK) {

							tipoint64 *old = c->trans;

							// resize and copy over the older stuff
							c->capacity = c->num_trans;
							c->trans = new tipoint64[c->num_trans];
							memcpy(c->trans, old, sizeof(tipoint64) * c->num_trans);

							// delete the older stuff
							delete[] old;
						}
					}
				}
			}
		}

		elapsed_time(prev_time);

		// done with all file reading.  close data file
		fclose(f);

		// output base cuboid
		/*
		for (int i = 0; i < n_dimensions; i++) {
		outfile << "* ";
		}
		outfile << ": " << n_tuples << endl;
		*/

		// print out the one-dimensional cells
		if (DEBUG) {

			for (tipoint64 i = 0; i < root->num_cells; i++) {

				for (tipoint64 j = 0; j < root->child[i]->num_cells; j++) {

					for (tipoint64 k = 0; k < root->child[i]->child[j]->num_cells; k++) {

						if (root->child[i]->child[j]->cells[k] == NULL)
							cout << endl;
						else {
							cout << "(" << i << ", " << j << ", " << k << ") -> ";
							for (tipoint64 u = 0; u < root->child[i]->child[j]->cells[k]->num_trans; u++)
								cout << root->child[i]->child[j]->cells[k]->trans[u] << " ";
							cout << endl;
						}
					}
				}
			}

			// ISAIAS cout << "----------------------------" << endl;
		}

		return;
	}

	void fragment(node *my_root, cell *parent_cell, vector< tipoint64 > cuboid_name,
		tipoint64 *cell_name, tipoint64 n_dim)
	{
		tipoint64 num_children, cell_offset, temp, last_name;
		cell *c;

		if (DEBUG) {
			cout << endl << "cut starting at" << endl;
			cout << "\tcuboid: ";
			for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++)
				cout << cuboid_name[i] << " ";
			cout << endl;
			cout << "\tcell  : ";
			for (tipoint64 i = 0; i < n_dimensions; i++) {
				if (cell_name[i] == -1)
					cout << "* ";
				else
					cout << cell_name[i] << " ";
			}
			cout << endl;
		}

		cell_offset = 0;
		last_name = cuboid_name[cuboid_name.size() - 1];

		if (parent_cell != NULL) {

			// allocate my cells if it hasn't been done yet.
			if (my_root->num_cells <= 0) {
				my_root->num_cells = 1;

				for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++)
					my_root->num_cells *= cardinality[cuboid_name[i]];

				my_root->cells = new cell*[my_root->num_cells];

				if (DEBUG)
					cout << "created " << my_root->num_cells << " cells." << endl;
			}

			// compute cell offset
			for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++) {

				if (cell_name[cuboid_name[i]] != -1) {

					temp = 1;
					for (unsigned tipoint64 j = i + 1; j < cuboid_name.size(); j++)
						temp *= cardinality[cuboid_name[j]];

					cell_offset += temp * cell_name[cuboid_name[i]];
				}
			}

			if (DEBUG)
				cout << "cell_offset = " << cell_offset << endl;

			// aggregate all my cells with my parent cell
			for (tipoint64 i = cell_offset + 1; i < cardinality[last_name] +
				cell_offset; i++) {

				my_root->cells[i] = new cell;

				c = root->child[frag[last_name]]->child[frag_index[last_name]]->cells[i - cell_offset];

				// intersect 
				intersect(my_root->cells[i], parent_cell, c);

				if (DEBUG) {
					cout << "intersecting:" << endl;
					cout << "[" << parent_cell->num_trans << "]: ";
					for (tipoint64 j = 0; j < parent_cell->num_trans; j++)
						cout << parent_cell->trans[j] << " ";
					cout << endl;
					cout << "[" << c->num_trans << "]: ";
					for (tipoint64 j = 0; j < c->num_trans; j++)
						cout << c->trans[j] << " ";
					cout << endl;
					cout << "result: ";
					cout << "[" << my_root->cells[i]->num_trans << "]: ";
					for (tipoint64 j = 0; j < my_root->cells[i]->num_trans; j++)
						cout << my_root->cells[i]->trans[j] << " ";
					cout << endl;
				}
			}

			// reset name
			cell_name[last_name] = -1;
		}

		num_children = 0;

		// not at a leaf yet, need to create my children first
		if ((last_name % frag_size) < root->child[frag[my_root->name]]->num_cells - 1) {

			num_children = root->child[frag[my_root->name]]->num_cells - 1 -
				my_root->name % frag_size;

			// create my children nodes (if necessary)
			if (my_root->child == NULL) {

				my_root->child = new node*[num_children];

				for (tipoint64 i = 0; i < num_children; i++) {
					my_root->child[i] = new node;
					my_root->child[i]->name = my_root->name + i + 1;
					my_root->child[i]->num_cells = 0;
					my_root->child[i]->child = NULL;
				}

				if (DEBUG)
					cout << "created " << num_children << " children." << endl;
			}
		}

		// for all my cells
		for (tipoint64 i = cell_offset + 1; i < cardinality[last_name] +
			cell_offset; i++) {

			c = my_root->cells[i];

			// check iceberg condition
			if (c != NULL && c->num_trans >= minsup) {

				cell_name[last_name] = i - cell_offset;

				// output it
				/*
				for (int u = 0; u < n_dimensions; u++) {
				if (cell_name[u] == -1)
				outfile << "* ";
				else
				outfile << cell_name[u] << " ";
				}
				outfile << ": " << c->num_trans << endl;
				*/

				if (DEBUG) {
					cout << "cell ";
					for (tipoint64 u = 0; u < n_dimensions; u++) {
						if (cell_name[u] == -1)
							cout << "* ";
						else
							cout << cell_name[u] << " ";
					}
					cout << endl;
				}

				// recursive call for all my children (if any)
				for (tipoint64 j = 0; j < num_children; j++) {
					cuboid_name.push_back(my_root->child[j]->name);
					fragment(my_root->child[j], c, cuboid_name, cell_name, n_dim);
					cuboid_name.pop_back();
				}
			}

			cell_name[last_name] = -1;
		}
	}

	void intersect(cell *target, cell *a, cell *b)
	{
		cell *smallM, *big;
		tipoint64 size;

		if (a->num_trans < b->num_trans) {
			size = a->num_trans;
			smallM = a;
			big = b;
		}
		else {
			size = b->num_trans;
			smallM = b;
			big = a;
		}

		// not enough for minsup
		if (size < minsup) {
			target->trans = NULL;
			target->num_trans = 0;
			target->capacity = 0;
			return;
		}

		target->trans = new tipoint64[size];
		target->num_trans = 0;
		target->capacity = size;

		tipoint64 j = 0;

		// intersect
		for (tipoint64 i = 0; i < smallM->num_trans; i++) {

			// look-ahead heuristic
			if (j + LOOK < big->num_trans &&
				big->trans[j + LOOK] < smallM->trans[i])
				j += LOOK;

			while (big->trans[j] < smallM->trans[i] && j < big->num_trans)
				j++;

			// found a match
			if (smallM->trans[i] == big->trans[j]) {

				// insert it into intersection set
				target->trans[target->num_trans] = smallM->trans[i];
				target->num_trans++;

				j++;
			}

			if (j == big->num_trans) break;
		}

		// if fails minsup, erase everything
		if (target->num_trans < minsup) {
			delete[] target->trans;
			target->trans = NULL;
			target->capacity = 0;
			target->num_trans = 0;
			return;
		}

		// shrink down transaction list if capacity is too big
		if (target->capacity - target->num_trans > SHRINK) {

			tipoint64 *old = target->trans;

			// resize and copy over the older stuff
			target->capacity = target->num_trans;
			target->trans = new tipoint64[target->num_trans];
			memcpy(target->trans, old, sizeof(tipoint64) * target->num_trans);

			// delete the older stuff
			delete[] old;
		}
	}

	void io_intersect(cell *target, cell *a, cell *b)
	{
		cell *smallM, *big;
		tipoint64 size = a->num_trans;

		smallM = b;
		big = a;

		// not enough for minsup
		if (size < minsup) {
			target->trans = NULL;
			target->num_trans = 0;
			target->capacity = 0;
			return;
		}

		target->trans = new tipoint64[size];
		target->num_trans = 0;
		target->capacity = size;

		tipoint64 j = 0;
		tipoint64 mod = -1;
		tipoint64 mod_other = -1;

		// intersect
		for (tipoint64 i = 0; i < smallM->num_trans; i++) {

			if (i / 1024 != mod_other) {
				_IOMAX++;
				mod_other = i / 1024;
			}

			// look-ahead heuristic
			if (j + LOOK < big->num_trans &&
				big->trans[j + LOOK] < smallM->trans[i])
				j += LOOK;

			while (big->trans[j] < smallM->trans[i] && j < big->num_trans)
				j++;

			// found a match
			if (smallM->trans[i] == big->trans[j]) {

				// insert it into intersection set
				target->trans[target->num_trans] = smallM->trans[i];
				target->num_trans++;

				j++;

				// I/O count
				if (i / 1024 != mod) {
					_IOCOUNT++;
					mod = i / 1024;
				}
			}

			if (j == big->num_trans) break;
		}

		// if fails minsup, erase everything
		if (target->num_trans < minsup) {
			delete[] target->trans;
			target->trans = NULL;
			target->capacity = 0;
			target->num_trans = 0;
			return;
		}

		// shrink down transaction list if capacity is too big
		if (target->capacity - target->num_trans > SHRINK) {

			tipoint64 *old = target->trans;

			// resize and copy over the older stuff
			target->capacity = target->num_trans;
			target->trans = new tipoint64[target->num_trans];
			memcpy(target->trans, old, sizeof(tipoint64) * target->num_trans);

			// delete the older stuff
			delete[] old;
		}
	}


	void multi_intersect(cell *target, cell **sources, tipoint64 num_sources)
	{
		tipoint64 *index = new tipoint64[num_sources];
		tipoint64 size = n_tuples + 1;
		tipoint64 smallest_source = -1;
		bool all_equal, done;

		for (tipoint64 i = 0; i < num_sources; i++) {
			if (sources[i] != NULL) {
				if (sources[i]->num_trans < size) {
					size = sources[i]->num_trans;
					smallest_source = i;
				}
			}
		}

		if (smallest_source == -1 || size < minsup) {
			target->trans = NULL;
			target->num_trans = 0;
			target->capacity = 0;
			delete[] index;
			return;
		}
		else {
			target->trans = new tipoint64[size];
			target->num_trans = 0;
			target->capacity = size;
		}

		memset(index, 0, sizeof(tipoint64) * num_sources);
		done = false;

		for (tipoint64 i = 0; i < size; i++) {

			for (tipoint64 j = 0; j < num_sources; j++) {
				if (sources[j] != NULL) {
					while (sources[j]->trans[index[j]] <
						sources[smallest_source]->trans[i] && index[j] <
						sources[j]->num_trans) {

						index[j]++;
					}
				}
			}

			all_equal = true;

			for (tipoint64 j = 0; j < num_sources; j++) {
				if (sources[j] != NULL && sources[j]->trans[index[j]] !=
					sources[smallest_source]->trans[i])
					all_equal = false;
			}

			if (all_equal) {
				target->trans[target->num_trans] =
					sources[smallest_source]->trans[i];
				target->num_trans++;

				for (tipoint64 j = 0; j < num_sources; j++) {
					if (sources[j] != NULL)
						index[j]++;
				}
			}

			for (tipoint64 j = 0; j < num_sources; j++) {
				if (sources[j] != NULL && index[j] == sources[j]->num_trans)
					done = true;
			}

			if (done)
				break;
		}

		// if fails minsup, erase everything
		if (target->num_trans < minsup) {
			delete[] target->trans;
			target->capacity = 0;
			target->num_trans = 0;
		}

		// shrink down transaction list if capacity is too big
		if (target->capacity - target->num_trans > SHRINK) {

			//cout << "shrinking by " << target->capacity - target->num_trans << endl;

			tipoint64 *old = target->trans;

			// resize
			target->capacity = target->num_trans;

			// copy over the older stuff
			target->trans = new tipoint64[target->num_trans];
			memcpy(target->trans, old, sizeof(tipoint64) * target->num_trans);

			// delete the older stuff
			delete[] old;
		}

		delete[] index;
	}


	cell* cell_of(tipoint64 *cell_name, tipoint64 name_size)
	{
		tipoint64 j, offset, temp;
		bool at_root = true;
		node *n;

		j = 0;

		while (cell_name[j] == -1)
			j++;

		n = root->child[frag[j]];

		j = 0;
		offset = 0;

		for (tipoint64 i = 0; i < n_dimensions; i++) {

			// traverse down the tree and locate the node
			if (cell_name[i] != -1) {

				if (at_root) {
					n = n->child[frag_index[i]];
					at_root = false;
				}
				else {
					if (n->child == NULL)
						return NULL;

					n = n->child[i - n->name - 1];
				}

				j++;
			}

			// found the node, now find the cell
			if (j == name_size) {

				// calculate the offset of the cell within this node
				for (tipoint64 u = 0; u < n_dimensions; u++) {

					if (cell_name[u] > 0) {

						temp = 1;

						for (tipoint64 v = u + 1; v < n_dimensions; v++) {

							if (cell_name[v] != -1) {
								temp *= cardinality[v];
							}
						}

						offset += temp * cell_name[u];
					}
				}

				return n->cells[offset];
			}
		}

		// should never get here
		cout << "Invalid cell_of()!" << endl;
		for (tipoint64 i = 0; i < n_dimensions; i++) {
			if (cell_name[i] == -1)
				cout << "* ";
			else
				cout << cell_name[i] << " ";
		}
		cout << "{" << name_size << "}" << endl;
		return NULL;
	}

	cell* online_cell_of(node *local_root, tipoint64 *cell_name, tipoint64 name_size,
		tipoint64 *local_name)
	{
		tipoint64 j, offset, temp;
		node *n = local_root;
		bool at_root = true;

		j = 0;
		offset = 0;

		// for all characters in the name
		for (tipoint64 i = 0; i < n_dimensions; i++) {

			// traverse down the tree if the character is non-empty
			if (cell_name[i] != -1) {

				// traverse down the tree
				if (at_root) {
					at_root = false;
					n = n->child[local_name[i]];
				}
				else {
					n = n->child[local_name[i] - local_name[n->name] - 1];
				}

				// record how many non-empty characters we've seen
				j++;
			}

			// found the node, now find the cell
			if (j == name_size) {

				// calculate the offset of the cell within this node
				for (tipoint64 u = 0; u < n_dimensions; u++) {

					if (cell_name[u] > 0) {

						temp = 1;

						for (tipoint64 v = u + 1; v < n_dimensions; v++) {

							if (cell_name[v] != -1) {
								temp *= cardinality[v];
							}
						}

						offset += temp * cell_name[u];
					}
				}

				return n->cells[offset];
			}
		}

		// should never get here
		cout << "Invalid online_cell_of()!" << endl;
		for (tipoint64 i = 0; i < n_dimensions; i++) {
			if (cell_name[i] == -1)
				cout << "* ";
			else
				cout << cell_name[i] << " ";
		}
		cout << "{" << name_size << "}" << endl;
		return NULL;
	}

	void elapsed_time(timeval *st)
	{
		timeval *et = new timeval();
		gettimeofdayM(et, NULL);
		// ISAIAS long time = 1000 * (et->tv_sec - st->tv_sec) + (et->tv_usec - st->tv_usec) / 1000;
		// ISAIAS cout << " ... used time: " << time << " ms." << endl;
		memcpy(st, et, sizeof(timeval));
		delete et;
	}

	void mark_time(timeval *st)
	{
		gettimeofdayM(st, NULL);
	}

	void time_taken(timeval *st)
	{
		timeval *et = new timeval();
		gettimeofdayM(et, NULL);
		long time = 1000 * (et->tv_sec - st->tv_sec) + (et->tv_usec -
			st->tv_usec) / 1000;
		cout << time << " ms." << endl;

		memcpy(st, et, sizeof(timeval));
		delete et;
	}

	long final_time_taken(timeval *st)
	{
		timeval *et = new timeval();
		gettimeofdayM(et, NULL);
		long time = 1000 * (et->tv_sec - st->tv_sec) + (et->tv_usec -
			st->tv_usec) / 1000;
		delete et;

		return time;
	}

	void exp_time_taken(timeval *st)
	{
		timeval *et = new timeval();
		gettimeofdayM(et, NULL);
		long time = 1000 * (et->tv_sec - st->tv_sec) + (et->tv_usec -
			st->tv_usec) / 1000;
		cout << time << " ms." << endl;
		memcpy(st, et, sizeof(timeval));
		delete et;
	}

	void interfaceM(string query, map<string, long>* mapa, long* tempo_execute_query)
	{
		string word;
		bool verbose, file_out;

		verbose = false;
		file_out = false;

		while (true) {
			// ISAIAS cout << "> ";
			// ISAIAS getline(cin, query);

			istringstream iss(query);

			iss >> word;

			/* ISAIAS if (word == "quit" || word == "exit" || word == "x")
				break;
			else if (word == "m" || word == "minsup") {
				iss >> minsup;
				cout << "Minimum support set to " << minsup << "." << endl;
			}
			else if (word == "v" || word == "verbose") {
				if (verbose) {
					verbose = false;
					cout << "Verbose mode off." << endl;
				}
				else {
					verbose = true;
					cout << "Verbose mode on." << endl;
				}
			}
			// file output?
			else if (word == "f" || word == "file") {
				if (file_out) {
					file_out = false;
					cout << "File output off." << endl;
				}
				else {
					file_out = true;
					cout << "File output on." << endl;
				}
			}
			// query
			else ISAIAS */
			if (word == "q" || word == "query") {
				timeval* inicio_execute = new timeval();
				gettimeofdayM(inicio_execute, NULL);

				mark_time(prev_time);

				num_rows = 0;
				execute_query(query, verbose, file_out, mapa);

				cout << endl << "    " << num_rows << " rows returned.";
				cout << "  Query executed in ";
				time_taken(prev_time);

				*tempo_execute_query = final_time_taken(inicio_execute);
				delete inicio_execute;

			}
			// experiments
			/* ISAIAS else if (word == "exp") {
				tipoint64 *query_value = new tipoint64[NUM_DIMS];
				tipoint64 num_ins, num_inq;
				tipoint64 r;
				tipoint64 value;
				string que;
				tipoint64 total_rows = 0;

				srand(1);

				mark_time(prev_time);
				_IOCOUNT = 0;
				_IOMAX = 0;
				_CUBEIO = 0;

				// do NUM_QUERIES queries
				for (tipoint64 i = 0; i < NUM_QUERIES; i++) {

					num_ins = 0;
					num_inq = 0;

					// all irrelevant dimensions to begin with
					for (tipoint64 i = 0; i < NUM_DIMS; i++)
						query_value[i] = -2;

					// mark the instantiated dimensions
					while (num_ins < NUM_INS) {

						r = (tipoint64)((double)NUM_DIMS * rand() / RAND_MAX);
						if (r == NUM_DIMS) r = 0;

						if (query_value[r] == -2) {
							value = (tipoint64)((double)CARD * rand() / RAND_MAX);
							if (value == 0) value = 10;
							query_value[r] = value;
							num_ins++;
						}
					}

					// mark the inquired dimensions
					while (num_inq < NUM_INQ) {

						r = (tipoint64)((double)NUM_DIMS * rand() / RAND_MAX);
						if (r == NUM_DIMS) r = 0;

						if (query_value[r] == -2) {
							query_value[r] = -1;
							num_inq++;
						}
					}

					tipoint64 _CUBEIO_TEMP = 951483 * 8 / 1024;
					bool stop = false;

					// make the query string
					que = "q";

					// make the query string
					for (tipoint64 i = 0; i < NUM_DIMS; i++) {

						// irrelevant
						if (query_value[i] == -2) {
							que += " *";

							if (!stop)
								_CUBEIO_TEMP /= 10;
						}

						// query
						else if (query_value[i] == -1) {
							que += " ?";
							stop = true;
						}

						// instantiated
						else {
							que += " " + itos(query_value[i]);

							if (!stop)
								_CUBEIO_TEMP /= 10;
						}
					}

					_CUBEIO += _CUBEIO_TEMP;

					//cout << "  " << que << endl;
					num_rows = 0;
					execute_query(que, false, true, mapa);
					//cout << "    " << num_rows << " rows returned." << endl << endl;

					total_rows += num_rows;
				}

				cout << (double)total_rows / (double)NUM_QUERIES << " rows per query." << endl;

				cout << NUM_QUERIES << " query executed in ";
				exp_time_taken(prev_time);
				cout << _IOMAX << " I/Os counted." << endl;
				cout << _CUBEIO << " cube I/Os counted." << endl;
				cout << endl;
			}
			else if (word == "");
			else {
				cout << "Invalid command." << endl;
			} ISAIAS */

			word = "";

			// ISAIAS - Saida forçada do while
			break;
		}
	}

	// converts a string to an integer 
	tipoint64 stoiM(const string &s)
	{
		tipoint64 result;
		istringstream(s) >> result;
		return result;
	}


	void execute_query(string query, bool verbose, bool file_out, map<string, long>* mapa)
	{
		string word;
		tipoint64 n_dims, *query_cell;
		bool empty_result;

		istringstream iss(query);

		// the first word is presumably "q" or "query"
		iss >> word;

		n_dims = 0;
		query_cell = new tipoint64[n_dimensions];

		// read in the query values
		while (iss >> word) {
			if (word == "*")
				query_cell[n_dims] = -1;
			else if (word == "?")
				query_cell[n_dims] = -2;
			else
				query_cell[n_dims] = stoiM(word);

			n_dims++;
		}

		// not the correct number of dimensions
		if (n_dims != n_dimensions) {
			cout << "    ";
			cout << "Invalid query of " << n_dims << " dimensions.  Loaded "
				<< "database has " << n_dimensions << " dimensions." << endl;
			delete[] query_cell;
			return;
		}

		// re-iterate the query
		if (file_out) {
			//outfile << endl << "    Executing query: ";
			outfile << "Executing query: ";
			for (tipoint64 i = 0; i < n_dimensions; i++) {
				if (query_cell[i] == -1)
					outfile << "* ";
				else if (query_cell[i] == -2)
					outfile << "? ";
				else
					outfile << query_cell[i] << " ";
			}
			//outfile << endl << endl;
			outfile << endl;
		}
		else {
			cout << endl << "    Executing query: ";
			for (tipoint64 i = 0; i < n_dimensions; i++) {
				if (query_cell[i] == -1)
					cout << "* ";
				else if (query_cell[i] == -2)
					cout << "? ";
				else
					cout << query_cell[i] << " ";
			}
			cout << endl << endl;
		}

		empty_result = false;

		// make sure the query is not asking for a value that is too big or
		// too small
		for (tipoint64 i = 0; i < n_dimensions; i++) {
			if (query_cell[i] >= cardinality[i] || query_cell[i] < -2) {
				cout << "Invalid query value at dimension " << i << "." << endl;
				empty_result = true;
			}
		}

		// nothing
		if (empty_result) {
			delete[] query_cell;
			return;
		}

		//===========================================//
		//  the actual query processing starts here  //
		//===========================================//

		tipoint64 *cell_name, cell_size, num_variables, num_instantiated;
		cell **group_cells, *fixed_results;

		cell_name = new tipoint64[n_dimensions];
		group_cells = new cell*[n_frags];
		fixed_results = new cell;
		fixed_results->num_trans = 0;
		num_variables = 0;
		num_instantiated = 0;

		// count the number of inquired and instantiated dimensions
		for (tipoint64 i = 0; i < n_dimensions; i++) {
			if (query_cell[i] == -2) num_variables++;
			if (query_cell[i] >= 0) num_instantiated++;
		}

		// instantiated dimensions exist
		if (num_instantiated > 0) {

			// locate the fixed shell fragments
			for (tipoint64 i = 0; i < n_frags; i++) {
				memset(cell_name, -1, sizeof(tipoint64) * n_dimensions);
				cell_size = 0;

				// for all the queries in this fragment
				for (tipoint64 j = 0; j < root->child[i]->num_cells; j++) {

					// if this fragment has a instantiated value, mark it
					if (query_cell[root->child[i]->child[j]->name] >= 0) {
						cell_name[root->child[i]->child[j]->name] =
							query_cell[root->child[i]->child[j]->name];
						cell_size++;
					}
				}

				// if this fragment has a instantiated value, retrieve the
				// tidlist
				if (cell_size > 0) {

					// grab its transaction list from the pre-computed
					// shell-fragments
					group_cells[i] = cell_of(cell_name, cell_size);

					if (ONLINE_DEBUG) {
						cout << "instantiated cell: ";
						for (tipoint64 u = 0; u < n_dimensions; u++) {
							if (cell_name[u] == -1)
								cout << "* ";
							else
								cout << cell_name[u] << " ";
						}
						cout << "[" << group_cells[i]->num_trans << "]";
						cout << " -> " << group_cells[i]->num_trans / 1024.0
							<< " I/O" << endl;
					}

					// count I/O of fetching the tidlist
					double io = group_cells[i]->num_trans / 1024.0;
					_IOMAX += (tipoint64)ceil(io);
				}
				else
					// this fragment has no instantiated values
					group_cells[i] = NULL;
			}

			// intersect the instantiated shell fragments
			multi_intersect(fixed_results, group_cells, n_frags);

			if (ONLINE_DEBUG) {
				cout << "fixed multi-intersect result: ";
				//for (int u = 0; u < fixed_results->num_trans; u++) {
				//	cout << fixed_results->trans[u] << " ";
				//}
				cout << "[" << fixed_results->num_trans << "]" << endl << endl;
			}
		}

		// the original query was not the entire cube and it didn't have any
		// transactions.  basically, no results.
		if (fixed_results->num_trans == 0 && num_variables < n_dimensions
			&& num_instantiated > 0) {
			delete[] query_cell;
			delete[] cell_name;
			delete[] group_cells;
			delete fixed_results;
			return;
		}

		// point query, just output the results here.  no need for further
		// local cubing.
		if (num_variables == 0 && fixed_results->num_trans > 0
			&& num_instantiated > 0) {
			if (file_out) {
				/*
				outfile << "    ";
				for (int i = 0; i < n_dimensions; i++) {
				if (query_cell[i] == -1)
				outfile << "* ";
				else if (query_cell[i] == -2)
				outfile << "? ";
				else
				outfile << query_cell[i] << " ";
				}
				outfile << ": " << fixed_results->num_trans << endl;
				*/
			}
			else {
				cout << "    ";

				string label = "";
				for (tipoint64 i = 0; i < n_dimensions; i++) {
					if (query_cell[i] == -1)
						label += "*";
					else if (query_cell[i] == -2)
						label += "?";
					else
						label += to_string(query_cell[i]);
					if (i + 1 < n_dimensions)
						label += " ";
				}

				cout << label << " : " << fixed_results->num_trans << endl;
				(*mapa)[label] = fixed_results->num_trans;
			
			}

			num_rows++;

			delete[] query_cell;
			delete[] cell_name;
			delete[] group_cells;
			delete fixed_results;
			return;
		}


		//===============================================================//
		// fixed results are non-empty and the number of variables is at
		// least 1, we have a local cube to mine. 
		//===============================================================//

		node *local_root;
		cell *c, *temp;
		tipoint64 *local_name, *local_order, index;
		vector< tipoint64 > cuboid_name;

		// create local cube's root
		local_root = new node;
		local_root->name = -500;
		local_root->child = new node*[num_variables];
		local_root->num_cells = num_variables;
		local_root->cells = new cell*[1];
		local_root->cells[0] = fixed_results;

		local_name = new tipoint64[n_dimensions];
		local_order = new tipoint64[num_variables];

		// do the base cuboid if computing full cube
		if (num_variables == n_dimensions || num_instantiated == 0) {
			local_root->cells[0]->num_trans = n_tuples;
		}

		num_rows++;

		// output the base cuboid
		if (file_out) {
			/*
			outfile << "    ";
			for (int i = 0; i < n_dimensions; i++) {
			if (query_cell[i] == -1)
			outfile << "* ";
			else if (query_cell[i] == -2)
			outfile << "* ";
			else
			outfile << query_cell[i] << " ";
			}
			outfile << ": " << local_root->cells[0]->num_trans << endl;
			*/
		}
		else {
			cout << "    ";
			string label = "";
			for (tipoint64 i = 0; i < n_dimensions; i++) {
				if (query_cell[i] == -1)
					label += "*";
				else if (query_cell[i] == -2)
					label += "*";
				else
					label += to_string(query_cell[i]);
				if (i + 1 < n_dimensions)
					label += " ";
			}

			cout << label << " : " << local_root->cells[0]->num_trans << endl;
			(*mapa)[label] = local_root->cells[0]->num_trans;

		}

		index = 0;
		memset(local_name, -1, sizeof(tipoint64) * n_dimensions);

		//_IOCOUNT = 0;
		//_IOMAX = 0;

		// create the local cube's 1-dimensional cuboids
		for (tipoint64 i = 0; i < n_dimensions; i++) {

			// found a variable, make its node under local root
			if (query_cell[i] == -2) {

				// local name of dimension i is index
				local_name[i] = index;

				// the index_th local dimension is i
				local_order[index] = i;

				local_root->child[index] = new node;
				local_root->child[index]->name = i;
				local_root->child[index]->num_cells = cardinality[i];
				local_root->child[index]->cells = new cell*[cardinality[i]];
				local_root->child[index]->cells[0] = NULL;
				local_root->child[index]->child = NULL;

				// populate the cells for this node
				for (tipoint64 k = 1; k < cardinality[i]; k++) {
					local_root->child[index]->cells[k] = new cell;
					c = local_root->child[index]->cells[k];
					temp = root->child[frag[i]]->child[frag_index[i]]->cells[k];

					// in case we're computing the entire cube or there are
					// no instantiated dimensions, just copy over the
					// original transaction list
					if (num_variables == n_dimensions ||
						num_instantiated == 0) {
						c->capacity = temp->num_trans;
						c->num_trans = c->capacity;
						c->trans = new tipoint64[c->capacity];
						memcpy(c->trans, temp->trans, sizeof(tipoint64) * c->capacity);
					}
					// intersect with local base cuboid
					else {
						// have to count I/O here
						io_intersect(c, local_root->cells[0], temp);
					}
				}

				index++;
			}
		}

		/*
		cout << "I/O Count = " << _IOCOUNT << endl << endl;
		cout << "I/O Max = " << _IOMAX << endl << endl;
		*/

		if (ONLINE_DEBUG) {
			cout << "local_name: ";
			for (tipoint64 u = 0; u < n_dimensions; u++) {
				if (local_name[u] == -1)
					cout << "* ";
				else
					cout << local_name[u] << " ";
			}
			cout << endl;
			cout << "local_order: ";
			for (tipoint64 u = 0; u < num_variables; u++) {
				cout << local_order[u] << " ";
			}
			cout << endl;
		}

		memset(cell_name, -1, sizeof(tipoint64) * n_dimensions);
		cuboid_name.reserve(num_variables);

		for (tipoint64 i = 0; i < n_dimensions; i++) {
			if (query_cell[i] >= 0)
				cell_name[i] = query_cell[i];
		}

		// do the local data cube
		for (tipoint64 i = 0; i < local_root->num_cells; i++) {
			cuboid_name.push_back(local_root->child[i]->name);
			online_fragment(local_root, local_root->child[i], NULL,
				cuboid_name, cell_name, num_variables, local_name,
				local_order, file_out, mapa);
			cuboid_name.pop_back();
		}

		// delete local cube
		// free_tree(local_root, local_name, n_dims);

		delete[] query_cell;
		delete[] cell_name;
		delete[] group_cells;
		delete[] local_name;
		delete[] local_order;

		return;
	}



	void online_fragment(node *local_root, node *my_root, cell *parent_cell,
		vector< tipoint64 > cuboid_name, tipoint64 *cell_name, tipoint64 n_dim, tipoint64
		*local_name, tipoint64 *local_order, bool file_out, map<string, long>* mapa)
	{
		tipoint64 num_children, cell_offset, temp, last_name;
		cell *c;
		bool meet_minsup = true;

		if (ONLINE_DEBUG) {
			cout << "online_fragment starting at" << endl;
			cout << "\tcuboid: ";
			for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++)
				cout << cuboid_name[i] << " ";
			cout << endl;
			cout << "\tcell  : ";
			for (tipoint64 i = 0; i < n_dimensions; i++) {
				if (cell_name[i] == -1)
					cout << "* ";
				else
					cout << cell_name[i] << " ";
			}
			cout << endl;
		}

		cell_offset = 0;
		last_name = cuboid_name[cuboid_name.size() - 1];

		if (parent_cell != NULL) {

			// allocate my cells if it hasn't been done yet.
			if (my_root->num_cells <= 0) {
				my_root->num_cells = 1;

				for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++)
					my_root->num_cells *= cardinality[cuboid_name[i]];

				my_root->cells = new cell*[my_root->num_cells];

				for (tipoint64 i = 0; i < my_root->num_cells; i++)
					my_root->cells[i] = NULL;

				if (ONLINE_DEBUG)
					cout << "created " << my_root->num_cells << " cells." << endl;
			}

			// compute cell offset within this node
			for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++) {

				if (cell_name[cuboid_name[i]] != -1) {

					temp = 1;
					for (unsigned tipoint64 j = i + 1; j < cuboid_name.size(); j++)
						temp *= cardinality[cuboid_name[j]];

					cell_offset += temp * cell_name[cuboid_name[i]];
				}
			}

			if (ONLINE_DEBUG) {
				cout << "cell_offset = " << cell_offset << endl;
				cout << "last_name = " << last_name << endl;
			}


			tipoint64 same_frag, *same_frag_cell_name, j, cell_size;
			cell *same_cell, *unsame_cell;

			same_cell = NULL;
			unsame_cell = NULL;
			cell_size = 0;
			same_frag = 0;
			same_frag_cell_name = new tipoint64[n_dimensions];

			// if we're computing the full data cube, we should use the
			// shell fragments to our advantage.
			if (n_dim == n_dimensions) {

				memset(same_frag_cell_name, -1, sizeof(tipoint64) * n_dimensions);

				for (unsigned tipoint64 i = 0; i < cuboid_name.size(); i++) {
					if (cell_name[cuboid_name[i]] != -1) {
						same_frag_cell_name[cuboid_name[i]] =
							cell_name[cuboid_name[i]];
						cell_size++;
					}
				}

				j = cuboid_name.size() - 2;

				while (j >= 0 && frag[last_name] == frag[cuboid_name[j]]) {
					same_frag_cell_name[cuboid_name[j]] = -1;
					j--;
					cell_size--;
				}

				same_frag = cuboid_name.size() - j - 2;

				if (ONLINE_DEBUG) {
					cout << "same_frag = " << same_frag << endl;
				}

				// locate the un-same cell
				if (same_frag > 0) {

					// this cell and the immediately preceding cells in the
					// local cube have already been aggregated in the shell
					// fragments.  fetch their stored answer and intersect
					// with the rest.

					if ((tipoint64)cuboid_name.size() - same_frag - 1 > 0) {

						if (ONLINE_DEBUG) {
							cout << "unsame cell: ";
							for (tipoint64 u = 0; u < n_dimensions; u++) {
								if (same_frag_cell_name[u] == -1)
									cout << "* ";
								else
									cout << same_frag_cell_name[u] << " ";
							}
							cout << ": name_size = " << cell_size << "; ";
						}

						unsame_cell = online_cell_of(local_root,
							same_frag_cell_name, cell_size, local_name);

						if (ONLINE_DEBUG) {
							cout << "[" << unsame_cell->num_trans << "]" <<
								endl;
						}
					}

				} // end same frag
			}

			meet_minsup = false;

			memset(same_frag_cell_name, -1, sizeof(tipoint64) * n_dimensions);

			j = cuboid_name.size() - 2;

			while (j >= 0 && frag[last_name] == frag[cuboid_name[j]]) {
				same_frag_cell_name[cuboid_name[j]] =
					cell_name[cuboid_name[j]];
				j--;
			}

			// aggregate all my cells with my parent cell
			for (tipoint64 i = cell_offset + 1; i < cardinality[last_name] +
				cell_offset; i++) {

				my_root->cells[i] = new cell;

				cell_name[last_name] = i - cell_offset;

				if (n_dim == n_dimensions && same_frag > 0) {

					same_frag_cell_name[last_name] = cell_name[last_name];

					same_cell = cell_of(same_frag_cell_name, same_frag + 1);

					if (ONLINE_DEBUG) {
						cout << "same cell: ";
						for (tipoint64 u = 0; u < n_dimensions; u++) {
							if (same_frag_cell_name[u] == -1)
								cout << "* ";
							else
								cout << same_frag_cell_name[u] << " ";
						}
						cout << "[" << same_cell->num_trans << "]" << endl;
					}

					if (unsame_cell == NULL) {
						if (ONLINE_DEBUG) {
							cout << "unsame_cell == NULL" << endl;
						}
						if (local_root->cells[0]->num_trans == n_tuples) {
							my_root->cells[i]->num_trans = same_cell->num_trans;
							my_root->cells[i]->capacity = same_cell->num_trans;
							my_root->cells[i]->trans = new tipoint64[same_cell->num_trans];
							memcpy(my_root->cells[i]->trans,
								same_cell->trans, sizeof(tipoint64) *
								same_cell->num_trans);
						}
						else {
							if (ONLINE_DEBUG) {
								cout << "intersecting local root data" << endl;
							}
							intersect(my_root->cells[i], local_root->cells[0], same_cell);
						}
					}
					else {
						intersect(my_root->cells[i], unsame_cell, same_cell);
					}
				}
				else {

					// this cell doesn't share any fragments with the
					// immediately preceding cells in the local cube, thus
					// just do the plain intersection
					c = local_root->child[local_name[last_name]]->cells[i - cell_offset];

					intersect(my_root->cells[i], parent_cell, c);
				}

				if (my_root->cells[i]->num_trans >= minsup)
					meet_minsup = true;

				/*
				if (ONLINE_DEBUG) {
				cout << "result: ";
				cout << "[" << my_root->cells[i]->num_trans << "]: ";
				for (int j = 0; j < my_root->cells[i]->num_trans; j++)
				cout << my_root->cells[i]->trans[j] << " ";
				cout << endl;
				}
				*/
			}

			// reset name
			cell_name[last_name] = -1;
			delete[] same_frag_cell_name;
		}


		if (!meet_minsup) {
			return;
		}

		if (ONLINE_DEBUG)
			cout << "n_dim = " << n_dim << endl;

		num_children = 0;

		// not at a leaf yet, need to create my children first
		if (local_name[last_name] < n_dim - 1) {

			num_children = n_dim - 1 - local_name[my_root->name];

			if (ONLINE_DEBUG)
				cout << "num_children = " << num_children << endl;

			if (my_root->child == NULL) {

				// create my children nodes
				my_root->child = new node*[num_children];

				for (tipoint64 i = 0; i < num_children; i++) {
					my_root->child[i] = new node;
					my_root->child[i]->name = local_order[local_name[my_root->name] + i + 1];
					my_root->child[i]->num_cells = 0;
					my_root->child[i]->child = NULL;
				}

				if (ONLINE_DEBUG)
					cout << "created " << num_children << " children." << endl;
			}
		}

		// for all my cells
		for (tipoint64 i = cell_offset + 1; i < cardinality[last_name] +
			cell_offset; i++) {

			c = my_root->cells[i];

			// check iceberg condition
			if (c != NULL && c->num_trans >= minsup) {

				cell_name[last_name] = i - cell_offset;

				if (ONLINE_DEBUG) {
					cout << "cell ";
					for (tipoint64 u = 0; u < n_dimensions; u++) {
						if (cell_name[u] == -1)
							cout << "* ";
						else
							cout << cell_name[u] << " ";
					}
					cout << endl;
				}

				// output it
				if (file_out) {
					/*
					outfile << "    ";
					for (int u = 0; u < n_dimensions; u++) {
					if (cell_name[u] == -1)
					outfile << "* ";
					else
					outfile << cell_name[u] << " ";
					}
					outfile << ": " << c->num_trans << endl;
					*/
				}
				else {
					cout << "    ";
					string label = "";
					for (tipoint64 u = 0; u < n_dimensions; u++) {
						if (cell_name[u] == -1)
							//cout << "* ";
							label += "*";
						else
							//cout << cell_name[u] << " ";
							label += to_string(cell_name[u]);
						if (u + 1 < n_dimensions)
							label += " ";
					}

					cout << label << " : " << c->num_trans << endl;
					(*mapa)[label] = c->num_trans;

				}

				num_rows++;

				// recursive call for all my children (if any)
				for (tipoint64 j = 0; j < num_children; j++) {
					cuboid_name.push_back(my_root->child[j]->name);
					online_fragment(local_root, my_root->child[j], c,
						cuboid_name, cell_name, n_dim, local_name,
						local_order, file_out, mapa);
					cuboid_name.pop_back();
				}
			}

			cell_name[last_name] = -1;
		}
	}

	void free_tree(node *n, tipoint64 *local_name, tipoint64 n_dim)
	{
		if (n == NULL)
			return;

		tipoint64 num_children = n_dim - 1 - local_name[n->name];

		if (n->name == -500) num_children = n_dim;

		for (tipoint64 i = 0; i < num_children; i++) {
			free_tree(n->child[i], local_name, n_dim);
		}

		if (num_children == n_dim) {
			if (n->cells[0] != NULL && n->cells[0]->capacity > 0) {
				delete[] n->cells[0]->trans;
				delete n->cells[0];
			}
		}
		else {
			for (tipoint64 j = 0; j < n->num_cells; j++) {
				if (n->cells[j] != NULL && n->cells[j]->num_trans > 0) {
					delete[] n->cells[j]->trans;
					delete n->cells[j];
				}
			}
		}

		if (n->num_cells > 0)
			delete[] n->cells;

		if (num_children > 0)
			delete[] n->child;

		delete n;
	}

	string itos(tipoint64 n)
	{
		ostringstream o;
		o << n;
		return o.str();
	}

	tipoint64 gettimeofdayM(struct timeval *tv, struct timezone2 *tz)
	{
		FILETIME ft;
		unsigned tipoint64 tmpres = 0;
		static tipoint64 tzflag;

		if (NULL != tv)
		{
			GetSystemTimeAsFileTime(&ft);

			tmpres |= ft.dwHighDateTime;
			tmpres <<= 32;
			tmpres |= ft.dwLowDateTime;

			/*converting file time to unix epoch*/
			tmpres -= DELTA_EPOCH_IN_MICROSECS;
			tmpres /= 10;  /*convert into microseconds*/
			tv->tv_sec = (long)(tmpres / 1000000UL);
			tv->tv_usec = (long)(tmpres % 1000000UL);
		}

		if (NULL != tz)
		{
			if (!tzflag)
			{
				_tzset();
				tzflag++;
			}
			tz->tz_minuteswest = _timezone / 60;
			tz->tz_dsttime = _daylight;
		}

		return 0;
	}

#endif
