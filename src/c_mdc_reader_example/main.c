#include "mdc_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct u32_vec
{
	uint32_t* data;
	uint32_t size;
	uint32_t capacity;
} u32_vec;

u32_vec u32_vec_create(void)
{
	u32_vec vec;
	vec.data = NULL;
	vec.size = 0;
	vec.capacity = 0;
	return vec;
}

void u32_vec_push_back(u32_vec* vec, uint32_t value)
{
	if (vec->size == vec->capacity)
	{
		vec->capacity = vec->capacity == 0 ? 1 : vec->capacity * 2;
		vec->data = (uint32_t*)realloc(vec->data, vec->capacity * sizeof(uint32_t));
	}
	vec->data[vec->size++] = value;
}

void u32_vec_free(u32_vec* vec)
{
	free(vec->data);
	vec->size = 0;
	vec->capacity = 0;
	vec->data = NULL;
}

typedef struct str_vec
{
	char** data;
	uint32_t size;
	uint32_t capacity;
} str_vec;

str_vec str_vec_create(void)
{
	str_vec vec;
	vec.data = NULL;
	vec.size = 0;
	vec.capacity = 0;
	return vec;
}

void str_vec_push_back(str_vec* vec, char* value)
{
	if (vec->size == vec->capacity)
	{
		vec->capacity = vec->capacity == 0 ? 1 : vec->capacity * 2;
		vec->data = (char**)realloc(vec->data, vec->capacity * sizeof(char*));
	}
	vec->data[vec->size++] = value;
}

void str_vec_free(str_vec* vec)
{
	for (uint32_t i = 0; i < vec->size; i++)
	{
		free(vec->data[i]);
	}
	free(vec->data);
	vec->size = 0;
	vec->capacity = 0;
	vec->data = NULL;
}

struct parameters
{
	const char* input_file;
	const char* output_file;
	str_vec segment_ids;
	u32_vec atom_ids;
};

struct parameters params_create(void)
{
	struct parameters res;
	res.input_file = NULL;
	res.output_file = NULL;
	res.segment_ids = str_vec_create();
	res.atom_ids = u32_vec_create();

	return res;
}

void free_params(struct parameters* params)
{
	str_vec_free(&params->segment_ids);
	u32_vec_free(&params->atom_ids);
}

str_vec split_string(const char* str, char delim)
{
	str_vec vec = str_vec_create();
	const char* start = str;
	const char* end = str;
	while (*end != '\0')
	{
		if (*end == delim)
		{
			char* token = (char*)malloc(end - start + 1);
			memcpy(token, start, end - start);
			token[end - start] = '\0';
			str_vec_push_back(&vec, token);
			start = end + 1;
		}
		end++;
	}
	if (end != start)
	{
		char* token = (char*)malloc(end - start + 1);
		memcpy(token, start, end - start);
		token[end - start] = '\0';
		str_vec_push_back(&vec, token);
	}
	return vec;
}

u32_vec parse_atom_ids(char* str)
{
	u32_vec vec = u32_vec_create();
	str_vec tokens = split_string(str, ',');
	for (uint32_t i = 0; i < tokens.size; i++)
	{
		str_vec range = split_string(tokens.data[i], '-');
		if (range.size == 1)
		{
			u32_vec_push_back(&vec, atoi(range.data[0]));
		}
		else if (range.size == 2)
		{
			uint32_t start = atoi(range.data[0]);
			uint32_t end = atoi(range.data[1]);
			if (end < start)
			{
				printf("Error: Invalid atom range: %s\n", tokens.data[i]);
				str_vec_free(&range);
				u32_vec_free(&vec);
				break;
			}
			for (uint32_t j = start; j <= end && j < UINT32_MAX ; j++)
			{
				u32_vec_push_back(&vec, j);
			}
		}
		else
		{
			printf("Error: Invalid atom range: %s\n", tokens.data[i]);
			str_vec_free(&range);
			u32_vec_free(&vec);
			break;
		}
		str_vec_free(&range);
	}
	str_vec_free(&tokens);
	return vec;
}

void help(char* prog)
{
	printf("Usage: %s -i <input> -o <output_file>  [--segments seg1,seg2,...] [--atoms id1,start_range1-end_range1,...]\n", prog);
}

//here are the parameters
// -i input file
// -o output file
int parse_parameters(int argc, char** argv, struct parameters* params)
{
	params->input_file = NULL;
	params->output_file = NULL;

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-i") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("Error: -i requires an argument\n");
				help(argv[0]);
				return 1;
			}
			params->input_file = argv[i + 1];
			i++;
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("Error: -o requires an argument\n");
				help(argv[0]);
				return 1;
			}
			params->output_file = argv[i + 1];
			i++;
		}
		else if (strcmp(argv[i], "--segments") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("Error: --segments requires an argument\n");
				help(argv[0]);
				return 1;
			}
			params->segment_ids = split_string(argv[i + 1], ',');
		}
		else if (strcmp(argv[i], "--atoms") == 0)
		{
			if (i + 1 >= argc)
			{
				printf("Error: --atoms requires an argument\n");
				help(argv[0]);
				return 1;
			}
			params->atom_ids = parse_atom_ids(argv[i + 1]);
			
		}
		else
		{
			printf("Error: Unknown argument: %s\n", argv[i]);
			help(argv[0]);
			return 1;
		}
	}

	if (params->input_file == NULL)
	{
		printf("Error: -i is required\n");
		help(argv[0]);
		return 1;
	}

	if (params->output_file == NULL)
	{
		printf("Error: -o is required\n");
		help(argv[0]);
		return 1;
	}

	return 0;
}

//currently it just writes binary data to the output file
//not to xtc, but it is possible to convert xtc to binary with to_fp32 of main application
int main(int argc, char** argv)
{
	int ret = EXIT_FAILURE;
	struct parameters params = params_create();
	mdc_reader* reader = NULL;
	mdc_query_engine* engine = NULL;
	FILE* output = NULL;
	mdc_query_result* result = NULL;
	if (parse_parameters(argc, argv, &params) != 0)
	{
		goto NO_RAII_IN_C;
	}
	reader = mdc_reader_open(params.input_file);
	if (mdc_reader_get_error(reader))
	{
		printf("Error: %s\n", mdc_reader_get_error(reader));
		goto NO_RAII_IN_C;
	}

	engine = mdc_get_query_engine(reader, params.segment_ids.data, params.segment_ids.size, params.atom_ids.data, params.atom_ids.size);
	if (mdc_query_engine_get_error(engine))
	{
		printf("Error: %s\n", mdc_query_engine_get_error(engine));
		goto NO_RAII_IN_C;
	}

//	uint32_t no_segments = mdc_get_no_segments(reader);
	uint32_t no_frames = mdc_get_no_frames(reader);

	output = fopen(params.output_file, "wb");
	if (!output)
	{
		printf("Error: Cannot open output file\n");
		goto NO_RAII_IN_C;
	}

	result = mdc_create_query_result();

	for (uint32_t frame_id = 0; frame_id < no_frames; frame_id++)
	{
		int tot_atoms = 0;
		
		if (mdc_query(engine, &frame_id, 1, result))
		{
			printf("Error: %s\n", mdc_query_engine_get_error(engine));
			goto NO_RAII_IN_C;
		}
		assert(result->n_frames == 1);

		uint32_t n_atoms = result->frames[0].n_coords;
		tot_atoms += n_atoms;
		mdc_atom_coords* coords = result->frames[0].coords;
		fwrite(coords, sizeof(mdc_atom_coords), n_atoms, output);
		//for (uint32_t atom_id = 0; atom_id < n_atoms; atom_id++)
		//{
		//	fwrite(&coords[atom_id].x, sizeof(float), 1, output);
		//	fwrite(&coords[atom_id].y, sizeof(float), 1, output);
		//	fwrite(&coords[atom_id].z, sizeof(float), 1, output);
		//}
		
		printf("Frame time: %g\tNo. atoms: %d\n", result->frames[0].time, tot_atoms);
	}

	ret = EXIT_SUCCESS;

NO_RAII_IN_C:
	if (output)
		fclose(output);

	free_params(&params);
	mdc_free_query_result(result);
	mdc_free_query_engine(engine);
	mdc_reader_close(reader);

	return ret;
}