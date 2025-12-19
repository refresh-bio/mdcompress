#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "mdc_reader.h"

using std::vector;
using std::string;

struct Params
{
	string input_file;
	string output_file;
	vector<string> segment_ids;
	vector<uint32_t> atom_ids;
};

void help(char* prog)
{
	printf("Usage: %s -i <input> -o <output_file>  [--segments seg1,seg2,...] [--atoms id1,start_range1-end_range1,...]\n", prog);
}

std::vector<std::string> split_string(const std::string& s, char sep)
{
	auto p = s.begin();
	std::vector<std::string> components;

	while (true)
	{
		auto q = std::find(p, s.end(), sep);
		components.push_back(std::string(p, q));
		if (q == s.end())
			break;

		p = q + 1;
	}

	return components;
}

bool parse_args(int argc, char* argv[], Params& params)
{
	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]) == "-i")
		{
			if (i + 1 < argc)
			{
				params.input_file = argv[i + 1];
				i++;
			}
			else
			{
				std::cout << "Error: -i option requires one argument" << std::endl;
				help(argv[0]);
				return false;
			}
		}
		else if (string(argv[i]) == "-o")
		{
			if (i + 1 < argc)
			{
				params.output_file = argv[i + 1];
				i++;
			}
			else
			{
				std::cout << "Error: -o option requires one argument" << std::endl;
				help(argv[0]);
				return false;
			}
		}
		else if (string(argv[i]) == "--segments")
		{
			if (i + 1 < argc)
			{
				params.segment_ids = split_string(argv[i + 1], ',');
				i++;
			}
			else
			{
				std::cout << "Error: --segments option requires one argument" << std::endl;
				help(argv[0]);
				return false;
			}
		}
		else if (string(argv[i]) == "--atoms")
		{
			if (i + 1 < argc)
			{
				auto splitted = split_string(argv[++i], ',');
				for (auto& s : splitted)
				{
					auto dash = s.find('-');
					if (dash != std::string::npos)
					{
						int start = atoi(s.substr(0, dash).c_str());
						int end = atoi(s.substr(dash + 1).c_str());
						for (int j = start; j <= end; ++j)
							params.atom_ids.push_back(j);
					}
					else
					{
						params.atom_ids.push_back(atoi(s.c_str()));
					}
				}
			}
			else
			{
				std::cout << "Error: --atoms option requires one argument" << std::endl;
				help(argv[0]);
				return false;
			}
		}
		else
		{
			std::cout << "Error: unknown option " << argv[i] << std::endl;
			help(argv[0]);
			return false;
		}
	}

	if (params.input_file.empty())
	{
		printf("Error: -i is required\n");
		help(argv[0]);
		return false;
	}

	if (params.output_file.empty())
	{
		printf("Error: -o is required\n");
		help(argv[0]);
		return false;
	}
	return true;
}

string segment_type_to_str(mdc::segment_type type)
{
	switch (type)
	{
		case mdc::segment_type::unknown: return "unknown";
		case mdc::segment_type::molecule: return "molecule";
		case mdc::segment_type::other: return "other";
		case mdc::segment_type::water: return "water";
		case mdc::segment_type::none: return "none";
	}
	return ""; // this should never happen
}

int main(int argc, char* argv[])
{
	Params params;
	if (!parse_args(argc, argv, params))
		return 1;

	mdc::reader reader(params.input_file);
	if (const auto& err = reader.get_current_error(); !err.empty())
	{
		std::cerr << err;
		return 1;
	}
	//query for all segments and all atoms
	auto query_engine_ptr = reader.get_query_engine(params.segment_ids, params.atom_ids);
	auto& query_engine = *query_engine_ptr;
	if (const auto& err = query_engine.get_current_error(); !err.empty())
	{
		std::cerr << err;
		return 1;
	}
	uint32_t no_segments = reader.get_segments().size();
	uint32_t no_frames = reader.get_no_frames();

	std::cout << "Segments: " << no_segments<< ", Frames: " << no_frames << "\n";

	for (const auto& segment : reader.get_segments()) {
		std::cout << "Segment: " << segment.name << ", size: " << segment.size << ", type: " << segment_type_to_str(segment.type) << "\n";
	}

	std::ofstream output(params.output_file, std::ios::binary);
	if (!output)
	{
		std::cerr << "Error: Cannot open output file\n";
		return 1;
	}
	mdc::query_result result;
	for (uint32_t frame_id = 0; frame_id < no_frames; frame_id++)
	{
		if (!query_engine.query({ frame_id }, result))
		{
			std::cerr << "Error: " << query_engine.get_current_error();
			return 1;
		}
		assert(result.frames.size() == 1);
		auto coords = result.frames[0].coords;
		auto n_atoms = coords.size();
		output.write(reinterpret_cast<const char*>(coords.data()), coords.size()*3*sizeof(float));
		//for (const auto & coord : coords)
		//{
		//	//warning: this is not endian aware
		//	output.write(reinterpret_cast<const char*>(&coord.x), sizeof(coord.x));
		//	output.write(reinterpret_cast<const char*>(&coord.y), sizeof(coord.y));
		//	output.write(reinterpret_cast<const char*>(&coord.z), sizeof(coord.y));
		//}
		std::cerr << "Frame time: " << result.frames[0].time << "\tNo. atoms: " << n_atoms << "\n";
	
		
	}
	return 0;
}