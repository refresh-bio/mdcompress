#include <iostream>

#include "params.h"
#include "application.h"
#include "../common/version.h"

//using namespace std;

//std::string fn = "minimal.xtc";

CParams params;

bool parse_params(int argc, char** argv);
bool parse_params_compress(int argc, char** argv);
bool parse_params_decompress(int argc, char** argv);
bool parse_params_select(int argc, char** argv);
bool parse_make_desc(int argc, char** argv);
bool parse_convert(int argc, char** argv);
bool parse_params_info(int argc, char** argv);
bool parse_params_to_fp32(int argc, char** argv);
bool check_preset(int argc, char** argv);
void usage(int argc, char** argv);


// *******************************************************************************
bool check_preset(int argc, char** argv, CParams& params_to_update)
{
    const std::string preset_param = "--preset";
    auto p = std::find_if(argv, argv + argc, [&preset_param](const char* s) { return s == preset_param; });
    if (p == argv + argc)
        return true; //no preset, its ok

    ++p;
    if (p == argv + argc)
    {
        std::cerr << "Preset value must be provided\n";
        return false;
    }

    auto preset = preset_from_string(*p);
    if (!preset)
    {
        std::cerr << "Unknown preset: " << *p << std::endl;
        return false;
    }
    params_to_update.preset_values = preset_values_t::get(*preset);
    return true;
}

// *******************************************************************************
void usage(int argc, char** argv)
{
    CParams defaults;
    if (!check_preset(argc, argv, defaults)) //if user specified preset, we want to show how it sets defaults
        return;

    MDCOMPRESS_VER_PRINT(std::cerr);
    std::cerr << "Usage:\n";
    std::cerr << "mdcompress <mode> [options]\n";
    std::cerr << "Modes:\n";
    std::cerr << "  compress    - compress a trajectory (XTC/TRR/DCD/...) into MDC format\n";
    std::cerr << "  decompress  - decompress MDC file into a trajectory (XTC/TRR/DCD/...)\n";
    std::cerr << "  select      - decompress some frames from MDC file into a trajectory (XTC/TRR/DCD/...)\n";
    std::cerr << "  info        - info about contents of MDC file\n";
    std::cerr << "  make_desc   - create description file (for -d) from a topology (TPR/PSF/...)\n";

	std::cerr << "Options - compress mode\n";
    std::cerr << "  -i <file_name>               - input trajectory (XTC, TRR, DCD, ...)\n";
    std::cerr << "  -o <file_name>               - output file name (.mdc)\n";
    std::cerr << "  -d <file_name>               - description of the segments of a frame.\n"
                 "                                 May be omitted if --topology is given: mdcompress then\n"
                 "                                 builds the description from the topology automatically.\n";
    std::cerr << "  --topology <file_name>       - topology file (TPR, PSF, PDB, ...). It is stored inside the .mdc.\n"
                 "                                 If -d is not given, mdcompress infers the description from it.\n"
                 "                                 NOTE: coordinate-only trajectories (e.g. DCD) carry no topology,\n"
                 "                                 so either -d or --topology must be provided for them.\n"
                 "                                 NOTE 2: It may happen that mdcompress cannot build a description\n"
                 "                                 from the file specified with the --topology switch;\n"
                 "                                 in that case, -d must be provided.\n";
    std::cerr << "  --only-mol                   - compress only the molecule segments, skipping water/ions/other\n"
                 "                                 from the description/topology. Usually not needed: if the\n"
                 "                                 trajectory contains only molecules, mdcompress detects this and\n"
                 "                                 enables --only-mol automatically (with a warning).\n";
    std::cerr << "  (description check)          - whether -d or --topology is used, the description is validated\n"
                 "                                 against the trajectory's atom count before compressing:\n"
                 "                                   * exact match            -> used as-is;\n"
                 "                                   * only molecules match   -> --only-mol enabled automatically;\n"
                 "                                   * no match               -> compression is aborted with an\n"
                 "                                     explanation, and (when the description was inferred from a\n"
                 "                                     topology) a candidate is written to <output>.candidate.desc\n"
                 "                                     for you to edit and reuse with -d.\n";
    std::cerr << "  -l <int>                     - compression level (default: " << defaults.compression_level() << ")\n";
    std::cerr << "  --preset <preset>            - use one of preset parameters for -b, --subsegment-size,\n"
				 "                                 available presets:\n"
	             "                                   * default,\n"
			     "                                   * archive - use if the main decompression use scheme is decompress all,\n"
		         "                                   * trajectory - use if the main decompression use scheme is tracking trajectory\n"
		         "                                                  of some atoms,\n"
		         "                                   * frames - use if the main decompression use scheme is selecting subset of frames.\n"
		         "                                Current usage message prints defaults in regards to configured preset.\n";
    std::cerr << "  -t <int>                     - no. threads\n";
    //std::cerr << "  -a|--anchor-separation <int> - no. of frames between anchors (default: " << defaults.preset_values.anchor_separation() << ")\n"; 2025-04-02 -> changing this in parameters to batch size which is -a value + 1
    std::cerr << "  -b|--batch-size <int>        - no. of frames in a batch (default: " << defaults.preset_values.anchor_separation() + 1 << ")\n";
    std::cerr << "  -h|--max-history-size <int>  - no. of previous framed used to predict the current one (default: " << defaults.max_history_size() << "; max: " << defaults.max_history_size.max_val() << ")\n";
    std::cerr << "  --res <int>                  - min. resolution in fm (default: " << defaults.resolution() << "; 1000 fm = 0.01 Angstrom)\n";

    std::cerr << "  --subsegment <id1,id2,...>   - list of segments to split into subsegments (if no specified subsegment all)\n";
    std::cerr << "  --subsegment-size <int>      - number of atoms in a single subsegment (default: " << defaults.preset_values.subsegment_size() << ", 0 means don't use subsegments)\n";

    std::cerr << "  --n-frames-for-model         - no. of frames to build model (default: " << defaults.n_frames_for_model() << ")\n";
    std::cerr << "  --max-dist-in-model          - max. distance in segment of reference atoms (default: " << defaults.max_dist_in_model() << ")\n";

    std::cerr << "Options - decompress mode\n";
    std::cerr << "  -i <file_name>               - input file name (.mdc)\n";
    std::cerr << "  -o <file_name>               - output trajectory file name (format chosen by extension: .xtc/.trr/.dcd/...)\n";
    std::cerr << "  --topology <file_name>       - if a topology was stored at compression time, write it out to\n"
                 "                                 this path (the extension must match the original topology file)\n";

    std::cerr << "Options - select mode\n";
    std::cerr << "  -i <file_name>               - input file name (.mdc)\n";
    std::cerr << "  -o <file_name>               - output trajectory file name (format chosen by extension: .xtc/.trr/.dcd/...)\n";
    std::cerr << "  --topology <file_name>       - if a topology was stored at compression time, write it out to\n"
                 "                                 this path (the extension must match the original topology file)\n";
    std::cerr << "  --fid <int>                  - frame id (0-based)\n";
    std::cerr << "  --fr <int> <int|MAX>         - range of frame ids (0-based) (`MAX` or " << std::numeric_limits<int>::max() << " means last frame) \n";
    std::cerr << "  --stride <int>               - stride size (for range of frames) (default: " << defaults.stride() << ")\n";
    std::cerr << "  --segments <id1,id2,...>     - list of segment ids\n";
    std::cerr << "  --atoms <id1,id2,...>        - list of atoms to tracks, `id<n>` may be single id or closed interval in format start-end\n";

    std::cerr << "Options - info mode\n";
    std::cerr << "  -i <file_name>               - input file name (.mdc)\n";
    std::cerr << "  --full                       - print full info (include subsegments)\n";

    std::cerr << "Options - make_desc mode\n";
    std::cerr << "  -i <file_name>               - input topology file (TPR, PSF, ...)\n";
    std::cerr << "  -o <file_name>               - output desc file name\n";
    std::cerr << "  --only-mol                   - include only molecule (skip water and 'other')\n";

}

// *******************************************************************************
bool parse_params(int argc, char** argv)
{
    if (argc < 2)
    {
        usage(argc, argv);
        return false;
    }

    if (argv[1] == std::string("compress"))
        return parse_params_compress(argc, argv);
    else if (argv[1] == std::string("decompress"))
        return parse_params_decompress(argc, argv);
    else if (argv[1] == std::string("select"))
        return parse_params_select(argc, argv);
    else if (argv[1] == std::string("make_desc"))
        return parse_make_desc(argc, argv);
    else if (argv[1] == std::string("convert"))
        return parse_convert(argc, argv);
    else if (argv[1] == std::string("info"))
        return parse_params_info(argc, argv);
    else if(argv[1] == std::string("to_fp32"))
        return parse_params_to_fp32(argc, argv);

    usage(argc, argv);

    return false;
}

// *******************************************************************************
bool parse_params_compress(int argc, char** argv)
{
    params.mode = CParams::mode_t::compress;

    if (!check_preset(argc, argv, params))
		return false;

    for (int i = 2; i < argc; ++i)
    {
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.input_fn = argv[++i];
        else if (argv[i] == std::string("-o") && i + 1 < argc)
            params.output_fn = argv[++i];
        else if (argv[i] == std::string("-d") && i + 1 < argc)
            params.description_fn = argv[++i];
        else if (argv[i] == std::string("--topology") && i + 1 < argc)
            params.topology_fn = argv[++i];
        else if (argv[i] == std::string("--only-mol"))
            params.only_mol = true;
		else if (argv[i] == std::string("-t") && i + 1 < argc)
		{
			if (!params.no_threads.assign(std::stoi(argv[++i])))
			{
				std::cerr << "-t must be in range " << params.no_threads.range() << std::endl;
				return false;
			}
		}
		else if (argv[i] == std::string("-l") && i + 1 < argc)
		{
			if (!params.compression_level.assign(std::stoi(argv[++i])))
			{
				std::cerr << "-l must be in range " << params.compression_level.range() << std::endl;
				return false;
			}
		}
        else if(argv[i] == std::string("--native-lib-mode"))
			params.native_lib_mode = true;
		else if (argv[i] == std::string("--subsegment") && i + 1 < argc)
			params.subsegment_ids = params.split_string(argv[++i], ',');
		else if (argv[i] == std::string("--subsegment-size") && i + 1 < argc)
		{
            uint32_t subsegment_size{};

			if (!(std::istringstream(argv[++i]) >> subsegment_size))
			{
				std::cerr << "Invalid subsegment size: " << argv[i] << std::endl;
				return false;
			}

            if (!params.preset_values.subsegment_size.assign(subsegment_size))
            {
                std::cerr << "--subsegment-size must be in range " << params.preset_values.subsegment_size.range() << std::endl;
				return false;
            }
		}
        //else if ((argv[i] == std::string("-a") || argv[i] == std::string("--anchor-separation")) && i + 1 < argc)
        else if ((argv[i] == std::string("-b") || argv[i] == std::string("--batch-size")) && i + 1 < argc)
        {
            if (!params.preset_values.anchor_separation.assign(std::stoi(argv[++i]) - 1))
            {
				//std::cerr << "--anchor-separation must be in range " << params.preset_values.anchor_separation.range() << std::endl;
				std::cerr << "--batch-size must be in range " << params.preset_values.anchor_separation.range(1) << std::endl;
				return false;
            }
        }
        else if ((argv[i] == std::string("-h") || argv[i] == std::string("--max-history-size")) && i + 1 < argc)
        {
            if(!params.max_history_size.assign(std::stoi(argv[++i])))
            {
                std::cerr << "--max-history-size must be in range " << params.max_history_size.range() << std::endl;
                return false;
            }
        }
        else if (argv[i] == std::string("--res") && i + 1 < argc)
        {
            if (!params.resolution.assign(std::stoi(argv[++i])))
            {
                std::cerr << "--res must be in range " << params.resolution.range() << std::endl;
                return false;
            }
			params.resolution_set_by_user = true;
        }
        else if (argv[i] == std::string("--n-frames-for-model") && i + 1 < argc)
        {
            if(!params.n_frames_for_model.assign(std::stoi(argv[++i])))
			{
				std::cerr << "--n-frames-for-model must be in range " << params.n_frames_for_model.range() << std::endl;
				return false;
			}
        }
        else if (argv[i] == std::string("--max-dist-in-model") && i + 1 < argc)
        {
            if (!params.max_dist_in_model.assign(std::stoi(argv[++i])))
            {
				std::cerr << "--max-dist-in-model must be in range " << params.max_dist_in_model.range() << std::endl;
				return false;
            }
        }
        else if (argv[i] == std::string("--preset") && i + 1 < argc)
        {
            //we just skip this from parsing because it was handled in check_preset()
            //still we need to skipt this to not report as unknown option
            ++i;
        }
        else
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            return false;
        }
    }

    if (params.output_fn.empty())
    {
        std::cerr << "-o must be provided\n";
        return false;
    }

    if (params.output_fn.empty())
    {
        std::cerr << "-i must be provided\n";
        return false;
    }

    if (params.description_fn.empty() && params.topology_fn.empty())
    {
        std::cerr << "-d or --topology must be provided\n";
        return false;
    }

    return true;
}

// *******************************************************************************
bool parse_params_decompress(int argc, char** argv)
{
    params.mode = CParams::mode_t::decompress;

    for (int i = 2; i < argc; ++i)
    {
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.input_fn = argv[++i];
        else if (argv[i] == std::string("-o") && i + 1 < argc)
            params.output_fn = argv[++i];
        else if (argv[i] == std::string("--topology") && i + 1 < argc)
            params.topology_fn = argv[++i];
        else if (argv[i] == std::string("--native-lib-mode"))
            params.native_lib_mode = true;
        else if (argv[i] == std::string("-t") && i + 1 < argc)
        {
            if (!params.no_threads.assign(std::stoi(argv[++i])))
            {
                std::cerr << "-t must be in range " << params.no_threads.range() << std::endl;
                return false;
            }
        }
        else
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            return false;
        }
    }

    if (params.input_fn.empty() || params.output_fn.empty())
    {
        std::cerr << "Input and output file names must be provided\n";
        return false;
    }

    return true;
}

// *******************************************************************************
bool parse_params_info(int argc, char** argv)
{
    params.mode = CParams::mode_t::info;

    for (int i = 2; i < argc; ++i)
    {
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.input_fn = argv[++i];
        else if (argv[i] == std::string("--full"))
            params.full_info = true;
        else
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            return false;
        }
    }

    if (params.input_fn.empty())
    {
        std::cerr << "Input file name must be provided\n";
        return false;
    }

    return true;
}

// *******************************************************************************
bool parse_params_select(int argc, char** argv)
{
    params.mode = CParams::mode_t::select;

    for (int i = 2; i < argc; ++i)
    {
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.input_fn = argv[++i];
        else if (argv[i] == std::string("-o") && i + 1 < argc)
            params.output_fn = argv[++i];
        else if (argv[i] == std::string("--topology") && i + 1 < argc)
            params.topology_fn = argv[++i];
        else if (argv[i] == std::string("--native-lib-mode"))
            params.native_lib_mode = true;
        else if (argv[i] == std::string("-t") && i + 1 < argc)
        {
            if (!params.no_threads.assign(std::stoi(argv[++i])))
            {
                std::cerr << "-t must be in range " << params.no_threads.range() << std::endl;
                return false;
            }
        }
        else if (argv[i] == std::string("--fid") && i + 1 < argc)
        {
            params.frame_id = std::stoi(argv[++i]);
            if (params.frame_id < 0)
            {
                std::cerr << "--fid must be >= 0\n";
                return false;
            }
        }
        else if (argv[i] == std::string("--fr") && i + 2 < argc)
        {
            params.range_id_from = std::stoi(argv[++i]);

            if (params.range_id_from < 0)
            {
                std::cerr << "first value of --fr must be >=0\n";
                return false;
            }

            std::string_view tmp = argv[++i];
            if (tmp == "MAX" or tmp == "max" or tmp == "Max")
                params.range_id_to = std::numeric_limits<int>::max();
            else
                params.range_id_to = std::stoi(tmp.data());

            if (params.range_id_to < 0)
            {
                std::cerr << "second value of --fr must be >= 0 or `MAX`\n";
                return false;
            }

            if (params.range_id_to < params.range_id_from)
            {
                std::cerr << "second value of --fr must be >= first value\n";
                return false;
            }
        }
        else if (argv[i] == std::string("--stride") && i + 1 < argc)
        {
            if (!params.stride.assign(std::stoi(argv[++i])))
            {
				std::cerr << "--stride must be in range " << params.stride.range() << std::endl;
				return false;
            }
        }
        else if (argv[i] == std::string("--segments") && i + 1 < argc)
        {
			params.segment_ids = params.split_string(argv[++i], ',');
        }
        else if (argv[i] == std::string("--atoms") && i + 1 < argc)
        {
            auto splitted = params.split_string(argv[++i], ',');
            for (auto& s : splitted)
			{
				auto dash = s.find('-');
				if (dash != std::string::npos)
				{
					int start = std::stoi(s.substr(0, dash).c_str());
					int end = std::stoi(s.substr(dash + 1).c_str());
					for (int j = start; j <= end; ++j)
						params.atom_ids.push_back(j);
				}
				else
				{
					params.atom_ids.push_back(std::stoi(s.c_str()));
				}
			}
        }
        else
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            return false;
        }
    }

    if (params.input_fn.empty() || params.output_fn.empty())
    {
        std::cerr << "Input and output file names must be provided\n";
        return false;
    }

    return true;
}

// *******************************************************************************
bool parse_make_desc(int argc, char** argv)
{
    params.mode = CParams::mode_t::make_desc;
	for (int i = 1 ; i < argc ; ++i)
	{
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.params_make_desc.input_path = argv[++i];
        else if (argv[i] == std::string("-o") && i + 1 < argc)
            params.params_make_desc.output_path = argv[++i];
        else if (argv[i] == std::string("--only-mol"))
            params.only_mol = true;
	}
    if (params.params_make_desc.input_path.empty() || params.params_make_desc.output_path.empty())
    {
        std::cerr << "Input and output file names must be provided\n";
        return false;
    }
    return true;
}

// *******************************************************************************
bool parse_convert(int argc, char** argv)
{
    params.mode = CParams::mode_t::convert;
	for (int i = 1 ; i < argc ; ++i)
	{
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.input_fn = argv[++i];
        else if (argv[i] == std::string("-o") && i + 1 < argc)
            params.output_fn = argv[++i];
	}
    if (params.input_fn.empty() || params.output_fn.empty())
    {
        std::cerr << "Input and output file names must be provided\n";
        return false;
    }
    return true;
}

// *******************************************************************************
bool parse_params_to_fp32(int argc, char** argv)
{
    params.mode = CParams::mode_t::to_fp32;

    for (int i = 2; i < argc; ++i)
    {
        if (argv[i] == std::string("-i") && i + 1 < argc)
            params.input_fn = argv[++i];
        else if (argv[i] == std::string("-o") && i + 1 < argc)
            params.output_fn = argv[++i];
        else if (argv[i] == std::string("--native-lib-mode"))
            params.native_lib_mode = true;
        else if (argv[i] == std::string("-d") && i + 1 < argc)
            params.description_fn = argv[++i];
        else
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            return false;
        }
    }

    if (params.input_fn.empty() || params.output_fn.empty())
    {
        std::cerr << "Input and output file names must be provided\n";
        return false;
    }

    return true;
}

// *******************************************************************************
int main(int argc, char **argv)
{
    if (!parse_params(argc, argv))
    {
        return 1;
    }

    CApplication app;

    if (!app.run(params))
        return 1;

    return 0;
}
