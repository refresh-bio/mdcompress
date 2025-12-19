#include <cstdint>
#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>
#include <thread>

#ifdef _MSC_VER
#include <mimalloc-override.h>
#endif

#include "application.h"
#include "../libs/refresh/parallel_queues/lib/thread-control.h"

#include <ranges>
#include <refresh/compression/lib/lzss_compress.h>
#include <refresh/compression/lib/lzss_decompress.h>

#include "../common/mdc_file_version.h"

#include "input_reader.h"
#include "output_writer.h"
#include "desc_builder.h"

// *******************************************************************************
void CApplication::reduce_segment_desc_to_molecules_only()
{
    std::vector<segment_desc_t> reduced;
    for (const auto& desc : segment_desc)
        if (desc.type == mdc::segment_type::molecule)
            reduced.push_back(desc);
    segment_desc = std::move(reduced);
}

// *******************************************************************************
bool CApplication::load_desc_from_trajectory_file()
{
    try
    {
        segment_desc = build_desc(params.topology_fn);
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what();
        return false;
    }
    catch (...)
    {
        std::cerr << "Unknown error cautch at " << __FILE__ << ":" << __LINE__ << "\n";
        return false;
    }
    return true;
}


class TopologyInfo
{
public:
    enum class CompressionType {
        plain, //topology is stored without compression
        refresh_LZSS //topology is stored with refresh_lzss
    };
private:
    bool store_topology{};
    CompressionType compression_type{};
    std::string filename;
public:
    TopologyInfo(const std::string& topology_fn, CompressionType compr_type = CompressionType::refresh_LZSS)
    {
        if (topology_fn.empty())
        {
            store_topology = false;
            return;
        }
        store_topology = true;
        filename = std::filesystem::path(topology_fn).filename().string();
        compression_type = compr_type;
    }
    void Serialize(std::vector<uint8_t>& res) const
    {
        Serializer serializer;
        if (not store_topology)
        {
            serializer.append_str("no_topology");
            res = serializer.get_data();
            return;
        }

        serializer.append_str("with_topology");

        switch (compression_type)
        {
        case CompressionType::plain:
            serializer.append_str("plain");
            break;
        case CompressionType::refresh_LZSS:
            serializer.append_str("refresh_LZSS");
            break;
        }

        serializer.append_str(filename);
        res = serializer.get_data();
    }
    TopologyInfo(const std::vector<uint8_t>& serialized)
    {
        Serializer serializer;
        serializer.set_data(serialized);
        std::string tmp = serializer.load_str();

        if (tmp == "no_topology")
        {
            store_topology = false;
            return;
        }

        if (tmp != "with_topology")
        {
            std::cerr << "Error: expected `with_topology` or `no_topology` but get " << tmp << "\n";
            exit(1);
        }
        store_topology = true;

        std::string _compression_type = serializer.load_str();

        if (_compression_type == "plain")
        {
            compression_type = CompressionType::plain;
        }
        else if (_compression_type == "refresh_LZSS")
        {
            compression_type = CompressionType::refresh_LZSS;
        }
        else
        {
            std::cerr << "Error: unknown compression type: " << _compression_type << "\n";
            exit(1);
        }
        filename = serializer.load_str();
    }

    bool IsTopologyStored() const
    {
        return store_topology;
    }
    CompressionType GetCompressionType() const
    {
        return compression_type;
    }
    const std::string& GetFilename() const
    {
        return filename;
    }
};

void store_topology(const std::string& topology_fn, refresh::archive_output& ar_out)
{
    //ar_out must be opened
    TopologyInfo topology_info(topology_fn, TopologyInfo::CompressionType::refresh_LZSS);

    std::vector<uint8_t> serialized_topology_info;
    topology_info.Serialize(serialized_topology_info);
    ar_out.add_part(ar_out.register_stream("topology-info"), serialized_topology_info);

    if (not topology_info.IsTopologyStored())
        return;

    auto topology_file_size = std::filesystem::file_size(topology_fn);
    std::vector<uint8_t> file_data(topology_file_size);

    FILE* f = fopen(topology_fn.c_str(), "rb");

    if (!f)
    {
        std::cerr << "Error: cannot open file " << topology_fn << "\n";
        exit(1);
    }
    if (topology_file_size != fread(file_data.data(), 1, topology_file_size, f))
    {
        std::cerr << "Error: reading topology file: " << topology_fn << "\n";
        fclose(f);
        exit(1);
    }
    fclose(f);
    std::vector<uint8_t> to_store;
    if (topology_info.GetCompressionType() == TopologyInfo::CompressionType::plain)
    {
        to_store = std::move(file_data);
    }
    else if (topology_info.GetCompressionType() == TopologyInfo::CompressionType::refresh_LZSS)
    {
        if (not refresh::lzss_rc_compress(file_data, to_store, 6))
        {
            std::cerr << "Error: cannot compress file " << topology_fn << "\n";
            exit(1);
        }
    }
    else
    {
        std::cerr << "Error: unsupported compression type, code need to be updated\n";
        exit(1);
    }

    ar_out.add_part(ar_out.register_stream("topology"), to_store);
}

std::optional<TopologyInfo> read_topology_info(refresh::archive_input& ar_in)
{
    int topology_info_id = ar_in.get_stream_id("topology-info");
    if (topology_info_id < 0)
    {
        std::cerr << "Corrupted input file, cannot find `topology-info` stream" << std::endl;
        return std::nullopt;
    }
    std::vector<uint8_t> topology_info_serialized;
    uint64_t metadata;
    if (not ar_in.get_part(topology_info_id, topology_info_serialized, metadata))
    {
        std::cerr << "Error: cannot read `topology-info` stream\n";
        return std::nullopt;
    }
    TopologyInfo topology_info(topology_info_serialized);

    return topology_info;
}

bool extract_topology(refresh::archive_input& ar_in, const std::string& topology_fn)
{
    auto topology_info = read_topology_info(ar_in);
    if (not topology_info)
        return false;

    if (not topology_info->IsTopologyStored())
    {
        std::cerr << "Error: topology was not stored during compression\n";
        return false;
    }

    std::vector<uint8_t> stored_data;
    int topology_id = ar_in.get_stream_id("topology");
    if (topology_id < 0)
    {
        std::cerr << "Corrupted input file, cannot find `topology` stream" << std::endl;
        return false;
    }

    uint64_t metadata;

    if (not ar_in.get_part(topology_id, stored_data, metadata))
    {
        std::cerr << "Error: cannot read `topology` stream\n";
        return false;
    }

    std::vector<uint8_t> to_store;
    if (topology_info->GetCompressionType() == TopologyInfo::CompressionType::plain)
    {
        to_store = std::move(stored_data);
    }
    else if (topology_info->GetCompressionType() == TopologyInfo::CompressionType::refresh_LZSS)
    {
        if (not refresh::lzss_rc_decompress(stored_data, to_store))
        {
            std::cerr << "Error: cannot decompress topology file\n";
            return false;
        }
    }
    else
    {
        std::cerr << "Error: unsupported compression type, code need to be updated\n";
        return false;
    }

    auto stored_extension = std::filesystem::path(topology_info->GetFilename()).extension().string();
    auto wanted_extension = std::filesystem::path(topology_fn).extension().string();
    if (stored_extension != wanted_extension)
    {
        std::cerr << "Error: topology file was stored with extension: " << stored_extension << ", the same extension must be specified for output file name, but " << wanted_extension << " was used\n";
        return false;
    }
    FILE* f = fopen(topology_fn.c_str(), "wb");
    if (!f)
    {
        std::cerr << "Error: cannot open output file " << topology_fn << "\n";
        return false;
    }
    if (to_store.size() != fwrite(to_store.data(), 1, to_store.size(), f))
    {
        std::cerr << "Error: cannot write topology data to file " << topology_fn << "\n";
        fclose(f);
        return false;
    }
    fclose(f);

    return true;
}

// *******************************************************************************
bool CApplication::load_segment_description_file()
{
    std::ifstream ifs(params.description_fn);

    if (!ifs)
    {
        std::cerr << "Cannot open file: " << params.description_fn << std::endl;
        return false;
    }

    segment_desc.clear();

    std::string s_type;
    segment_type_t type;
    std::string name;
    uint32_t n_atoms;

    while (ifs >> s_type >> name >> n_atoms)
    {
        if (s_type == "MOLECULE" || s_type == "MOL")
            type = segment_type_t::molecule;
        else if (s_type == "WATER" || s_type == "WAT")
            type = segment_type_t::water;
        else if (s_type == "OTHER" || s_type == "OTH")
            type = segment_type_t::other;
        else if (s_type == "NONE" || s_type == "NON")
            type = segment_type_t::none;
        else
        {
            std::cerr << "Unknown segment type: " << s_type << std::endl;
            return false;
        }

        segment_desc.emplace_back(type, name, n_atoms);
    }

    return true;
}

// *******************************************************************************
//the general idea is that we split segments as defined in desc file into subsegments
//it is represented such that the first subsegment of a segment contains original segment name
//the remaining subsegments are named as empty string
//if user does not specify any subsegments then we keep the original segments but we still name them subsegments
//because this representation is general enough
void CApplication::expand_segments()
{
    //if segments to be splitted into subsegments were not specified take all except water
    //actually for parallelization it may be better to also subsegment water
    if (params.subsegment_ids.empty())
        for (const auto& segment : segment_desc)
            //if (segment.type != mdc::segment_type::water)
                params.subsegment_ids.push_back(segment.name);

    std::vector<segment_desc_t> expanded_segment_desc;
    for (const auto& segment : segment_desc)
    {
        //if not specified to be splitted into subsegments just keep a single subsegment (whole segment)
        if (std::find_if(params.subsegment_ids.begin(), params.subsegment_ids.end(), [&](const std::string& name) {return segment.name == name; }) == params.subsegment_ids.end())
        {
            expanded_segment_desc.emplace_back(segment);
            continue;
        }
        auto left = segment.size;
        bool first = true;
        while (left)
        {
            auto size = std::min(left, params.preset_values.subsegment_size());
            expanded_segment_desc.emplace_back(segment.type, first ? segment.name : "", size);
            first = false;
            left -= size;
        }
    }
    segment_desc = std::move(expanded_segment_desc);
}

// *******************************************************************************
const char* CApplication::segment_type_t_to_str(const segment_type_t type) const
{
    switch (type) {
    case segment_type_t::molecule:
        return "MOL";
    case segment_type_t::water:
        return "WAT";
    case segment_type_t::other:
        return "OTHER";
    case segment_type_t::none:
        return "NON";
    default:
        assert(0);              // Never should be here
        return "";
    }
}

// *******************************************************************************
bool CApplication::compress()
{
    auto t1 = std::chrono::high_resolution_clock::now();

    std::unique_ptr<InputReader> input_reader;
    auto traj_format = get_traj_format(params.input_fn, false);

#if defined(USE_CHEMFILES) && defined(USE_XDRFILE)
    if (params.native_lib_mode)
    {
        switch (traj_format)
        {
        case traj_file_format_t::trr:
            input_reader = std::make_unique<CTRRReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
        case traj_file_format_t::xtc:
            input_reader = std::make_unique<CXTCReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
        default:
            return false;
        }
    }
    else
        input_reader = std::make_unique<CTrajReader>(params.resolution_set_by_user ? params.resolution() : 0);
#elif defined(USE_CHEMFILES)
    input_reader = std::make_unique<CTrajReader>(params.resolution_set_by_user ? params.resolution() : 0);
#else
    switch (traj_format)
    {
    case traj_file_format_t::trr:
        input_reader = std::make_unique<CTRRReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
    case traj_file_format_t::xtc:
        input_reader = std::make_unique<CXTCReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
    default:
        return false;
    }
#endif

    CMDCCompress mdc_compress(params.preset_values.anchor_separation(), params.max_history_size(), params.compression_features[params.compression_level()]);
    refresh::archive_output ar_out;

    std::vector<uint64_t> stream_sizes;

    //if description file was provided use it, it has a predeecence over --topology file
    if (not params.description_fn.empty())
    {
        if (!load_segment_description_file())
            return false;
    }
    else
    {
        if (!load_desc_from_trajectory_file())
            return false;
    }
    if (params.only_mol) //reduce segment description to only point to molecules
        reduce_segment_desc_to_molecules_only();

    if (params.preset_values.subsegment_size() != 0)
        expand_segments();

    if (!input_reader->open(params.input_fn, segment_desc))
    {
        std::cerr << "Cannot open: " << params.input_fn << std::endl;
        return false;
    }

    // *** Model construction
    frame_t frame;
    frame_t first_frame;

    if (params.n_frames_for_model() <= 1)
    {
        if (!input_reader->get_frame(frame))
        {
            std::cerr << "Cannot read first frame" << std::endl;
            return 0;
        }

        first_frame = frame;

        // *** Model building
        for (uint32_t i = 0; i < (uint32_t)frame.segments.size(); ++i)
            mdc_compress.BuildModel(i, frame.segments[i], params);
    }
    else
    {
        std::vector<std::vector<segment_t>> segments;

        for (uint32_t i = 0; i < params.n_frames_for_model(); ++i)
        {
            if (!input_reader->get_frame(frame))
            {
                if (i == 0)
                {
                    std::cerr << "Cannot read first frame" << std::endl;
                    return 0;
                }
                else
                    break;
            }

            if (i == 0)
            {
                segments.resize(frame.segments.size());
                first_frame = frame;
            }

            for (uint32_t j = 0; j < (uint32_t)segments.size(); ++j)
                if (i == 0 || segment_desc[j].type == segment_type_t::molecule)
                    segments[j].push_back(frame.segments[j]);
        }

        // *** Model building
        for (uint32_t i = 0; i < (uint32_t)segments.size(); ++i)
            mdc_compress.BuildModel(i, segments[i], params);

        // *** Adding orientations
        for (uint32_t i = 0; i < (uint32_t)segments.size(); ++i)
            mdc_compress.AddOrientation(i, first_frame.segments[i], params);

        input_reader->close();
        if (!input_reader->open(params.input_fn, segment_desc))
        {
            std::cerr << "Cannot open: " << params.input_fn << std::endl;
            return false;
        }

        if (!input_reader->get_frame(frame))
        {
            std::cerr << "Cannot read first frame" << std::endl;
            return 0;
        }
    }

	stream_sizes.resize(frame.segments.size() + 1);

    // *** Model storing
    if (!ar_out.open_file_buffered(params.output_fn, false))
    {
        std::cerr << "Cannot create file: " << params.output_fn << std::endl;
        return false;
    }

    //default constructed mdc_file_version contains the version of the software that created the file and the version of the file format
    if (!mdc_file_version{}.store(ar_out))
    {
        std::cerr << "Cannot store file version" << std::endl;
		return false;
	}

    store_topology(params.topology_fn, ar_out);

    CMDCCompress::packed_frames_t packed_frames;
    packed_t packed;

    std::vector<int> stream_ids;
    int frame_id = ar_out.register_stream("frame-id");
//    size_t no_segments = frame.segments.size();

    for (uint32_t i = 0; i < (uint32_t) frame.segments.size(); ++i)
    {
        mdc_compress.SerializeModel(i, packed);
        ar_out.add_part(ar_out.register_stream("model-" + std::to_string(i)), packed);

        stream_sizes[i] = packed.size();

        stream_ids.emplace_back(ar_out.register_stream("stream-" + std::to_string(i)));
    }

    // *** Compression
    mdc_compress.StartCompression();

    refresh::thread_control tc(params.no_threads());

	refresh::parallel_queue<std::vector<frame_t>> q_batches(4, 1, "q_batches");

    std::thread t_readed([&] {
        tc.acquire();

        std::vector<frame_t> batch;


        std::cerr << "Frame time: " << frame.desc.time << "\tNo. atoms: " << frame.desc.n_atoms << std::endl;
        batch.emplace_back(std::move(frame));

        while (input_reader->get_frame(frame)) {

            std::cerr << "Frame time: " << frame.desc.time << "\tNo. atoms: " << frame.desc.n_atoms << std::endl;

            batch.emplace_back(std::move(frame));
            if (batch.size() > params.preset_values.anchor_separation())
            {
                tc.release();
                //q_frames.push(std::move(frame));
                q_batches.push(std::move(batch));
                tc.acquire();

                batch.clear();
            }
        }
        tc.release();
        if (!batch.empty())
            q_batches.push(std::move(batch));

		q_batches.mark_completed();
    });

    refresh::parallel_queue<compress_task_desc_t> q_workers(tc.max_running() * 2, 1, "q_workers");

    std::thread t_compressor([&] {
        std::vector<compress_task_desc_t> generated_tasks;

        std::vector<frame_t> batch;
        while (q_batches.pop(batch))
        {
            tc.acquire();
            mdc_compress.SplitBatch(batch, generated_tasks);
            tc.release();
            for (auto& task : generated_tasks)
                q_workers.push(std::move(task));
        }
        q_workers.mark_completed();

        frame_t frame;
        });

    struct storer_task_desc
    {
	    enum class type_t {segment, desc} type;
        size_t segment_id = std::numeric_limits<size_t>::max();
        packed_t data;
    };

    refresh::parallel_priority_queue<storer_task_desc> q_storer(tc.max_running() * 2, tc.max_running(), "q_storer");

    std::vector<std::thread> workers;
    workers.reserve(tc.max_running());
    for (uint32_t worker_id = 0; worker_id < tc.max_running(); ++worker_id)
        workers.emplace_back([&]
    {
        CMDCCompressWorker worker;
        mdc_compress.InitWorker(worker);
        compress_task_desc_t task_desc;
        packed_t packed;
        while (q_workers.pop(task_desc))
        {
            tc.acquire();
			storer_task_desc out_task;
	        switch (task_desc.work_type)
	        {
            case compress_task_desc_t::work_type_t::desc:
                worker.HandleDesc(*task_desc.batch, packed);
                out_task.type = storer_task_desc::type_t::desc;
                out_task.data = std::move(packed);
                break;
            case compress_task_desc_t::work_type_t::segment:
                worker.HandleSegment(task_desc.segment_id, *task_desc.batch, packed);
                out_task.type = storer_task_desc::type_t::segment;
                out_task.segment_id = task_desc.segment_id;
                out_task.data = std::move(packed);
                break;
            default:
                assert(false);
	        }
            tc.release();

            q_storer.push(task_desc.priority, std::move(out_task));
        }
        q_storer.mark_completed();
    });

    std::thread t_storer([&]
        {
        storer_task_desc task;
        while (q_storer.pop(task))
        {
            tc.acquire();
            switch (task.type)
            {
            case storer_task_desc::type_t::desc:
                ar_out.add_part(frame_id, task.data);
                stream_sizes.back() += task.data.size();
                break;
            case storer_task_desc::type_t::segment:
                ar_out.add_part(stream_ids[task.segment_id], task.data);
                stream_sizes[task.segment_id] += task.data.size();
                break;
            default:
                assert(false);
            }
            tc.release();
        }

     });

	t_readed.join();
	t_compressor.join();
	for (auto& th : workers)
		th.join();
	
    t_storer.join();

    packed_t packed_anchor_desc;
    mdc_compress.Finalize(packed_anchor_desc);

    packed_t segment_desc_packed;
	mdc_compress.SerializeSegmentDesc(segment_desc_packed);
    
    ar_out.add_part(ar_out.register_stream("segment-desc"), segment_desc_packed);
    ar_out.add_part(ar_out.register_stream("anchor-ids"), packed_anchor_desc);

    ar_out.close();

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = t2 - t1;

    std::cerr << "Compressing time: " << diff.count() << "s\n";

    //combines summaries by subsegments into summaries by segments
    auto get_summary = [this, &stream_sizes, &mdc_compress]()
    {
		size_t i{};
		auto next = [&i, this, &stream_sizes, &mdc_compress](size_t& stream_size, size_t& no_atoms)
		{
			if (i >= segment_desc.size())
				return false;

			assert(segment_desc[i].name != "");

			stream_size = stream_sizes[i];
			no_atoms = mdc_compress.GetNoAtomsInSubSegment(i);
			++i;
			for ( ;i < segment_desc.size(); ++i)
			{
				if (segment_desc[i].name != "")
					break;

				stream_size += stream_sizes[i];
				no_atoms += mdc_compress.GetNoAtomsInSubSegment(i);
			}
			return true;
		};

        std::vector<std::pair<size_t, size_t>> summary;
        size_t stream_size{};
        size_t no_atoms{};
        while (next(stream_size, no_atoms))
			summary.emplace_back(stream_size, no_atoms);

		return summary;
	};
	auto summary = get_summary();

	std::cerr << "Segment sizes:\n";
	for (size_t i = 0; i < summary.size(); ++i)
		std::cerr << "    Segment " << i << ": " << summary[i].first << " bytes\n";
	std::cerr << "    Frame desc: " << stream_sizes.back() << " bytes\n";

	std::cerr << "Segment ratios:\n";
	for (size_t i = 0; i < summary.size(); ++i)
		std::cerr << "    Segment " << i << ": " << 8.0 * summary[i].first / (mdc_compress.GetNoFrames() * summary[i].second) << " bits/atom\n";

	return true;
}

// *******************************************************************************
bool CApplication::decompress() const
{
	// At the moment decompress is just a copy of select with no limiting switches
	return select();
}

std::vector<uint32_t> CApplication::make_frame_ids(uint32_t tot_frames) const
{
    std::vector<uint32_t> frame_ids;
    if (params.frame_id >= 0)
    {
        if (static_cast<uint32_t>(params.frame_id) >= tot_frames)
        {
            std::cerr << "Error: frame " << params.frame_id << " was requested but there are only " << tot_frames << " frames\n";
            return {};
        }
        frame_ids.push_back(params.frame_id);
    }
    else if (params.range_id_from >= 0)
    {
        if (static_cast<uint32_t>(params.range_id_from) >= tot_frames)
        {
            std::cerr << "Error: the first frame requested is " << params.range_id_from << " but there are only " << tot_frames << " frames\n";
            return {};
        }
        auto range_id_to = params.range_id_to; //copy because we dont want to check content of Params

        if (range_id_to == std::numeric_limits<int>::max())
            range_id_to = static_cast<int>(tot_frames) - 1;

        if (static_cast<uint32_t>(range_id_to) >= tot_frames)
        {
            std::cerr << "Error: the last frame requested is " << range_id_to << " but there are only " << tot_frames << " frames\n";
            return {};
        }

        assert(range_id_to >= params.range_id_from);

        uint32_t stride = params.stride() > 0 ? params.stride() : 1;

        for (uint32_t i = static_cast<uint32_t>(params.range_id_from); i <= static_cast<uint32_t>(range_id_to); i += stride)
            frame_ids.emplace_back(i);
    }
    else
    {
        frame_ids.resize(tot_frames);
        std::iota(frame_ids.begin(), frame_ids.end(), 0ul);
    }
    return frame_ids;
}

template<typename FUNC_T>
class AtExit
{
    FUNC_T func;
public:
    template<typename T>
    AtExit(T&& func) : func(std::forward<FUNC_T>(func)) {}
    AtExit(const AtExit<FUNC_T>& rhs) = delete;
    AtExit(AtExit<FUNC_T>&& rhs) = delete;
    AtExit<FUNC_T>& operator=(const AtExit<FUNC_T>& rhs) = delete;
    AtExit<FUNC_T>& operator=(AtExit<FUNC_T>&& rhs) = delete;
    ~AtExit()
    {
        func();
    }
};
template<typename F>
AtExit(F&&) -> AtExit<std::decay_t<F>>;


// *******************************************************************************
bool CApplication::select() const
{
    //mkokot_TODO: at this point I'm doing it this way just to test
    //but it must be moved to the API and this code here fixed to use API
    if (not params.topology_fn.empty())
    {
        refresh::archive_input ar_in;
        if (not ar_in.open_file_unbuffered(params.input_fn, false))
        {
            std::cerr << "Error: cannot open " << params.input_fn << "\n";
            return false;
        }
        if (not extract_topology(ar_in, params.topology_fn))
            return false;
    }

    auto start = std::chrono::high_resolution_clock::now();
    mdc::reader reader(params.input_fn);

    if (const auto& err = reader.get_current_error(); !err.empty())
    {
        std::cerr << err;
        return false;
    }

    const auto& anchor_ids = reader.get_anchor_ids();
    if (const auto& err = reader.get_current_error(); !err.empty())
    {
        std::cerr << err;
        return false;
    }

    std::vector<uint32_t> frame_ids = make_frame_ids(reader.get_no_frames());

    //means some error which was reported inside make_frame_ids
	if (frame_ids.empty())
        return false;

    struct task_desc
    {
        std::span<const uint32_t> frame_ids;
        uint32_t priority;
    };

    refresh::thread_control tc(params.no_threads());

    auto make_decompress_tasks_descs = [this] (
        const std::vector<uint32_t>& frame_ids,
        const std::vector<uint32_t>& anchor_ids)
        {
            // ! ! ! frames id must be sorted ! ! !

            std::vector<task_desc> tasks_descs;
            bool need_to_locate = true;
            uint32_t last_frame_id_in_batch{};
            int32_t start{};
            uint32_t priority{};
            for (int32_t i = 0; i < std::ssize(frame_ids); ++i)
            {
                if (need_to_locate)
                {
                    need_to_locate = false;
                    const auto cur_anchor_id = locate_frame(anchor_ids, frame_ids[i]).first;

                    if (cur_anchor_id + 1 == anchor_ids.size())
                        last_frame_id_in_batch = frame_ids.back();
                    else
                        last_frame_id_in_batch = anchor_ids[cur_anchor_id + 1] - 1;
                }

                if (frame_ids[i] > last_frame_id_in_batch)
                {
                    tasks_descs.emplace_back(task_desc{ std::span(frame_ids.begin() + start, i - start), priority++ });
                    start = i;
                    need_to_locate = true;
                    --i;
                }
            }
            tasks_descs.emplace_back(task_desc{ std::span(frame_ids.begin() + start, frame_ids.size() - start), priority++ });

            return tasks_descs;
        };

    auto task_descs = make_decompress_tasks_descs(frame_ids, anchor_ids);
    std::unique_ptr<OutputWriter> output_writer;

    auto traj_format = get_traj_format(params.output_fn, true);

#if defined(USE_CHEMFILES) && defined(USE_XDRFILE)
    if (params.native_lib_mode)
    {
        switch (traj_format)
        {
        case traj_file_format_t::trr:
            output_writer = std::make_unique<CTRRWriter>(); break;
        case traj_file_format_t::xtc:
            output_writer = std::make_unique<CXTCWriter>(); break;
        default:
            return false;
        }
    }
    else
        output_writer = std::make_unique<CTrajWriter>();
#elif defined(USE_CHEMFILES)
    output_writer = std::make_unique<CTrajWriter>();
#else
    switch (traj_format)
    {
    case traj_file_format_t::trr:
        output_writer = std::make_unique<CTRRWriter>(); break;
    case traj_file_format_t::xtc:
        output_writer = std::make_unique<CXTCWriter>(); break;
    default:
        return false;
    }
#endif

    if (!output_writer->open(params.output_fn))
    {
        std::cerr << "Cannot create: " << params.output_fn << std::endl;
        return false;
    }

    using q_frames_t = refresh::parallel_priority_queue<mdc::query_result>;
    q_frames_t q_frames(tc.max_running() * 3, tc.max_running(), "q_frames");

    std::thread t_storer([&] {
        mdc::query_result result;
        while (q_frames.pop(result))
        {
            tc.acquire();

            output_writer->add_frames(result);
            for (uint32_t frame_id = 0; frame_id < result.frames.size(); ++frame_id)
                std::cerr << "Frame time: " << result.frames[frame_id].time << "\tNo. atoms: " << result.frames[frame_id].coords.size() << std::endl;

            tc.release();
        }
        });

    std::vector<std::thread> decompr_threads;
    decompr_threads.reserve(tc.max_running());
    std::atomic<size_t> task_id{};
    std::atomic<bool> was_error = false;
    for (std::ptrdiff_t tid = 0 ; tid < tc.max_running() ; ++tid)
    {
        decompr_threads.emplace_back([&]
        {
            AtExit raii_for_mark_completed([&q_frames]
                {
                    q_frames.mark_completed();
                });

            auto query_engine_ptr = reader.get_query_engine(params.segment_ids, params.atom_ids);
            auto& query_engine = *query_engine_ptr;
            if (const auto& err = query_engine.get_current_error(); !err.empty())
            {
                std::cerr << err;

                was_error = true;
                task_id = task_descs.size(); // effectively stop all decmpr threads
                return;
            }

            mdc::query_result result;

            while (true)
            {
                const size_t my_task_id = task_id++;
				if (my_task_id >= task_descs.size())
                    break;

                tc.acquire();
                const bool query_ok = query_engine.query(task_descs[my_task_id].frame_ids, result);
                tc.release();

                if (!query_ok)
                {
                    std::cerr << std::string("Error: ") + query_engine.get_current_error() + "\n";

                    was_error = true;
                    task_id = task_descs.size(); // effectively stop all decmpr threads
                    return;
                }

                q_frames.push(task_descs[my_task_id].priority, std::move(result));
            }
        });
    }

    for (auto& dt : decompr_threads)
        dt.join();

    t_storer.join();

    if (was_error)
        return false;

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Decompress time: " << std::chrono::duration<double>(end - start).count() << "s\n";

    return true;
}
// *******************************************************************************
bool CApplication::convert() const
{
      try {
        // Open input, let chemfiles guess format from extension (.xtc, .trr, .pdb, …)
        chemfiles::Trajectory input(params.input_fn);

        // Same for output. The format is chosen by the extension of output_fn.
        chemfiles::Trajectory output(params.output_fn, 'w');

        chemfiles::Frame frame;
        while (!input.done()) {
            frame = input.read();   // read next frame
            output.write(frame); // write to output in the new format
        }
    } catch (const chemfiles::Error& e) {
        // TODO: plug into your logging / error handling
        std::cerr << "Conversion failed: " << e.what() << '\n';
        return false;
    }

    return true;
}
// *******************************************************************************
bool CApplication::make_desc() const
{
    std::ofstream out(params.params_make_desc.output_path, std::ios::binary);
    if (!out)
    {
        std::cerr << "Error: cannot open file " << params.params_make_desc.output_path << "\n";
        return false;
    }
    try
    {
        auto desc = build_desc(params.params_make_desc.input_path);
        for (const auto& segment : desc)
        {
            if (params.only_mol and segment.type != mdc::segment_type::molecule)
                continue;
            out << segment_type_t_to_str(segment.type) << "\t" << segment.name << "\t" << segment.size << "\n";
        }
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what();
        return false;
    }
    catch (...)
    {
        return false;
    }
    return true;
}

// *******************************************************************************
bool CApplication::info()
{
    CMDCDecompress mdc_decompress(params.preset_values.anchor_separation(), params.max_history_size());
    refresh::archive_input ar_in;

    if (!ar_in.open_file_buffered(params.input_fn, false, 8 << 20))
    {
        std::cerr << "Cannot open: " << params.input_fn << std::endl;
        return false;
    }

    auto topology_info = read_topology_info(ar_in);
    if (not topology_info)
        return false;

    packed_t packed;
    uint64_t metadata;

    mdc_file_version file_version;
    file_version.load(ar_in);

    std::cerr << "mdc format ver. " << file_version.format_major << "." << file_version.format_minor << std::endl;
    std::cerr << "produced by mdc v. " << file_version.produced_by_major << "." << file_version.produced_by_minor << "." << file_version.produced_by_patch << std::endl;

    if (not topology_info->IsTopologyStored())
    {
        std::cerr << "topology was not stored during compression\n";
    }
    else
    {
        std::cerr << "topology was stored during compression, original topology file name: " << topology_info->GetFilename() << "\n";
    }

    int id = ar_in.get_stream_id("segment-desc");

    if (id < 0)
    {
        std::cerr << "Corrupted input file" << std::endl;
        return false;
    }

    ar_in.get_part(id, packed, metadata);

    if (!mdc_decompress.DeserializeSegmentDesc(packed))
    {
        std::cerr << "Corrupted input file" << std::endl;
        return false;
    }

    // to have no frames
    id = ar_in.get_stream_id("anchor-ids");

    if (id < 0)
    {
        std::cerr << "Corrupted input file" << std::endl;
        return false;
    }

    ar_in.get_part(id, packed, metadata);

    if (!mdc_decompress.DeserializeAnchorDesc(packed))
    {
        std::cerr << "Corrupted input file" << std::endl;
        return false;
    }

    const auto& segment_description = mdc_decompress.GetSegmentDescription();
    const auto& subsegment_description = mdc_decompress.GetSubSegmentDescription();

    std::cerr << "no. frames: " << mdc_decompress.GetNoFrames() << std::endl;
    std::cerr << "no. atoms: "  << std::accumulate(segment_description.begin(),
                                                    segment_description.end(),
                                                    0, [](uint32_t sum, const segment_desc_t& desc) {
                                                        return sum + desc.size;
                                                    }) << std::endl;
    uint32_t subsegment_i{};
    for (uint32_t i = 0; i < segment_description.size(); ++i)
    {
        std::cerr << i << " - " << segment_description[i].name << " - ";
        std::cerr << segment_type_t_to_str(segment_description[i].type) << " - ";
        std::cerr << segment_description[i].size << std::endl;
        if (params.full_info)
        {
            do {
                std::cerr << "    " << subsegment_i << " - " << subsegment_description[subsegment_i].name << " - ";
                std::cerr << segment_type_t_to_str(subsegment_description[subsegment_i].type) << " - ";
                std::cerr << subsegment_description[subsegment_i].size << std::endl;
                ++subsegment_i;
            } while (subsegment_i < subsegment_description.size() && subsegment_description[subsegment_i].name.empty());
        }
    }

    return true;
}

// *******************************************************************************
bool CApplication::to_fp32()
{
    auto t1 = std::chrono::high_resolution_clock::now();

    std::unique_ptr<InputReader> input_reader;

    auto traj_format = get_traj_format(params.input_fn, false);

#if defined(USE_CHEMFILES) && defined(USE_XDRFILE)
    if (params.native_lib_mode)
    {
        switch (traj_format)
        {
        case traj_file_format_t::trr:
            input_reader = std::make_unique<CTRRReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
        case traj_file_format_t::xtc:
            input_reader = std::make_unique<CXTCReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
        default:
            return false;
        }
    }
    else
        input_reader = std::make_unique<CTrajReader>(params.resolution_set_by_user ? params.resolution() : 0);
#elif defined(USE_CHEMFILES)
    input_reader = std::make_unique<CTrajReader>(params.resolution_set_by_user ? params.resolution() : 0);
#else
    switch (traj_format)
    {
    case traj_file_format_t::trr:
        input_reader = std::make_unique<CTRRReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
    case traj_file_format_t::xtc:
        input_reader = std::make_unique<CXTCReader>(params.resolution_set_by_user ? params.resolution() : 0); break;
    default:
        return false;
    }
#endif

    FILE *f_fp32;
    refresh::archive_output ar_out;

    if (!load_segment_description_file())
        return false;

    if (!input_reader->open(params.input_fn, segment_desc))
    {
        std::cerr << "Cannot open: " << params.input_fn << std::endl;
        return false;
    }

    f_fp32 = fopen(params.output_fn.c_str(), "wb");
    if (!f_fp32)
    {
        std::cerr << "Cannot open: " << params.output_fn << std::endl;
        return false;
    }

    frame_t frame;

    while (input_reader->get_frame(frame))
    {
        std::cerr << "Frame time: " << frame.desc.time << "\tNo. atoms: " << frame.desc.n_atoms << std::endl;

        for (uint32_t i = 0; i < (uint32_t)frame.segments.size(); ++i)
        {
            for (auto& coords : frame.segments[i].coords)
            {
                float fcoord[3] = { (float)coords.x / frame.desc.prec, (float)coords.y / frame.desc.prec , (float)coords.z / frame.desc.prec };
                fwrite((char*) fcoord, 4, 3, f_fp32);
            }
        }
    }

    fclose(f_fp32);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = t2 - t1;

    std::cerr << "Compressing time: " << diff.count() << "s\n";

    return true;
}
