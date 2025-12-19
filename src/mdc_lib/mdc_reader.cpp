#include "mdc_reader.h"
#include <refresh/archive/lib/archive_input.h>
#include "../common/mdc_decompress.h"
#include "../common/mdc_file_version.h"
#include "../common/version.h"

namespace mdc
{
	query_result::query_result(const query_result& rhs) : all_frames_coords(rhs.all_frames_coords)
	{
		frames.reserve(rhs.frames.size());
		for (const auto& f : rhs.frames)
		{
			frames.push_back(f);
			frames.back().coords = std::span(std::span(all_frames_coords.data() + (f.coords.data() - rhs.all_frames_coords.data()), f.coords.size()));
		}
	}

	query_result& query_result::operator=(const query_result& rhs)
	{
		if (this != &rhs)
		{
			all_frames_coords = rhs.all_frames_coords;
			frames.clear();
			frames.reserve(rhs.frames.size());
			for (const auto& f : rhs.frames)
			{
				frames.push_back(f);
				frames.back().coords = std::span(std::span(all_frames_coords.data() + (f.coords.data() - rhs.all_frames_coords.data()), f.coords.size()));
			}
		}
		return *this;
	}

	class query_engine_impl
	{
		CMDCDecompress mdc_decompress;

		refresh::archive_input* ar_in;
		std::vector<int> streams_segments_ids; //set only for selected segments, for other segments it is -1
		int stream_frame_id{ -1 };

		frame_t out_frame; //class field to avoid reallocations

		//for cacheing query state:
		CMDCDecompress::packed_frames_t packed_frames;
		uint32_t prev_anchor_idx = std::numeric_limits<uint32_t>::max();
		uint32_t prev_frame_offset;

		//mkokot_TODO: copied from CApplication -> refactor
		std::pair<uint32_t, uint32_t> locate_frame(const std::vector<uint32_t>& anchor_ids, const uint32_t frame_id)
		{
			auto p = std::upper_bound(anchor_ids.begin(), anchor_ids.end(), frame_id) - 1;
			return std::make_pair(p - anchor_ids.begin(), frame_id - *p);
		}

		std::string current_error_message;
	public:
		query_engine_impl(
			refresh::archive_input* ar_in,
			const std::vector<std::string>& segment_ids,
			const std::vector<uint32_t>& atom_ids) :
			mdc_decompress(0, 1),
			ar_in(ar_in)
		{
			packed_t packed;
			uint64_t metadata;

			int id = ar_in->get_stream_id("segment-desc");

			if (id < 0)
			{
				current_error_message = "Corrupted input file (cannot load segment-desc)";
			}

			ar_in->get_part(id, 0, packed, metadata);

			if (!mdc_decompress.DeserializeSegmentDesc(packed))
			{
				current_error_message = "Corrupted input file (cannot deserialize segment desc)";
				return;
			}

			if (!mdc_decompress.SelectSegmentsAndAtoms(segment_ids, atom_ids))
			{
				current_error_message = "Error in segment and atoms selection";
				return;
			}

			id = ar_in->get_stream_id("anchor-ids");

			if (id < 0)
			{
				current_error_message = "Corrupted input file (cannot load anchor-ids)";
				return;
			}

			ar_in->get_part(id, 0, packed, metadata);

			if (!mdc_decompress.DeserializeAnchorDesc(packed))
			{
				current_error_message = "Corrupted input file (cannot deserialize anchor desc)";
				return;
			}

			stream_frame_id = ar_in->get_stream_id("frame-id");
			if (stream_frame_id == -1)
			{
				current_error_message = "Cannot find stream frame-id";
				return;
			}

			streams_segments_ids.assign(mdc_decompress.GetNoSubSegments(), -1);
			for (uint32_t i = 0; i < mdc_decompress.GetNoSubSegments(); ++i)
			{
				if (mdc_decompress.IsSubsegmentSelected(i))
				{
					streams_segments_ids[i] = ar_in->get_stream_id("stream-" + std::to_string(i));
					if (streams_segments_ids[i] == -1)
					{
						current_error_message = "Cannot find stream stream-" + std::to_string(i);
						return;
					}
					const int model_id = ar_in->get_stream_id("model-" + std::to_string(i));
					if (model_id == -1)
					{
						current_error_message = "Cannot find stream model-" + std::to_string(i);
						return;
					}
					if (!ar_in->get_part(model_id, 0, packed, metadata))
					{
						current_error_message = "Cannot read data of model-" + std::to_string(i);
						return;
					}
					if (!mdc_decompress.DeserializeModel(i, packed))
					{
						current_error_message = "Cannot deserialize model for " + std::to_string(i) + "-th subsegment";
						return;
					}
				}
			}

			if (!mdc_decompress.StartDecompression())
			{
				current_error_message = "Cannot start decompression";
				return;
			}

			//set initial state of query state
			prev_anchor_idx = std::numeric_limits<uint32_t>::max();
			prev_frame_offset = std::numeric_limits<uint32_t>::max();
			packed_frames.packed.clear();
			packed_frames.packed.resize(mdc_decompress.GetNoSubSegments(), packed_t{});
		}
		query_engine_impl(const query_engine_impl& rhs) = delete;
		query_engine_impl& operator=(const query_engine_impl& rhs) = delete;
		query_engine_impl(query_engine_impl&& rhs) = delete;
		query_engine_impl& operator=(query_engine_impl&& rhs) = delete;

		const std::vector<uint32_t>& get_original_atom_ids() const
		{
			return mdc_decompress.GetOriginalAtomIds();
		}

		bool query(std::span<const uint32_t> frame_ids, query_result& result)
		{
			if (frame_ids.empty())
			{
				current_error_message = "frame_ids is empty";
				return false;
			}

			result.all_frames_coords.clear();
			result.frames.resize(frame_ids.size());

			std::vector<uint32_t> n_atoms_per_frame(result.frames.size());

			for (uint32_t frame_idx = 0 ; frame_idx < frame_ids.size() ; ++frame_idx)
			{
				uint32_t frame_id = frame_ids[frame_idx];
				if (frame_id >= mdc_decompress.GetNoFrames())
				{
					current_error_message = "frame_id (" + std::to_string(frame_id) + ") is >= number of frames (" + std::to_string(mdc_decompress.GetNoFrames()) + ")";
					return false;
				}

				//mkokot_TODO: maybe we could avoild locating frames at each query by keeping the start and end of current anchor
				auto [anchor_idx, frame_offset] = locate_frame(mdc_decompress.GetAnchorIds(), frame_id);

				//query for new anchor?
				if (prev_anchor_idx != anchor_idx)
				{
					prev_anchor_idx = anchor_idx;

					//force start reading from anchor_idx
					prev_frame_offset = std::numeric_limits<uint32_t>::max();
				}

				if (prev_frame_offset > frame_offset)
				{
					uint64_t metadata;
					for (uint32_t subsegment_id = 0; subsegment_id < mdc_decompress.GetNoSubSegments(); ++subsegment_id)
					{
						if (mdc_decompress.IsSubsegmentSelected(subsegment_id))
						{
							if (!ar_in->get_part(streams_segments_ids[subsegment_id], anchor_idx, packed_frames.packed[subsegment_id], metadata))
							{
								current_error_message = "Cannot read data for subsegment id " + std::to_string(subsegment_id) + ", frame_id = " + std::to_string(frame_id) + ", anchor_idx = " + std::to_string(anchor_idx);
								return false;
							}
						}
					}
					if (!ar_in->get_part(stream_frame_id, anchor_idx, packed_frames.desc, metadata))
					{
						current_error_message = "Cannot read description for frame_id = " + std::to_string(frame_id) + ", anchor_idx = " + std::to_string(anchor_idx);
						return false;
					}

					if (!mdc_decompress.SetPackedFrames(packed_frames))
					{
						current_error_message = "Cannot set packed frames";
						return false;
					}

					prev_frame_offset = 0;
				}

				for (; prev_frame_offset < frame_offset + 1; ++prev_frame_offset)
				{
					if (!mdc_decompress.DecompressFrame(out_frame))
					{
						current_error_message = "Cannot decompress frame";
						return false;
					}
				}

				//the last decompressed frame was the one asked for...

				result.frames[frame_idx].step = out_frame.desc.step;
				result.frames[frame_idx].time = out_frame.desc.time;
				result.frames[frame_idx].prec = out_frame.desc.prec;

				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						result.frames[frame_idx].box[i][j] = out_frame.desc.box[i][j];

				auto prev_coords_size = result.all_frames_coords.size();

#ifdef FOLLOW_CHEMFILES
				float inv_precision = 1.0f / (result.frames[frame_idx].prec * 10);
#else
				float inv_precision = 1.0f / result.frames[frame_idx].prec;
#endif
				for (const auto& out_segment : out_frame.segments)
					for (const auto& coord : out_segment.coords)
						result.all_frames_coords.emplace_back(
#ifdef FOLLOW_CHEMFILES
							coord.x * inv_precision * 10.0,
							coord.y * inv_precision * 10.0,
							coord.z * inv_precision * 10.0
#else
							coord.x * inv_precision,
							coord.y * inv_precision,
							coord.z * inv_precision
#endif
						);
				n_atoms_per_frame[frame_idx] = result.all_frames_coords.size() - prev_coords_size;
			}

			auto ptr = result.all_frames_coords.data();

			for (uint32_t frame_idx = 0; frame_idx < frame_ids.size(); ++frame_idx) {
				result.frames[frame_idx].coords = std::span(ptr, n_atoms_per_frame[frame_idx]);
				ptr += n_atoms_per_frame[frame_idx];
			}

			return true;
		}

		//if returns empty string, then no error
		//in the oposite case returns last error message
		//errors are never cleared and are considered unrecoverable
		const std::string& get_current_error() const
		{
			return current_error_message;
		}

		~query_engine_impl()
		{
		}
	};

	class reader_impl
	{
		friend class query_engine;
		refresh::archive_input ar_in;

		//mkokot_INFO:
		//its here only for access to segment ids + no of frames + anchor ids
		//in fact it would be better to read these info here and pass to each query_engine
		//but it would require changes in CMCDecompress so lets leave it this way for now
		CMDCDecompress mdc_decompress;

		std::string current_error_message;
	public:
		reader_impl(const std::string& path):
			mdc_decompress(0, 1)
		{
			//mkokot_TODO: consider buffering! harmful for parallel decompression
			//if (!ar_in.open_file_buffered(path, false, 8 << 20))
			//if (!ar_in.open_file_buffered(path, false, 8 << 12))
			if (!ar_in.open_file_unbuffered(path, false))
			{
				current_error_message = "Cannot open: " + path;
				return;
			}
			mdc_file_version file_version;
			if (!file_version.load(ar_in))
			{
				current_error_message = "Cannot load mdc version from file " + path;
				return;
			}

			if (!file_version.is_supported())
			{
				current_error_message = "mdc file format is newer than this software - please upade software";
				return;
			}

			int id = ar_in.get_stream_id("segment-desc");

			if (id < 0)
			{
				current_error_message = "Corrupted input file (cannot load segment-desc)";
			}

			packed_t packed;
			uint64_t metadata;

			ar_in.get_part(id, 0, packed, metadata);

			if (!mdc_decompress.DeserializeSegmentDesc(packed))
			{
				current_error_message = "Corrupted input file (cannot deserialize segment desc)";
				return;
			}

			id = ar_in.get_stream_id("anchor-ids");

			if (id < 0)
			{
				current_error_message = "Corrupted input file (cannot load anchor-ids)";
				return;
			}

			ar_in.get_part(id, 0, packed, metadata);

			if (!mdc_decompress.DeserializeAnchorDesc(packed))
			{
				current_error_message = "Corrupted input file (cannot deserialize anchor desc)";
				return;
			}
		}

		const std::vector<segment_desc_t>& get_segments() const
		{
			return mdc_decompress.GetSegmentDescription();
		}

		uint32_t get_no_frames() const
		{
			return mdc_decompress.GetNoFrames();
		}

		const std::vector<uint32_t>& get_anchor_ids() const
		{
			return mdc_decompress.GetAnchorIds();
		}

		reader_impl(const reader_impl& rhs) = delete;
		reader_impl& operator=(const reader_impl& rhs) = delete;
		reader_impl(reader_impl&& rhs) = delete;
		reader_impl& operator=(reader_impl&& rhs) = delete;

		//if returns empty string, then no error
		//in the oposite case returns last error message
		//errors are never cleared and are considered unrecoverable
		const std::string& get_current_error() const
		{
			return current_error_message;
		}
	};

	query_engine::query_engine(
		reader_impl* reader,
		const std::vector<std::string>& segment_ids,
		const std::vector<uint32_t>& atom_ids) :
		pImpl(std::make_unique<query_engine_impl>(&reader->ar_in, segment_ids, atom_ids))
	{

	}

	query_engine::query_engine(query_engine&& rhs) noexcept = default;
	query_engine& query_engine::operator=(query_engine&& rhs) noexcept = default;

	const std::vector<uint32_t>& query_engine::get_original_atom_ids() const
	{
		return pImpl->get_original_atom_ids();
	}

	bool query_engine::query(std::span<const uint32_t> frame_ids, query_result& result)
	{
		return pImpl->query(frame_ids, result);
	}

	const std::string& query_engine::get_current_error() const
	{
		return pImpl->get_current_error();
	}

	query_engine::~query_engine() = default;

	reader::reader(const std::string& path) :
			pImpl(std::make_unique<reader_impl>(path))
	{
		
	}

	reader::reader(reader&& rhs) noexcept = default;
	reader& reader::operator=(reader&& rhs) noexcept = default;

	const std::vector<segment_desc_t>& reader::get_segments() const
	{
		return pImpl->get_segments();
	}

	uint32_t reader::get_no_frames() const
	{
		return pImpl->get_no_frames();
	}

	const std::vector<uint32_t>& reader::get_anchor_ids() const
	{
		return pImpl->get_anchor_ids();
	}

	std::unique_ptr<query_engine> reader::get_query_engine(
		const std::vector<std::string>& segment_ids,
		const std::vector<uint32_t>& atom_ids) const
	{
		return std::unique_ptr<query_engine>(new query_engine(pImpl.get(), segment_ids, atom_ids)); //cannot make_unique because of private ctor
	}

	const std::string& reader::get_current_error() const
	{
		return pImpl->get_current_error();
	}
	reader::~reader() = default;
}



namespace mdc {
	//check if can cast enum, if extended, this check should also be extended
	static_assert((int)MD_COMPRESS_SEGMENT_TYPE_UNKNOWN == (int)segment_type_t::unknown, "Enum mismatch");
	static_assert((int)MD_COMPRESS_SEGMENT_TYPE_MOLECULE == (int)segment_type_t::molecule, "Enum mismatch");
	static_assert((int)MD_COMPRESS_SEGMENT_TYPE_OTHER == (int)segment_type_t::other, "Enum mismatch");
	static_assert((int)MD_COMPRESS_SEGMENT_TYPE_WATER == (int)segment_type_t::water, "Enum mismatch");
	static_assert((int)MD_COMPRESS_SEGMENT_TYPE_NONE == (int)segment_type_t::none, "Enum mismatch");
}

extern "C" {

	mdc_reader* mdc_reader_open(const char* path) {
		return (mdc_reader*)new mdc::reader(path);
	}

	const char* mdc_reader_get_error(mdc_reader* reader) {
		const auto& err = ((mdc::reader*)reader)->get_current_error();
		if (err.empty())
			return nullptr;
		return err.c_str();
	}

	uint32_t mdc_get_no_segments(mdc_reader* reader) {
		return ((mdc::reader*)reader)->get_segments().size();
	}

	mdc_segment_desc mdc_get_segment_desc(mdc_reader* reader, uint32_t index) {
		mdc_segment_desc ret;
		auto& seg = ((mdc::reader*)reader)->get_segments()[index];
		ret.name = seg.name.c_str();
		ret.type = (mdc_segment_type)seg.type;
		ret.size = seg.size;
		return ret;
	}

	uint32_t mdc_get_no_anchors(mdc_reader* reader) {
		return (uint32_t)((mdc::reader*)reader)->get_anchor_ids().size();
	}

	const uint32_t* mdc_get_anchor_ids(mdc_reader* reader) {
		return ((mdc::reader*)reader)->get_anchor_ids().data();
	}

	uint32_t mdc_get_no_frames(mdc_reader* reader) {
		return ((mdc::reader*)reader)->get_no_frames();
	}

	void mdc_reader_close(mdc_reader* reader) {
		if (reader)
			delete (mdc::reader*)reader;
	}

	mdc_query_engine* mdc_get_query_engine(mdc_reader* reader, char** segment_ids, size_t num_segment_ids, uint32_t* atom_ids, size_t num_atom_ids) {
		std::vector<std::string> segment_ids_vec(segment_ids, segment_ids + num_segment_ids);
		std::vector<uint32_t> atom_ids_vec(atom_ids, atom_ids + num_atom_ids);
		auto eng = ((mdc::reader*)reader)->get_query_engine(segment_ids_vec, atom_ids_vec);
		return (mdc_query_engine*)eng.release();
	}

	const char* mdc_query_engine_get_error(mdc_query_engine* engine) {
		const auto& err = ((mdc::query_engine*)engine)->get_current_error();
		if (err.empty())
			return nullptr;
		return err.c_str();
	}

	uint32_t mdc_query_engine_get_no_original_atom_ids(mdc_query_engine* engine) {
		return ((mdc::query_engine*)engine)->get_original_atom_ids().size();
	}

	const uint32_t* mdc_query_engine_get_original_atom_ids(mdc_query_engine* engine) {
		return ((mdc::query_engine*)engine)->get_original_atom_ids().data();
	}

	struct mdc_query_result_impl
	{
		mdc::query_result query_result;
		std::vector<mdc_frame> frames;
	};


	int mdc_query(mdc_query_engine* engine, uint32_t* frame_ids, uint32_t n_frame_ids, mdc_query_result* result) {
		auto res_impl = result->impl;

		auto res = ((mdc::query_engine*)engine)->query(std::span(frame_ids, n_frame_ids), res_impl->query_result);
		if (!res)
			return 1;

		//rewrite from C++ representation (that contains std::span which is not available in C) to C representation
		res_impl->frames.resize(res_impl->query_result.frames.size());
		for (uint32_t frame_id = 0; frame_id < res_impl->query_result.frames.size(); ++frame_id)
		{
			res_impl->frames[frame_id].coords = res_impl->query_result.frames[frame_id].coords.data();
			res_impl->frames[frame_id].n_coords = res_impl->query_result.frames[frame_id].coords.size();
			res_impl->frames[frame_id].step = res_impl->query_result.frames[frame_id].step;
			res_impl->frames[frame_id].time = res_impl->query_result.frames[frame_id].time;
			static_assert(sizeof(res_impl->frames[frame_id].box) == sizeof(res_impl->query_result.frames[frame_id].box));
			std::memcpy(res_impl->frames[frame_id].box, res_impl->query_result.frames[frame_id].box, sizeof(res_impl->frames[frame_id].box));
			res_impl->frames[frame_id].prec = res_impl->query_result.frames[frame_id].prec;
		}

		result->frames = res_impl->frames.data();
		result->n_frames = static_cast<uint32_t>(res_impl->query_result.frames.size());

		return 0;
	}

	void mdc_free_query_engine(mdc_query_engine* engine) {
		if (engine) {
			delete (mdc::query_engine*)engine;
		}
	}

	mdc_query_result* mdc_create_query_result() {
		auto ret = new mdc_query_result;
		ret->impl = new mdc_query_result_impl;
		return ret;
	}

	void mdc_free_query_result(mdc_query_result* result) {
		if (result) {
			delete result->impl;
			delete result;
		}
	}

}
