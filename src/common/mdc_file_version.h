#pragma once
#include "serializer.h"
#include "../libs/refresh/archive/lib/archive_input.h"
#include "../libs/refresh/archive/lib/archive_output.h"
#include "version.h"
#include <vector>

struct mdc_file_version
{
	//currently use soft major and minor as format version
	int format_major = MDCOMPRESS_VERSION_MAJOR;
	int format_minor = MDCOMPRESS_VERSION_MINOR;

	int produced_by_major = MDCOMPRESS_VERSION_MAJOR;
	int produced_by_minor = MDCOMPRESS_VERSION_MINOR;
	int produced_by_patch = MDCOMPRESS_VERSION_PATCH;

	//currently require this but this may change if format versioning will differ from software versioning
	static constexpr int MAX_SUPPORTED_FORMAT_MAJOR = 1;

	bool load(refresh::archive_input& ar_in)
	{
		std::vector<uint8_t> packed;
		uint64_t metadata;

		int id = ar_in.get_stream_id("format-ver");
		if (id == -1)
			return false;

		if (!ar_in.get_part(id, 0, packed, metadata))
			return false;

		Serializer serializer;
		serializer.set_data(packed);
		format_major = serializer.load32i();
		format_minor = serializer.load32i();

		id = ar_in.get_stream_id("soft-ver");
		if (id == -1)
			return false;

		if (!ar_in.get_part(id, 0, packed, metadata))
			return false;

		serializer.set_data(packed);
		produced_by_major = serializer.load32i();
		produced_by_minor = serializer.load32i();
		produced_by_patch = serializer.load32i();

		return true;
	}

	bool is_supported() const
	{
		return format_major <= MAX_SUPPORTED_FORMAT_MAJOR;
	}

	bool store(refresh::archive_output& ar_out) const
	{
		Serializer serializer;
		serializer.append32i(format_major);
		serializer.append32i(format_minor);
		ar_out.add_part(ar_out.register_stream("format-ver"), serializer.get_data());

		serializer.clear();
		serializer.append32i(produced_by_major);
		serializer.append32i(produced_by_minor);
		serializer.append32i(produced_by_patch);
		ar_out.add_part(ar_out.register_stream("soft-ver"), serializer.get_data());

		return true;
	}
};