#pragma once

#define MDCOMPRESS_VERSION_MAJOR 2
#define MDCOMPRESS_VERSION_MINOR 0
#define MDCOMPRESS_VERSION_PATCH 1

#define MDCOMPRESS_DATE "2026-05-27"


inline void MDCOMPRESS_VER_PRINT(std::ostream& oss) {
	oss << "mdcompress version: " <<
		MDCOMPRESS_VERSION_MAJOR << "." <<
		MDCOMPRESS_VERSION_MINOR << "." <<
		MDCOMPRESS_VERSION_PATCH << "\n";
	oss << "Date: " << MDCOMPRESS_DATE << "\n";
	oss << "Authors: Marek Kokot, Sebastian Deorowicz\n";
}
