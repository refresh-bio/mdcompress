#pragma once

#ifdef __cplusplus

#include <string>
#include <vector>
#include <cinttypes>
#include <memory>
#include <span>

//if defined we try to follow chemfiles calculations (order of floating point aithmetic and double <-> float conversions) to produce exactly the same TRR
//In general the use case is: we have XTC file and we convert it to TRR. When this is doce by chemfiles (just reading and writing) and by gmx the results are identical.
//If we do this XTC -> MDC -> TRR we get slighly different TRR, but we believe this is acceptable, because the TRR from XTC isn't perfect as XTC looses some precision, so TRR from it is an approximation, and we beliewe our approximation is as good
//In general to have the same exact TRR as gmx we can do a couple of things
// 1. Change representation from angstrems to nm in our API - this is not implemented here, in such a case we just avoid one multiplication inside API by 10 in part of code available when following define is uncommented
// 2. The other approach requires changes not only in this code but also in chemfiles, i.e. not divide by 10 before storying - this is also not implemented here, but as the above I have tested and it seems to work, at least on one large TRR file I selected for tests (1xws.0.1ps.xtc)
// 3. This option is turned on when the define below is uncommented, unfortunatelly it changes the representation in API from floats to doubles (this breaks some bindings)
// In general this should not be used, I'm leaving it if this topic arise in the future to have a starting point
//#define FOLLOW_CHEMFILES

#else
#include <stdint.h>
#include <stddef.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef FOLLOW_CHEMFILES
    struct mdc_atom_coords
    {
        double x;
        double y;
        double z;
#ifdef __cplusplus
        mdc_atom_coords(double x, double y, double z) :
            x(x),
            y(y),
            z(z)
        {
        }
        bool operator==(const mdc_atom_coords& rhs) const
        {
            return x == rhs.x && y == rhs.y && z == rhs.z;
        }
#endif
    };
#else
    struct mdc_atom_coords
    {
        float x;
        float y;
        float z;
#ifdef __cplusplus
		mdc_atom_coords(float x, float y, float z) :
			x(x),
			y(y),
			z(z)
		{}
        bool operator==(const mdc_atom_coords& rhs) const
        {
            return x == rhs.x && y == rhs.y && z == rhs.z;
        }
#endif
    };
#endif

#ifndef __cplusplus
	typedef struct mdc_atom_coords mdc_atom_coords;
#endif

    //int values of this types are used in compressed archive
    //so for backward compatibility if this needs to be extended we should add new entries at the end
    typedef enum mdc_segment_type {
        MD_COMPRESS_SEGMENT_TYPE_UNKNOWN = 0,
        MD_COMPRESS_SEGMENT_TYPE_MOLECULE = 1,
        MD_COMPRESS_SEGMENT_TYPE_OTHER = 2,
        MD_COMPRESS_SEGMENT_TYPE_WATER = 3,
        MD_COMPRESS_SEGMENT_TYPE_NONE = 4
    } mdc_segment_type;

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
namespace mdc
{
    //!!!keep in sync with mdc_segment_type!!!
    enum class segment_type
    {
	    unknown = MD_COMPRESS_SEGMENT_TYPE_UNKNOWN,
        molecule = MD_COMPRESS_SEGMENT_TYPE_MOLECULE,
        other = MD_COMPRESS_SEGMENT_TYPE_OTHER,
    	water = MD_COMPRESS_SEGMENT_TYPE_WATER,
    	none = MD_COMPRESS_SEGMENT_TYPE_NONE
    };

    struct segment_desc_t
    {
        segment_type type;
        std::string name;
        uint32_t size;

        segment_desc_t() :
            type(segment_type::unknown),
            name(""),
            size(0)
        {}

        segment_desc_t(const segment_type type, const std::string& name, const uint32_t size) :
            type(type),
            name(name),
            size(size)
        {}

        segment_desc_t(const segment_desc_t&) = default;
        segment_desc_t(segment_desc_t&&) = default;

        segment_desc_t& operator=(const segment_desc_t&) = default;
        segment_desc_t& operator=(segment_desc_t&&) = default;
    };

    using atom_coords = mdc_atom_coords;
    
	//*************************************
    //query_result from which `frame` cames from must be kept alivie for `frame` to be usable
    struct frame
    {
        std::span<atom_coords> coords{}; //span in appropriate frame coords in all_frames_coords
        
        int step;
        float time;
        double box[3][3];
        float prec;
        float lambda{};
        int has_prop{};
    };

    struct query_result
    {
        //have its own copy ctor and assignment operator to handle span correctly
    private:
        friend class query_engine_impl;
        std::vector<atom_coords> all_frames_coords;
    public:
		std::vector<frame> frames;
        query_result() = default;
        query_result(const query_result& rhs);
        query_result& operator=(const query_result& rhs);
        query_result(query_result&&) = default;
        query_result& operator=(query_result&&) = default;
    };
    //*************************************

	//this class is not thread safe, and each thread should have its indivitual instance
    class query_engine
    {
        std::unique_ptr<class query_engine_impl> pImpl;

        friend class reader;
        //this may be created by mdc_reader only
        query_engine(class reader_impl* reader,
            const std::vector<std::string>& segment_ids,
            const std::vector<uint32_t>& atom_ids);

    public:
        query_engine(const query_engine& rhs) = delete;
        query_engine& operator=(const query_engine& rhs) = delete;
        query_engine(query_engine&& rhs) noexcept;
        query_engine& operator=(query_engine&& rhs) noexcept;

        const std::vector<uint32_t>& get_original_atom_ids() const;

        //mkokot_TODO: this probably should be marked const, and also in impl I should have const and use mutable for internal state (because it is all related to cacheing and does not modify the logical state?)
        bool query(std::span<const uint32_t> frame_ids, query_result& result);

        //mkokot_TODO: this is just a proxy to allow call like query({ 2 }, result);
        //consider if we need this...
        bool query(std::initializer_list<uint32_t> frame_ids, query_result& result) {
	        return query(std::span(frame_ids.begin(), frame_ids.end()), result);
        }

        //if returns empty string, then no error
        //in the oposite case returns last error message
        //errors are never cleared and are considered unrecoverable
        const std::string& get_current_error() const;

    	~query_engine();
    };

	//********************************************
    
    class reader
    {
        std::unique_ptr<class reader_impl> pImpl;
    public:
        reader(const std::string& path);
        
        reader(const reader& rhs) = delete;
        reader& operator=(const reader& rhs) = delete;
        reader(reader&& rhs) noexcept;
        reader& operator=(reader&& rhs) noexcept;

        //if returns empty string, then no error
        //in the oposite case returns last error message
        //errors are never cleared and are considered unrecoverable
        const std::string& get_current_error() const;

        const std::vector<segment_desc_t>& get_segments() const;

        uint32_t get_no_frames() const;

        const std::vector<uint32_t>& get_anchor_ids() const;

        std::unique_ptr<query_engine> get_query_engine(
            const std::vector<std::string>& segment_ids,
            const std::vector<uint32_t>& atom_ids) const;

        ~reader();
    };
}

extern "C" {
#endif
    typedef struct mdc_reader mdc_reader;
    typedef struct mdc_query_engine mdc_query_engine;
    typedef struct mdc_query_result_impl mdc_query_result_impl;
//
    typedef struct mdc_segment_desc {
        mdc_segment_type type;
        const char* name;
        uint32_t size;
    }mdc_segment_desc;


    typedef struct mdc_frame
    {
        mdc_atom_coords* coords;
        uint32_t n_coords;
        int step;
        float time;
        double box[3][3];
        float prec;
    }mdc_frame;

     //may be reused in subsequent queries (and its recommended)
    typedef struct mdc_query_result
    {
        mdc_frame* frames;
        uint32_t n_frames;

        mdc_query_result_impl* impl;
    } mdc_query_result;

    /*** mdc_reader ***/
    mdc_reader* mdc_reader_open(const char* path);

    const char* mdc_reader_get_error(mdc_reader* reader);

    uint32_t mdc_get_no_segments(mdc_reader* reader);

    mdc_segment_desc mdc_get_segment_desc(mdc_reader* reader, uint32_t index);

    uint32_t mdc_get_no_anchors(mdc_reader* reader);

    const uint32_t* mdc_get_anchor_ids(mdc_reader* reader);

    uint32_t mdc_get_no_frames(mdc_reader* reader);

    void mdc_reader_close(mdc_reader* reader);

    mdc_query_engine* mdc_get_query_engine(mdc_reader* reader, char** segment_ids, size_t num_segment_ids, uint32_t* atom_ids, size_t num_atom_ids);
    //*******************************************************

    /*** mdc_query_engine ***/
    const char* mdc_query_engine_get_error(mdc_query_engine* engine);

    uint32_t mdc_query_engine_get_no_original_atom_ids(mdc_query_engine* engine);

    const uint32_t* mdc_query_engine_get_original_atom_ids(mdc_query_engine* engine);

    int mdc_query(mdc_query_engine* engine, uint32_t* frame_ids, uint32_t n_frame_ids, mdc_query_result* result);

    void mdc_free_query_engine(mdc_query_engine* engine);
    //*******************************************************

    /*** mdc_query_result ***/
    mdc_query_result* mdc_create_query_result();

    void mdc_free_query_result(mdc_query_result* result);
    //*******************************************************
//
#ifdef __cplusplus
}
#endif
