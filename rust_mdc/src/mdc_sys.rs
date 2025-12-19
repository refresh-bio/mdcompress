#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int, c_float};

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct mdc_atom_coords {
    pub x: c_float,
    pub y: c_float,
    pub z: c_float,
}

#[repr(i32)] // or should it be repr(C)?
#[derive(Debug, Copy, Clone)]
pub enum mdc_segment_type {
    MD_COMPRESS_SEGMENT_TYPE_UNKNOWN = 0,
    MD_COMPRESS_SEGMENT_TYPE_MOLECULE,
    MD_COMPRESS_SEGMENT_TYPE_OTHER,
    MD_COMPRESS_SEGMENT_TYPE_WATER,
    MD_COMPRESS_SEGMENT_TYPE_NONE,
}

#[repr(C)]
#[derive(Debug)]
pub struct mdc_segment_desc {
    pub type_: mdc_segment_type,
    pub name: *const c_char,
    pub size: u32,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct mdc_frame {
    pub coords: *mut mdc_atom_coords,
    pub n_coords: u32,
    pub step: c_int,
    pub time: c_float,
    pub box_: [[c_float; 3]; 3],
    pub prec: c_float,
}

#[repr(C)]
pub struct mdc_query_result_impl {
    _private: [u8; 0],
}


#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct mdc_query_result {
    pub frames: *mut mdc_frame,
    pub n_frames: u32,

    _impl: *mut mdc_query_result_impl,
}

// Opaque structs from the C API
#[repr(C)]
pub struct mdc_reader {
    _private: [u8; 0],
}

#[repr(C)]
pub struct mdc_query_engine {
    _private: [u8; 0],
}

extern "C" {
    // -- mdc_reader --
    pub fn mdc_reader_open(path: *const c_char) -> *mut mdc_reader;

    pub fn mdc_reader_get_error(reader: *mut mdc_reader) -> *const c_char;

    pub fn mdc_get_no_segments(reader: *mut mdc_reader) -> u32;

    pub fn mdc_get_segment_desc(reader: *mut mdc_reader, index: u32) -> mdc_segment_desc;

    pub fn mdc_get_no_anchors(reader: *mut mdc_reader) -> u32;

    pub fn mdc_get_anchor_ids(reader: *mut mdc_reader) -> *const u32;

    pub fn mdc_get_no_frames(reader: *mut mdc_reader) -> u32;
    
    pub fn mdc_reader_close(reader: *mut mdc_reader);

    pub fn mdc_get_query_engine(reader: *mut mdc_reader, segment_ids: *const *const c_char, num_segment_ids : usize, atom_ids: *const u32, num_atom_ids: usize) -> * mut mdc_query_engine;


    // -- mdc_query_engine --
    pub fn mdc_query_engine_get_error(engine: *mut mdc_query_engine) -> *const c_char;

    pub fn mdc_query_engine_get_no_original_atom_ids(engine: *mut mdc_query_engine) -> u32;

    pub fn mdc_query_engine_get_original_atom_ids(engine: *mut mdc_query_engine) ->  *const u32;

    pub fn mdc_query(engine: *mut mdc_query_engine, frame_ids: *const u32, n_frame_ids: u32, result: *mut mdc_query_result) -> c_int;

    pub fn mdc_free_query_engine(engine: *mut mdc_query_engine);

    // -- mdc_query_result --
    pub fn mdc_create_query_result() -> *mut mdc_query_result;

    pub fn mdc_free_query_result(result: *mut mdc_query_result);
}
