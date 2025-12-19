
use std::os::raw::{c_char};
use std::ffi::{CStr, CString};
use std::ptr;
use std::fmt;

use crate::mdc_sys::*;

pub type MdcAtomCoords = mdc_atom_coords;

//mkokot: consideration: should i have references to all fields or copy
//I have ref to coords to avoid large copies
#[derive(Debug, Clone)]
pub struct Frame<'a> {
    pub coords: &'a [MdcAtomCoords],
    pub step: i32,
    pub time: f32,
    pub box_: [[f32; 3]; 3],
    pub prec: f32,
}

impl<'a> Frame<'a> {
    pub fn new(inner: &'a mdc_frame) -> Self {
        Self { 
            coords: unsafe {
                std::slice::from_raw_parts(inner.coords, inner.n_coords as usize)
            },
            step: inner.step,
            time: inner.time,
            box_: inner.box_,
            prec: inner.prec,
        }
    }
}

#[derive(Debug, Clone)]
pub struct QueryResult<'a> {
    inner: *mut mdc_query_result,
    pub frames: Vec<Frame<'a>>,
}

impl<'a> QueryResult<'a> {
    pub fn new() -> Result<Self, String> {
        let inner = unsafe { mdc_create_query_result() };
        if inner.is_null() {
            return Err("Failed to create query result".to_string());
        }
        Ok(Self { inner, frames : Vec::new()})
    }
    
}

impl Drop for QueryResult<'_> {
        fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe { mdc_free_query_result(self.inner) };
            self.inner = ptr::null_mut();
        }
    }
}

pub struct Reader {
    inner: *mut mdc_reader,
}


#[derive(Debug, Clone, Copy)]
pub enum SegmentType {
    Unknown,
    Molecule,
    Other,
    Water,
    None_,
}

impl From<mdc_segment_type> for SegmentType {
    fn from(value: mdc_segment_type) -> Self {
        match value {
            mdc_segment_type::MD_COMPRESS_SEGMENT_TYPE_UNKNOWN => Self::Unknown,
            mdc_segment_type::MD_COMPRESS_SEGMENT_TYPE_MOLECULE => Self::Molecule,
            mdc_segment_type::MD_COMPRESS_SEGMENT_TYPE_OTHER => Self::Other,
            mdc_segment_type::MD_COMPRESS_SEGMENT_TYPE_WATER => Self::Water,
            mdc_segment_type::MD_COMPRESS_SEGMENT_TYPE_NONE => Self::None_
        }
    }
}

impl fmt::Display for SegmentType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let name = match self {
            SegmentType::Unknown => "Unknown",
            SegmentType::Molecule => "Molecule",
            SegmentType::Other => "Other",
            SegmentType::Water => "Water",
            SegmentType::None_ => "None",
        };
        write!(f, "{}", name)
    }
}

/// High-level Rust struct for segment description
#[derive(Debug, Clone)]
pub struct SegmentDesc {
    pub type_: SegmentType,
    pub name: String,
    pub size: u32,
}

impl From<mdc_segment_desc> for SegmentDesc {
    fn from(desc: mdc_segment_desc) -> Self {
        Self {
            type_: SegmentType::from(desc.type_),
            name: unsafe { CStr::from_ptr(desc.name) }.to_string_lossy().into_owned(),
            size: desc.size,
        }
    }
}


pub struct QueryEngine {
    inner: *mut mdc_query_engine,
}

impl QueryEngine {
    //private because it should only be created by Reader
    fn new(inner: *mut mdc_query_engine) -> Result<Self, String> {
        if inner.is_null() {
            return Err("Failed to create query engine".to_string());
        }

        let err_ptr = unsafe { mdc_query_engine_get_error(inner) };
        if !err_ptr.is_null() {
            let err = unsafe { CStr::from_ptr(err_ptr) };
            unsafe { mdc_free_query_engine(inner) };
            return Err(err.to_string_lossy().into_owned());
        }

        Ok(Self { inner })
    }

    pub fn get_original_atom_ids(&self) -> Result<&[u32], String> {
        let count = unsafe { mdc_query_engine_get_no_original_atom_ids(self.inner) };

        let ids = unsafe { mdc_query_engine_get_original_atom_ids(self.inner) };

        let err_ptr = unsafe { mdc_query_engine_get_error(self.inner) };

        if !err_ptr.is_null() {
            let err = unsafe { CStr::from_ptr(err_ptr) };
            return Err(err.to_string_lossy().into_owned());
        }

        Ok (unsafe { std::slice::from_raw_parts(ids, count as usize) })
    }

    pub fn query(
        &self, 
        frame_ids: &[u32],
        result: &mut QueryResult
    ) -> Result<(), String> {
        let status = unsafe {
            mdc_query(
                self.inner,
                frame_ids.as_ptr() as *mut u32,
                frame_ids.len() as u32,
                result.inner
            )
        };
        if status != 0 {
            let err_ptr = unsafe { mdc_query_engine_get_error(self.inner) };
            if !err_ptr.is_null() {
                let err = unsafe { CStr::from_ptr(err_ptr) };
                return Err(err.to_string_lossy().into_owned());
            }
            return Err("Failed to query".to_string());
        }
        let inner = unsafe { &*result.inner }; 
        let frames_slice = unsafe {
            std::slice::from_raw_parts(inner.frames, inner.n_frames as usize)
        };
        result.frames = frames_slice.iter().map(Frame::new).collect();
        Ok(())
    }
}

impl Drop for QueryEngine {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe { mdc_free_query_engine(self.inner) };
            self.inner = ptr::null_mut();
        }
    }
}

impl Reader {
    pub fn new(path: &str) -> Result<Self, String> {
        let c_path = CString::new(path).unwrap();
        let reader = unsafe { mdc_reader_open(c_path.as_ptr()) };

        if reader.is_null() {
            return Err(format!("Failed to open mdc file: {}", path));
        }

        let err_ptr = unsafe { mdc_reader_get_error(reader) };
        if !err_ptr.is_null() {
            let err = unsafe { CStr::from_ptr(err_ptr) };
            unsafe { mdc_reader_close(reader) };
            return Err(err.to_string_lossy().into_owned());
        }

        Ok(Self { inner: reader })
    }

    pub fn get_no_frames(&self) -> u32 {
        unsafe { mdc_get_no_frames(self.inner) }
    }

    pub fn get_segments(&self) -> Result<Vec<SegmentDesc>, String> {
        let count = unsafe { mdc_get_no_segments(self.inner) };

        let err_ptr = unsafe { mdc_reader_get_error(self.inner) };

        if !err_ptr.is_null() {
            let err = unsafe { CStr::from_ptr(err_ptr) };
            return Err(err.to_string_lossy().into_owned());
        }

        let mut segments = Vec::with_capacity(count as usize);
        for i in 0..count {
            let desc = unsafe { mdc_get_segment_desc(self.inner, i) };
            segments.push(SegmentDesc::from(desc));
        }

        Ok(segments)
    }

    pub fn get_anchor_ids(&self) -> Result<&[u32], String> {
        let count = unsafe { mdc_get_no_anchors(self.inner) };

        let err_ptr = unsafe { mdc_reader_get_error(self.inner) };

        if !err_ptr.is_null() {
            let err = unsafe { CStr::from_ptr(err_ptr) };
            return Err(err.to_string_lossy().into_owned());
        }

        let anchor_ids = unsafe { mdc_get_anchor_ids(self.inner) };

        let err_ptr = unsafe { mdc_reader_get_error(self.inner) };

        if !err_ptr.is_null() {
            let err = unsafe { CStr::from_ptr(err_ptr) };
            return Err(err.to_string_lossy().into_owned());
        }

        Ok(unsafe { std::slice::from_raw_parts(anchor_ids, count as usize) })
    }

    pub fn get_query_engine(
        &self,
        segment_ids: &[String],
        atom_ids: &[u32]
    ) -> Result<QueryEngine, String> {
        let c_strings: Vec<CString> = segment_ids
        .iter()
        .map(|s| CString::new(s.as_str()).unwrap())
        .collect();

        let c_string_ptrs: Vec<*const c_char> = c_strings
        .iter()
        .map(|s| s.as_ptr() as *const c_char)
        .collect();

        let atom_ptr = atom_ids.as_ptr();

        let engine_ptr = unsafe {
            mdc_get_query_engine(self.inner,
            c_string_ptrs.as_ptr(),
            segment_ids.len(),
            atom_ptr as *mut u32,
            atom_ids.len())
        };

        if engine_ptr.is_null() {
            let err_ptr = unsafe { mdc_reader_get_error(self.inner) };
            if !err_ptr.is_null() {
                let err = unsafe { CStr::from_ptr(err_ptr) };
                return Err(err.to_string_lossy().into_owned());
            }
            return Err("Failed to create query engine".to_string());
        }

        QueryEngine::new(engine_ptr)
    }
}

impl Drop for Reader {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe { mdc_reader_close(self.inner) };
            self.inner = ptr::null_mut();
        }
    }
}
