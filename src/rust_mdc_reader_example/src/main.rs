use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

use rust_mdc::{Reader, QueryResult};

struct Params {
    input_file: Option<String>,
    output_file: Option<String>,
    segment_ids: Vec<String>,
    atom_ids: Vec<u32>,
}

fn parse_segment_ids(s: &str) -> Vec<String> {
    s.split(',').map(|x| x.to_string()).collect()
}

fn parse_atom_ids(s: &str) -> Vec<u32> {
    let mut atom_ids = vec![];
    for part in s.split(',') {
        if part.contains('-') {
            let range: Vec<&str> = part.split('-').collect();
            if range.len() != 2 {
                eprintln!("Invalid range: {}", part);
                continue;
            }
            let start = range[0].parse::<u32>();
            let end = range[1].parse::<u32>();
            if start.is_err() || end.is_err() {
                panic!("Invalid range: {}", part);
            }
            let start = start.unwrap();
            let end = end.unwrap();
            if start > end {
                panic!("Invalid range: {}", part);
            }
            for i in start..=end {
                atom_ids.push(i);
            }
        } else {
            let id = part.parse::<u32>();
            if id.is_err() {
                panic!("Invalid atom id: {}", part);
            }
            atom_ids.push(id.unwrap());
        }
    }
    atom_ids
}

fn parse_params() -> Result<Params, String> {
    let mut params = Params {
        input_file: None,
        output_file: None,
        segment_ids: vec![],
        atom_ids: vec![],
    };
    
    let args: Vec<String> = env::args().collect();
    let mut i = 1; 
    while i < args.len() {
        match args[i].as_str() {
            "-i" => {
                i += 1;
                if i >= args.len() {
                    return Err("Error: -i requires an argument".to_string());
                }
                params.input_file = Some(args[i].clone());
            }
            "-o" => {
                i += 1;
                if i >= args.len() {
                    return Err("Error: -o requires an argument".to_string());
                }
                params.output_file = Some(args[i].clone());
            }
            "--segments" => {
                i += 1;
                if i >= args.len() {
                    return Err("Error: --segments requires an argument".to_string());
                }
                params.segment_ids = parse_segment_ids(args[i].as_str());
            }
            "--atoms" => {
                i += 1;
                if i >= args.len() {
                    return Err("Error: --atoms requires an argument".to_string());
                }
                params.atom_ids = parse_atom_ids(args[i].as_str());
            }
            _ => {
                return Err(std::format!("Unknown argument: {}", args[i]));
            }
        }
        i += 1;
    }

    if params.input_file.is_none() || params.output_file.is_none() {
        return Err(std::format!("Usage: {} -i <input> -o <output> [--segments seg1,seg2,...] [--atoms id1,start_range1-end_range1,...]", args[0]));
    }
    Ok(params)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    
    let params = parse_params().map_err(|e| {
        format!("Failed to parse params: {}", e)
    })?;

    let reader = Reader::new(&params.input_file.unwrap()).map_err(|e| {
        format!("Failed to open reader: {}", e)
    })?;
    
    let segments = reader.get_segments().map_err(|e| {
        format!("Failed to get segments: {}", e)
    })?;

    let no_frames = reader.get_no_frames();

    println!("Segments: {}, Frames: {}", segments.len(), no_frames);

    for segment in segments {
        println!("Segment: {}, size: {}, type: {}", segment.name, segment.size, segment.type_);
    }

    let engine = reader.get_query_engine(&params.segment_ids, &params.atom_ids).map_err(|e| {
        format!("Failed to get query engine: {}", e)
    })?;


    let file = File::create(&params.output_file.unwrap())?;
    let mut writer = BufWriter::new(file);
    
    let mut result = QueryResult::new().map_err(|e| {
        format!("Failed to create query result: {}", e)
    })?;

    let mut frame_ids: Vec<u32> = vec![0];
    
    for frame_id in 0..no_frames {
        frame_ids[0] = frame_id;
       
        engine.query(&frame_ids, &mut result).map_err(|e| {
            format!("Failed to query: {}", e)
        })?;

        let frame = &result.frames[0];
        let total_atoms = frame.coords.len();
        
        let coords = frame.coords;
        // Each atom has 3 coords
        // Just write them raw:
        if !coords.is_empty() {
            let bytes = unsafe {
                std::slice::from_raw_parts(
                    coords.as_ptr() as *const u8,
                    coords.len() * std::mem::size_of::<f32>() * 3
                )
            };
            writer.write_all(bytes)?;
        }

        println!("Frame time: {}, No. atoms: {}", frame.time, total_atoms);
    }

    writer.flush()?;
    Ok(())
}
