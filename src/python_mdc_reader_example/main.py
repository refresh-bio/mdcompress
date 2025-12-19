#!/bin/env python3
import mdc
import sys
import numpy as np
import argparse

def parse_segments(desc: str):
    if len(desc) == 0:
        return []
    return desc.split(",")

def parse_atoms(desc: str):
    atoms = []

    for atom_desc in desc.split(","):
        if len(atom_desc) == 0:
            continue
        if "-" in atom_desc:
            start, end = atom_desc.split("-")
            atoms.extend(range(int(start), int(end)+1))
        else:
            atoms.append(int(atom_desc))
    return atoms

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Reads an MDC file and stores coords as binary file with raw floats")
    argparser.add_argument("-i", "--input", type=str, help="Input MDC file", required=True)
    argparser.add_argument("-o", "--output", type=str, help="Output file", required=True)
    argparser.add_argument("--segments", type=str, help="comma separated list of segments to read", default="")
    argparser.add_argument("--atoms", type=str, help="comma separated list of atoms to read or ranges (e.g. 1-10,50,200-300)", default="")
    args = argparser.parse_args()

    print("----------------  Configuration  ---------------")
    for arg, value in vars(args).items():
        print("%s: %s" % (arg, value))
    print("------------------------------------------------", flush=True)
    
    segments = parse_segments(args.segments)
    atoms = parse_atoms(args.atoms)

    reader = mdc.Reader(args.input)
    error = reader.get_current_error()
    if error:
        print(f"Error: {error}")
        sys.exit(1)

    query_engine = reader.get_query_engine(segments, atoms)
    error = query_engine.get_current_error()
    if error:
        print(f"Error: {error}")
        sys.exit(1)
    
    no_segments = len(reader.get_segments())
    no_frames = reader.get_no_frames()

    print(f"Segments: {no_segments}, Frames: {no_frames}")

    for segment in reader.get_segments():
        print(f"Segment: {segment.name}, size: {segment.size}, type: {segment.type}")

    result = mdc.QueryResult()
    frame_ids = np.array([0])

    # Open output file to write binary data
    with open(args.output, "wb") as output_file:
        for frame_id in range(no_frames):
            frame_ids[0] = frame_id
            if not query_engine.query(frame_ids, result):
                print("Error: {}", query_engine.get_current_error())
                sys.exit(1)

            assert(len(result.frames) == 1)
            print(f"Frame time: {result.frames[0].time}\tNo. atoms: {result.frames[0].coords.shape[0]}")
            output_file.write(result.frames[0].coords.tobytes())

