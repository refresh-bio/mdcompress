#!/bin/env python3
import mdc
import sys
import numpy as np

reader = mdc.Reader("example.mdc")

no_frames = reader.get_no_frames()

print(f"Segments: {reader.get_no_frames()}, Frames: {no_frames}")

for segment in reader.get_segments():
    print(f"Segment: {segment.name}, size: {segment.size}, type: {segment.type}")

#query_engine to get all atoms
query_engine = reader.get_query_engine([], [])

query_result = mdc.QueryResult()

# queries require np.array
frame_ids = np.array([0])

for frame_id in range(no_frames):
    frame_ids[0] = frame_id
    if not query_engine.query(frame_ids, query_result):
        print("Error: {}", query_engine.get_current_error())
        sys.exit(1)

    print(query_result.frames[0].coords)

    print(f"Frame time: {query_result.frames[0].time}\tNo. atoms: {query_result.frames[0].coords.shape[0]}")
