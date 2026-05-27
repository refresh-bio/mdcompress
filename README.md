MDCompress
=
MDCompress is a compressor of Molecular Dynamics files stored in the XTC and other formats (XTC, TRR, DCD, ...). MDCompress uses [chemfiles](https://github.com/chemfiles/chemfiles) for reading during compression and for writing during decompression, so each format that is supported by chemfiles is also supported by MDCompress.

## Quick start
To compile mdcompress g++ compiler (at least g++-12 is required) for both Linux and mac os (for mac os we advise to use g++ installed with [homebrew](https://brew.sh/).
```
git clone --recurse-submodules https://github.com/refresh-bio/mdcompress
cd mdcompress

# Linux compilation
make -j

# MacOS compilation (must specify g++ compiler)
gmake mdcompress CXX=g++-13 CC=gcc-13 -j

# Example: compression of XTC file using TPR file
bin/mdcompress compress -i examples/data/example.xtc --topology examples/data/example.tpr -o example.mdc

# Alternative using a simple text file describing topology (-d)
# bin/mdcompress compress -i examples/data/example.xtc -d examples/data/example.desc -o example.mdc

# Show info 
bin/mdcompress info -i example.mdc

# Decompress MDC file
bin/mdcompress decompress -i example.mdc -o aaa.xtc

# Select frame 1000
bin/mdcompress select -i example.mdc -o aaa_1000.xtc --fid 1000

# Select frames from 1000 to 1500
bin/mdcompress select -i example.mdc -o aaa_1000_1500.xtc --fr 1000 1500

# Select frames from 1000 to 1500 with a stride of 15
bin/mdcompress select -i example.mdc -o aaa_1000_1500_15.xtc --fr 1000 1500 --stride 15

# Select only MOL parts (no water and 'other')
bin/mdcompress select -i example.mdc -o aaa_MOL.xtc --segments protein/DNA/RNA,small_molecule

# Select only MOL and OTHER parts (no water)
bin/mdcompress select -i prod.part0001w.mdc -o aaa_MOL_OTHER.xtc --segments protein/DNA/RNA,small_molecule,ion

# INFO: names "protein/DNA/RNA", "small_molecule", "ion" are defined in the description file defined during compression with -d parameter (examples/data/example.desc)
```

## Dependencies
For the above to work correctly, some packages must be installed.
For Debina/Ubuntu they may be installed with the following command:
```
sudo apt install make cmake build-essential python3 python3-dev rustc cargo git
```

## Usage
`mdcompress <mode> [options]`

Modes:
* `compress`         &ndash; compress a trajectory (XTC/TRR/DCD/...) into MDC format
* `decompress`       &ndash; decompress MDC file into a trajectory (XTC/TRR/DCD/...)
* `select`           &ndash; decompress some frames from MDC file into a trajectory (XTC/TRR/DCD/...)
* `info`             &ndash; info about contents of MDC file
* `make_desc`        &ndash; create description file (for -d) from a topology (TPR/PSF/...)

Options - compress mode:
* `-i <file_name>`               &ndash; input trajectory (XTC, TRR, DCD, ...)
* `-o <file_name>`               &ndash; output file name (`.mdc`)
* `-d <file_name>`               &ndash; description of the segments of a frame. May be omitted if `--topology` is given, in which case mdcompress builds the description from the topology automatically.
* `--topology <file_name>`       &ndash; topology file (TPR, PSF, PDB, ...). It is stored inside the `.mdc`. If `-d` is not given, mdcompress infers the description from it. Note: coordinate-only trajectories (e.g. DCD) carry no topology, so either `-d` or `--topology` must be provided for them. It may happen that mdcompress cannot build a description from the file specified with the `--topology` switch; in that case, `-d` must be provided.
* `--only-mol`                   &ndash; compress only the molecule segments, skipping water/ions/other from the description/topology. Usually not needed: if the trajectory contains only molecules, mdcompress detects this and enables `--only-mol` automatically (with a warning).
* _description check_            &ndash; whether `-d` or `--topology` is used, the description is validated against the trajectory's atom count before compressing: an exact match is used as-is; if only the molecules match, `--only-mol` is enabled automatically; if nothing matches, compression is aborted with an explanation and (when the description was inferred from a topology) a candidate is written to `<output>.candidate.desc` for you to edit and reuse with `-d`.
* `-l <int>`                     &ndash; compression level
* `--preset <preset>`            &ndash; use one of preset parameters for `-b`, `--subsegment-size`,
available presets:
  * `default`,
  * `archive` - use if the main decompression use scheme is decompress all,
  * `trajectory` - use if the main decompression use scheme is tracking trajectory of some atoms,
  * `frames` - use if the main decompression use scheme is selecting subset of frames.
* `-t <int>`                     &ndash; no. threads
* `-b|--batch-size <int>`        &ndash; no. of frames in a batch (default: 20)
* `-h|--max-history-size <int>`  &ndash; no. of previous framed used to predict the current one (default: 1; max: 3)
* `--n-frames-for-model`         &ndash; no. of frames to build model (default: 50)
* `--res <int>`                  &ndash; min. resolution in fm (default: 1000; 1000 fm = 0.01 Angstrom)
* `--max-dist-in-model`          &ndash; max. distance in segment of reference atoms (default: 100)
* `--subsegment <id1,id2,...>`   &ndash; list of segments to split into subsegments (if no specified subsegment all except water)
* `--subsegment-size <int>`      &ndash; number of atoms in a single subsegment (default: 100, 0 means don't use subsegments)

Options - decompress mode
* `-i <file_name>`               &ndash; input file name (`.mdc`)
* `-o <file_name>`               &ndash; output trajectory file name (format chosen by extension: `.xtc`, `.trr`, `.dcd`, ...)
* `--topology <file_name>`       &ndash; if a topology was stored at compression time, write it out to this path (the extension must match the original topology file)

Options - select mode
* `-i <file_name>`               &ndash; input file name (`.mdc`)
* `-o <file_name>`               &ndash; output trajectory file name (format chosen by extension: `.xtc`, `.trr`, `.dcd`, ...)
* `--topology <file_name>`       &ndash; if a topology was stored at compression time, write it out to this path (the extension must match the original topology file)
* `--fid <int>`                  &ndash; frame id (0-based)
* `--fr <int> <int|MAX>`         &ndash; range of frame ids (0-based) (`MAX` or 2147483647 means last frame) 
* `--stride <int>`               &ndash; stride size (for range of frames) (default: 1)
* `--segments <id1,id2,...>`     &ndash; list of segment ids
* `--atoms <id1,id2,...>`        &ndash; list of atoms to tracks, `id<n>` may be single id or closed interval in format start-end

Options &ndash; info mode
* `-i <file_name>`               &ndash; input file name

Options &ndash; make_desc mode
* `-i <file_name>`               &ndash; input topology file (TPR, PSF, ...)
* `-o <file_name>`               &ndash; output desc file name
* `--only-mol`                   &ndash; include only molecule (skip water and 'other')

## How to get the description of the segments of the frame
MDCompress needs to distinguish the segments of a frame (molecule, water, other).
The most direct way to pass this info is to use `-d` parameter.
This, however, is not very convenient, and MDCompress is able (in most cases) to infer the required information from the topology file (like TPR, PSF).
In this case `--topology` parameter pointing to such a file should be used instead of `-d`.

In other words, you can hand MDCompress both the trajectory and the topology in a single call, and it will build the description, check that it makes sense, and use it:
```
bin/mdcompress compress -i trajectory.xtc --topology topology.tpr -o out.mdc
```

**Sanity check.** Whether the description comes from `-d` or is inferred from `--topology`, MDCompress validates it against the actual trajectory before compressing:
* If the total number of atoms in the description matches the trajectory, it is used as-is.
* If only the *molecule* atoms match the trajectory (a common case where the topology also describes water/ions that are not present in the trajectory), MDCompress prints a prominent warning and automatically enables `--only-mol`, compressing only the molecule segments.
* If neither matches, MDCompress stops, prints the inferred segments and atom counts, and writes a best-effort candidate description to `<output>.candidate.desc` so you can edit it (segment types `MOL`/`WAT`/`OTHER`/`NONE`) and re-run with `-d`.

It happens that the topology file describes the full topology, but the XTC/TRR file contains only molecules; in such a case `--only-mol` flag may be added explicitly (or left to the auto-detection described above).

**Important note**: `--only-mol` is designed to inform MDCompress that XTC/TRR file contains only molecules, so it knows to skip other segments from the topology file.
**It is not intended to make MDCompress skip non-molecules from XTC/TRR.** If this is needed, the description file with type `NONE` should be used (read below).

Flag `--only-mol` may not be applied if XTC/TRR contains all the atoms described in the topology (or description) file, but we want to have only molecule atoms in the compressed file.
If `--topology` is provided, it will also be stored in the compressed archive.

If there is no topology file (TPR/PSF), the only option is to use `-d` parameter to describe frame segments.
In such a case, `--topology` may also be provided, but it will be only stored in the compressed archive, and not used to infer the segments of the frame.
Description file (`-d`) is a text file. Each component of the frame is described in one line, and the order of lines should reflect the order of segments in the frame.
Each line is in the format:
```
<segment type> <segment name> <number of atoms in segment>
```

An example of a description file is [here](https://github.com/refresh-bio/mdcompress-dev/blob/master/examples/data/example.desc).
`<segment type>` may be one of:
 * `MOLECULE` (or `MOL`)  &ndash; molecule segment
 * `WATER` (or `WAT`)     &ndash; water segment
 * `NONE` (or `NONE`)     &ndash; MDCompress will skip these atoms
 * `OTHER` (or `OTH`)     &ndash; use it for ions and other types 

#### Use case: I only have TPR, but I want to compress only molecules, and I don't know how to make a description file
Don't worry! MDCompress may prepare a description file for you, such that you may edit it. 
For example, you may set different names or exclude some segments by setting their type to `NONE`.
Assuming your TPR file name is `topology.tpr`, you may just:
```
bin/mdcompress make_desc -i topology.tpr -o topology.desc
```
and edit `topology.desc` as needed and pass it for compression with `-d` switch.


## A note on DCD (and other topology-less formats)
DCD trajectories store only coordinates (and the box) - they contain no topology, atom names, or bonds.
Because of this, the description cannot be inferred from a DCD file itself; you must provide segment information either with `-d` (a description file) or with `--topology` pointing to a real topology file (e.g. TPR or PSF) that describes the same atoms.
```
# Using a description file
bin/mdcompress compress -i trajectory.dcd -d topology.desc -o out.mdc

# Or letting MDCompress infer (and validate) the description from a topology file
bin/mdcompress compress -i trajectory.dcd --topology structure.psf -o out.mdc

# Decompress back to DCD (or any other supported format)
bin/mdcompress decompress -i out.mdc -o trajectory.dcd
```
DCD coordinates are stored as 32-bit floats in Angstroms. By default MDCompress uses a resolution of 0.01 A; use `--res` to change it if you need more (or less) precision.

## Tracking specific atoms
It is possible to select only specific atoms using the `--atoms` switch. 
Example:
```
mdcompress select -i example.mdc -o atoms2000-2100_and_3000.xtc --atoms 2000-2100,3000
```
