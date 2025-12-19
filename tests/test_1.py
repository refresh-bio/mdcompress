#!/usr/bin/env python3
import hashlib
import os
import sys
import subprocess
from pathlib import Path

# Expected md5 of example.mdc
EXPECTED_MD5 = "50ee763db617fb29b4c790758bf9e6bf"

# Directory where this script lives
HERE = Path(__file__).resolve().parent

# Paths (adjust if your layout is different)
if os.name == "nt":
    # Windows: ..\x64\release\mdcompress.exe
    mdcompress = HERE / ".." / "x64" / "release" / "mdcompress.exe"
else:
    # Linux/macOS: ../bin/mdcompress
    mdcompress = HERE / ".." / "bin" / "mdcompress"

data_dir = HERE / ".." / "examples" / "data"
input_xtc = data_dir / "example.xtc"
desc_file = data_dir / "example.desc"
tpr_file = data_dir / "example.tpr"

compressed_file = HERE / "example.mdc"
decompressed_xtc = HERE / "example.xtc"


def run_mdcompress(args):
    """Run mdcompress with given arguments, abort on error."""
    cmd = [str(mdcompress)] + [str(a) for a in args]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def md5sum(path, chunk_size=8192):
    """Compute md5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def files_equal(a, b, chunk_size=8192):
    """Compare two files byte-by-byte."""
    if os.path.getsize(a) != os.path.getsize(b):
        return False
    with open(a, "rb") as fa, open(b, "rb") as fb:
        while True:
            ca = fa.read(chunk_size)
            cb = fb.read(chunk_size)
            if not ca and not cb:
                return True
            if ca != cb:
                return False


def main():
    # Check that mdcompress exists
    if not mdcompress.is_file():
        print(f"Error: {mdcompress} does not exist. Please compile the mdcompress binary first.")
        sys.exit(1)

    # 1) compress using description file
    run_mdcompress([
        "compress",
        "-i", input_xtc,
        "-d", desc_file,
        "-o", compressed_file,
    ])

    # 2) md5 check
    md5_out = md5sum(compressed_file)
    if md5_out != EXPECTED_MD5:
        print(f"Error: md5sum does not match expected ({EXPECTED_MD5}) value. Got {md5_out}.")
        print("This may happen because something in the compressor was updated and I forgot to update this script.")
        print("In such a case I should update the expected md5sum in this script to the new value.")
        print("This check is to assure that the compressor is deterministic and gives the same output on all platforms.")
        print("It is used as a github action test.")
        sys.exit(1)

    # 3) decompress back to xtc
    run_mdcompress([
        "decompress",
        "-i", compressed_file,
        "-o", decompressed_xtc,
    ])

    if not files_equal(input_xtc, decompressed_xtc):
        print("Error: Decompressed file does not match original.")
        sys.exit(1)

    # 4) compress using topology (tpr) instead of description
    run_mdcompress([
        "compress",
        "-i", input_xtc,
        "--topology", tpr_file,
        "-o", compressed_file,
    ])

    # 5) decompress back to xtc again
    run_mdcompress([
        "decompress",
        "-i", compressed_file,
        "-o", decompressed_xtc,
    ])

    if not files_equal(input_xtc, decompressed_xtc):
        print("Error: Decompressed file does not match original.")
        sys.exit(1)

    print("All tests passed.")


if __name__ == "__main__":
    main()
