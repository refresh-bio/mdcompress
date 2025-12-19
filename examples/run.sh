#!/bin/bash

if [ ! -f ../bin/mdcompress ]; then
	echo "Error: ../bin/mdcompress does not exist. Please compile the mdcompress binary first."
	exit 1
fi
# simple command to just compress xtc using description file
../bin/mdcompress compress -i data/example.xtc -d data/example.desc -o example.mdc

md5sum_out=$(md5sum example.mdc | awk '{print $1}')
md5sum_expected="50ee763db617fb29b4c790758bf9e6bf"
if [ "$md5sum_out" != "$md5sum_expected" ]; then
	echo "Error: md5sum does not match expected ($md5sum_expected) value. Got $md5sum_out". 
	echo "This may happen because something in the compressor was updated and I forgot to update this script."
	echo "In such a case I should update the expected md5sum in this script to the new value."
	echo "This check is to assure that the compressor is deterministic and gives the same output on all platforms."
	echo "It is used as a github action test."
	exit 1
fi
# now lets decompress back to xtc
../bin/mdcompress decompress -i example.mdc -o example.xtc

diff data/example.xtc example.xtc
if [ $? -ne 0 ]; then
	echo "Error: Decompressed file does not match original."
	exit 1
fi

# very often there is also a tpr file, so mdcompress can use it instead desc file, such that user does not need to prepare desc file
../bin/mdcompress compress -i data/example.xtc --topology data/example.tpr -o example.mdc
# now lets decompress back to xtc
../bin/mdcompress decompress -i example.mdc -o example.xtc

diff data/example.xtc example.xtc
if [ $? -ne 0 ]; then
	echo "Error: Decompressed file does not match original."
	exit 1
fi
