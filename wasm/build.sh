#!/bin/bash

set -e  # Exit on error

# Check if em++ is available
if command -v em++ &>/dev/null; then
    echo "em++ is already installed."
else
    echo "em++ not found. Installing Emscripten..."

    # Check if emsdk already exists
    if [ ! -d "emsdk" ]; then
        echo "Cloning Emscripten SDK..."
        git clone https://github.com/emscripten-core/emsdk.git
    fi
    cd emsdk
    ./emsdk install latest
    ./emsdk activate latest
    source ./emsdk_env.sh
    cd ..
fi

#emsdk is needed to compile, instructions here https://emscripten.org/docs/getting_started/downloads.html
# but also here :)

# git clone https://github.com/emscripten-core/emsdk.git
# cd emsdk
# ./emsdk install latest
# ./emsdk activate latest
# source ./emsdk_env.sh

# Verify em++ is available now
if ! command -v em++ &>/dev/null; then
    echo "Error: em++ is still not found after sourcing. Something went wrong."
    exit 1
fi


mkdir -p bin
mkdir -p build

em++ -c -O3 -std=c++20 -I ../src ../src/common/engine_decompression.cpp -o build/engine_decompression.cpp.o
em++ -c -O3 -std=c++20 -I ../src ../src/common/serializer.cpp -o build/serializer.cpp.o
em++ -c -O3 -std=c++20 -I ../src ../src/common/mdc_decompress.cpp -o build/mdc_decompress.cpp.o
#mkokot_TODO: there is some problem with static_assert in refres::archive fseek, but its not used here anyway so this ugly workaround is used
em++ -c -O3 -std=c++20 -DSKIP_FSEEK_STATIC_ASSERT -I ../src -I ../libs ../src/mdc_lib/mdc_reader.cpp -o build/mdc_reader.cpp.o

em++ -c -O3 -std=c++20 -I ../src/mdc_lib src/wasm_bind.cpp -o build/wasm_bind.cpp.o
emar rcs bin/libmdc.a build/*.o

emcc -O3 -lembind -o bin/libmdc.js \
    -Wl,--whole-archive bin/libmdc.a -Wl,--no-whole-archive \
    -s ALLOW_MEMORY_GROWTH=1 \
	-s FORCE_FILESYSTEM=1 \
    -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap", "UTF8ToString", "FS"]'

mkdir -p test_project
cp bin/libmdc.js test_project
cp bin/libmdc.wasm test_project

# to run local server
# (cd test_project && python3 -m http.server 8080)
