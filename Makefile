all: mdcompress libmdc c_mdc_reader_example cpp_mdc_reader_example rust_mdc_reader_example python_mdc_reader_example python_mdc

#I Just took it from splash Makefile...
PLATFORM?=generic

# *** REFRESH makefile utils
include refresh.mk

$(call INIT_SUBMODULES)
$(call INIT_GLOBALS)
$(call CHECK_OS_ARCH, $(PLATFORM))

# *** Project directories
$(call SET_SRC_OBJ_BIN,src,obj,bin)
$(call SET_3RD_PARTY,./libs)

# *** Project configuration
$(call ADD_MIMALLOC, $(3RD_PARTY_DIR)/mimalloc)
$(call ADD_CHEMFILES, $(3RD_PARTY_DIR)/chemfiles)
$(call ADD_REFRESH_LIB, $(3RD_PARTY_DIR))
$(call SET_STATIC, $(STATIC_LINK))
$(call SET_C_CPP_STANDARDS, c11, c++20)
$(call SET_GIT_COMMIT)

$(call SET_FLAGS, $(TYPE))

$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_x86_64, 12, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_aarch64, 12, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_x86_64, 12, 14)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_arm64, 12, 14)

ifneq ($(MAKECMDGOALS),clean)
$(call CHECK_COMPILER_VERSION)
endif

# DEFINE_FLAGS +=-DTEST_MDCOMPRESS_API

INCLUDE_DIRS+=-Isrc/mdc_lib

# *** Source files and rules

$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,COMMON,common))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,MDC_LIB,mdc_lib))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,MAIN,app))
$(eval $(call PREPARE_DEFAULT_C_COMPILE_RULE,C_MDC_READER_EXAMPLE,c_mdc_reader_example))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,CPP_MDC_READER_EXAMPLE,cpp_mdc_reader_example))

$(eval $(call PREPARE_DEFAULT_LIBC_COMPILE_RULE,LIBXDRFILE,libxdrfile))


#CPP_FLAGS+=-pg
#LINKER_FLAGS+=-pg

# *** Targets
mdcompress: $(OUT_BIN_DIR)/mdcompress

$(OUT_BIN_DIR)/mdcompress: mimalloc_obj \
	$(OBJ_MAIN) $(OBJ_COMMON) $(OBJ_LIBXDRFILE) $(OBJ_MDC_LIB)
	-mkdir -p $(OUT_BIN_DIR)	
	$(CXX) -o $@  \
	$(MIMALLOC_OBJ) \
	$(OBJ_MAIN) $(OBJ_COMMON) $(OBJ_LIBXDRFILE) $(OBJ_MDC_LIB) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

libmdc: $(OUT_BIN_DIR)/libmdc.a

$(OUT_BIN_DIR)/libmdc.a: $(OBJ_COMMON) $(OBJ_MDC_LIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(AR) $(AR_OPT) $@ $(OBJ_COMMON) $(OBJ_MDC_LIB)

c_mdc_reader_example: $(OUT_BIN_DIR)/c_mdc_reader_example

$(OUT_BIN_DIR)/c_mdc_reader_example: $(OUT_BIN_DIR)/libmdc.a $(OBJ_C_MDC_READER_EXAMPLE)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ \
	$(OBJ_C_MDC_READER_EXAMPLE) $(OUT_BIN_DIR)/libmdc.a -lstdc++ \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

cpp_mdc_reader_example: $(OUT_BIN_DIR)/cpp_mdc_reader_example

$(OUT_BIN_DIR)/cpp_mdc_reader_example: $(OUT_BIN_DIR)/libmdc.a $(OBJ_CPP_MDC_READER_EXAMPLE)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) -o $@ \
	$(OBJ_CPP_MDC_READER_EXAMPLE) $(OUT_BIN_DIR)/libmdc.a -lstdc++ \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

rust_mdc_reader_example: $(OUT_BIN_DIR)/rust_mdc_reader_example

#cargo build in the in the example will call cargo build for rust_mdc which will call this makefile to build libmcd.a (if not exist)
$(OUT_BIN_DIR)/rust_mdc_reader_example: libmdc
	-mkdir -p $(OUT_BIN_DIR)
	(cd src/rust_mdc_reader_example && cargo build --release && cp target/release/rust_mdc_reader_example ../../$(OUT_BIN_DIR)/rust_mdc_reader_example)

python_mdc: $(OUT_BIN_DIR)/python_mdc

$(OUT_BIN_DIR)/python_mdc: libmdc
	-mkdir -p $(OUT_BIN_DIR)
	(cd python && cmake $(CMAKE_OSX_FIX) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON && cmake --build build -- -j$(nproc) && cp build/mdc*.so ../bin)

python_mdc_reader_example: $(OUT_BIN_DIR)/python_mdc_reader_example.py

$(OUT_BIN_DIR)/python_mdc_reader_example.py: python_mdc
	-mkdir -p $(OUT_BIN_DIR)
	cp src/python_mdc_reader_example/main.py bin/python_mdc_reader_example.py

# *** Cleaning
.PHONY: clean init
clean: clean-mimalloc_obj
	-rm -r $(OBJ_DIR)
	-rm -r $(OUT_BIN_DIR)

init:
	$(call INIT_SUBMODULES)

