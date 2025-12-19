rem %1 - $(SolutionDir)
rem %2 - $(Configuration)

cd %1\libs

@echo "chemfiles"
cd chemfiles
cmake -DCHEMFILES_BUILD_MMTF=ON -B build_vs
cmake --build build_vs --config %2
