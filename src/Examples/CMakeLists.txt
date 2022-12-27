project(raw_pdb-Examples)

set(SOURCES
	ExampleContributions.cpp
	ExampleFunctionSymbols.cpp
	ExampleLines.cpp
	ExampleMain.cpp
	ExampleMemoryMappedFile.cpp
	ExampleMemoryMappedFile.h
	Examples_PCH.cpp
	Examples_PCH.h
	ExampleSymbols.cpp
	ExampleTypes.cpp
	ExampleTimedScope.cpp
	ExampleTimedScope.h
)

source_group(src FILES
    ${SOURCES}
)

add_executable(raw_pdb-Examples
    ${SOURCES}
)

target_link_libraries(raw_pdb-Examples
  PUBLIC
    raw_pdb
)

target_precompile_headers(raw_pdb-Examples
  PUBLIC
    Examples_PCH.h
)