get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})

# Find lp-engine
find_package(lp-engine REQUIRED)
include_directories(${lp-engine_INCLUDE_DIRS})

add_executable( ${ProjectId} main.cpp)
target_link_libraries(${ProjectId}
                      lp-engine
                      )
