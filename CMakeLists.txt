cmake_minimum_required(VERSION 3.20)
project(Live2DMeshAlgorithms LANGUAGES CXX)

# The pure geometry / texture-analysis library has no Qt or Cubism dependency.
# MeshAutoGen.cpp is the full pipeline reference; its host-document adapters
# are deliberately not included in this algorithm-only repository.
add_library(live2d_mesh_algorithms STATIC
    src/core/MeshTextureFeatures.cpp
    src/edit/DelaunayTri.cpp
    src/edit/MeshTopology.cpp)
target_compile_features(live2d_mesh_algorithms PUBLIC cxx_std_20)
target_include_directories(live2d_mesh_algorithms PUBLIC src)
if(MSVC)
    target_compile_options(live2d_mesh_algorithms PRIVATE /utf-8)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(delaunay_test tools/delaunay_test.cpp)
    target_link_libraries(delaunay_test PRIVATE live2d_mesh_algorithms)
    if(MSVC)
        target_compile_options(delaunay_test PRIVATE /utf-8)
    endif()
    add_test(NAME constrained_triangulation COMMAND delaunay_test)
endif()
