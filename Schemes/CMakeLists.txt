# TestCase
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase2d.cpp)

# HHO
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> HHO-diffusion/DiffusionEquation2d.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh2d quadrature hybridcore2d plot ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE})
TARGET_INCLUDE_DIRECTORIES(hho-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})
