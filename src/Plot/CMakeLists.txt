# Library: plot

FILE(GLOB source_list vtu_writer2d.*)
ADD_LIBRARY(plot ${source_list})
TARGET_LINK_LIBRARIES(plot mesh2d)
