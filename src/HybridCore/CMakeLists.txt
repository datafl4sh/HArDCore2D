# Library: hybridcore2d

FILE(GLOB source_list hybridcore2d.*)
ADD_LIBRARY(hybridcore2d ${source_list})
TARGET_LINK_LIBRARIES(hybridcore2d mesh2d quadrature)


