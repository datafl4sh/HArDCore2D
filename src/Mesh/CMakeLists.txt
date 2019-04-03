# Library: mesh2d

FILE(GLOB source_list mesh2d.* mesh2d_builder.* import_mesh2d.* cell2d.* edge2d.* vertex2d.*)
ADD_LIBRARY(mesh2d ${source_list})

