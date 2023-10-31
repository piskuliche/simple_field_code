add_library(field_module "")

target_sources(field_module
    PRIVATE
        HDF5::HDF5
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/field_module.f90
  )