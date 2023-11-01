add_executable(get_field_calc fieldcode.f90)

add_subdirectory(field_module)

target_link_libraries(field_module HDF5::HDF5)

target_link_libraries(get_field_calc
        field_module
)