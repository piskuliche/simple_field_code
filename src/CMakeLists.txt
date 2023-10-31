add_executable(get_field_calc fieldcode.f90)

add_subdirectory(field_module)
target_link_libraries(field_module
    PRIVATE
     HDF5
)
target_link_libraries(get_field_calc 
    PRIVATE
        field_module
)