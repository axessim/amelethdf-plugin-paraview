
INCLUDE_DIRECTORIES(${AMELETHDF_INCLUDE_DIRS})

# build fortran tools
FILE(GLOB tools_src *.cxx)
FOREACH(tool_file ${tools_src})
    FILE(RELATIVE_PATH rel_tool_file ${CMAKE_CURRENT_SOURCE_DIR} ${tool_file})
    STRING(REPLACE ".cxx"  "" tool_exec ${rel_tool_file})
    ADD_EXECUTABLE(${tool_exec} ${tool_file})
    
    TARGET_LINK_LIBRARIES(${tool_exec} ${VTKAH5_LIBRARIES})
    TARGET_LINK_LIBRARIES(${tool_exec} ${VTKAH5_EXTRA_TOOLS_LIBRARIES})
    INSTALL(TARGETS ${tool_exec} DESTINATION bin)
ENDFOREACH(tool_file)