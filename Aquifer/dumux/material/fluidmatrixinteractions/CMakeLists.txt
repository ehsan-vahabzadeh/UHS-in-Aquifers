add_subdirectory(dispersiontensors)
add_subdirectory(2p)
file(GLOB DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions)

