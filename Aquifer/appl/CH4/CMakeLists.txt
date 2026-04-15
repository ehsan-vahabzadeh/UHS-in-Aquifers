add_input_file_links(FILES params.input)
find_package(OpenMP REQUIRED)

dumux_add_test(NAME appl_2pnc_tpfa_CH4
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=mixingtpfa ENABLECACHING=0 SCHEIDEGGER=0
              CMD_ARGS params.input)

dumux_add_test(NAME appl_2pnc_box_CH4
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=mixingbox ENABLECACHING=0 SCHEIDEGGER=1
              CMD_ARGS params.input)

# Link UGGrid library (required for CakeGridManager with UGGrid<3>)
if(HAVE_DUNE_UGGRID)
  target_link_libraries(appl_2pnc_tpfa_CH4 PRIVATE Dune::UGGrid)
  target_link_libraries(appl_2pnc_box_CH4 PRIVATE Dune::UGGrid)
endif()
