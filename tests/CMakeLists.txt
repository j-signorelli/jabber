add_executable(unit-tests
                        test_compute_1d.cpp
                        test_compute_2d.cpp
                        test_compute_3d.cpp
                        test_app_common.cpp
                        test_app_config_input.cpp
                        test_interpolant.cpp
                        test_psd.cpp
                        test_wave.cpp
                        test_transfer_functions.cpp)

set(dependencies Catch2::Catch2WithMain jabber)

target_link_libraries(unit-tests PRIVATE ${dependencies})


include(CTest)
include(Catch)
catch_discover_tests(unit-tests)