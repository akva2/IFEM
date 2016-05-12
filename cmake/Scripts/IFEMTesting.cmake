# Get GTest tests as CMake tests.
# Copied from FindGTest.cmake
# Thanks to Daniel Blezek <blezek@gmail.com> for the GTEST_ADD_TESTS code
function(gtest_add_tests executable working_dir)
    if(NOT ARGN)
        message(FATAL_ERROR "Missing ARGN: Read the documentation for GTEST_ADD_TESTS")
    endif()
    if(NOT UNIT_TEST_NUMBER)
      set(UNIT_TEST_NUMBER 0 CACHE INT "" FORCE)
    endif()
    foreach(source ${ARGN})
        file(READ "${source}" contents)
        string(REGEX MATCHALL "TEST_?[F]?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
        foreach(hit ${found_tests})
            string(REGEX REPLACE ".*\\( *([A-Za-z_0-9]+), *([A-Za-z_0-9]+) *\\).*" "\\1.\\2" test_name ${hit})
            math(EXPR UNIT_TEST_NUMBER "${UNIT_TEST_NUMBER}+1")
            set(UNIT_TEST${UNIT_TEST_NUMBER} ${test_name} ${working_dir} ${executable} --gtest_filter=${test_name} CACHE STRING "" FORCE)
        endforeach()
        # Groups parametrized tests under a single ctest entry
        string(REGEX MATCHALL "INSTANTIATE_TEST_CASE_P\\(([^,]+), *([^,]+)" found_tests2 ${contents})
        foreach(hit ${found_tests2})
          string(SUBSTRING ${hit} 24 -1 test_name)
          string(REPLACE "," ";" test_name "${test_name}")
          list(GET test_name 0 filter_name)
          list(GET test_name 1 test_prefix)
          string(STRIP ${test_prefix} test_prefix)
          math(EXPR UNIT_TEST_NUMBER "${UNIT_TEST_NUMBER}+1")
          set(UNIT_TEST${UNIT_TEST_NUMBER} ${test_prefix}.${filter_name} ${working_dir} ${executable} --gtest_filter=${filter_name}* CACHE STRING "" FORCE)
        endforeach()
    endforeach()
    set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
endfunction()


macro(IFEM_add_test_app path workdir name)
  if("${path}" MATCHES "\\*")
    file(GLOB TEST_SRCS ${path})
  else()
    set(TEST_SRCS ${path})
  endif()
  add_executable(${name}-test EXCLUDE_FROM_ALL ${IFEM_PATH}/src/IFEM-test.C ${TEST_SRCS})
  gtest_add_tests($<TARGET_FILE:${name}-test> ${workdir} ${TEST_SRCS})
  list(APPEND TEST_APPS ${name}-test)
  target_link_libraries(${name}-test ${ARGN} gtest pthread)
endmacro()

macro(IFEM_add_unittests IFEM_PATH)
  IFEM_add_test_app("${IFEM_PATH}/src/Utility/Test/*.C;${IFEM_PATH}/src/ASM/Test/*.C;${IFEM_PATH}/src/LinAlg/Test/*.C;${IFEM_PATH}/src/SIM/Test/*.C"
                    ${IFEM_PATH}
                    IFEM
                    ${IFEM_LIBRARIES} ${IFEM_DEPLIBS})

  # Parallel unit tests. These are grouped under one CTest
  if(MPI_FOUND)
    set(TEST_SRCS ${IFEM_PATH}/src/ASM/Test/MPI/TestDomainDecomposition.C)
    add_executable(IFEM-MPI-test EXCLUDE_FROM_ALL ${IFEM_PATH}/src/IFEM-test.C ${TEST_SRCS})
    target_link_libraries(IFEM-MPI-test ${IFEM_LIBRARIES} ${IFEM_DEPLIBS} gtest)
    add_test(NAME MPI-tests WORKING_DIRECTORY ${IFEM_PATH} COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:IFEM-MPI-test>)
    list(APPEND TEST_APPS IFEM-MPI-test)
  endif()
endmacro()

function(IFEM_add_test name binary)
  separate_arguments(MEMCHECK_COMMAND)
  if(IFEM_TEST_EXTRA)
    set(test-name "${binary}+${IFEM_TEST_EXTRA}+${name}")
  else()
    set(test-name "${binary}+${name}")
  endif()
  if(IFEM_TEST_MEMCHECK)
    add_test("${test-name}" regtest.sh "${MEMORYCHECK_COMMAND} --log-file=${CMAKE_BINARY_DIR}/valgrindlog ${EXECUTABLE_OUTPUT_PATH}/${binary}" ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${name} ${ARGN})
  else(IFEM_TEST_MEMCHECK)
    add_test("${test-name}" regtest.sh ${EXECUTABLE_OUTPUT_PATH}/${binary} ${PROJECT_SOURCE_DIR}/${TEST_SUBDIR}/Test/${name} ${ARGN})
  endif(IFEM_TEST_MEMCHECK)
endfunction()

macro(add_check_target)
  add_custom_target(check ${CMAKE_CTEST_COMMAND} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  add_custom_command(TARGET check PRE_BUILD COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/failed.log)
  if(IFEM_AS_SUBMODULE OR IFEM_LIBRARY_BUILD)
    ifem_add_unittests(${IFEM_PATH})
  endif()
  if(NOT TARGET gtest)
    add_subdirectory(${IFEM_PATH}/3rdparty/gtest gtest EXCLUDE_FROM_ALL)
  endif()
  if (${UNIT_TEST_NUMBER} GREATER 0)
    foreach(test_number RANGE 1 ${UNIT_TEST_NUMBER})
      list(GET UNIT_TEST${test_number} 0 name)
      list(GET UNIT_TEST${test_number} 1 dir)
      list(GET UNIT_TEST${test_number} 2 -1 cmd)
      add_test(NAME ${name} WORKING_DIRECTORY ${dir} COMMAND ${cmd})
    endforeach()
  endif()
  add_dependencies(check ${TEST_APPS})
  add_custom_target(testapps DEPENDS ${TEST_APPS})

  if(NOT TARGET check-commits)
    add_custom_target(check-commits
                      COMMAND ${CMAKE_COMMAND}
                              -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                              -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
                              -P ${IFEM_CHECKCOMMITS_SCRIPT})
  endif()
endmacro()

set(IFEM_TESTING_INCLUDED 1)
if(IFEM_INTREE_BUILD)
  include_directories(${IFEM_PATH}/3rdparty/gtest/include)
elseif(NOT IFEM_AS_SUBMODULE AND NOT IFEM_LIBRARY_BUILD
       AND NOT TARGET gtest)
  if (EXISTS ${PROJECT_SOURCE_DIR}/../gtest)
    add_subdirectory(${PROJECT_SOURCE_DIR}/../gtest gtest)
    include_directories(${PROJECT_SOURCE_DIR}/../gtest/include)
  else()
    include(DownloadGTest)
    include_directories(${GTEST_INCLUDE_DIRS})
  endif()
endif()
