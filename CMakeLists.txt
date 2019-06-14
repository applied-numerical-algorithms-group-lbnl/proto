cmake_minimum_required(VERSION 3.9)
project(proto)

set(CMAKE_CXX_STANDARD 11)
set(PROTO_INC ./include  ./examples/Euler/src ./out)

if (NOT EXISTS(pdfl))
    include(ExternalProject)
    ExternalProject_Add(pdfl
        GIT_REPOSITORY https://github.com/BoiseState/PolyhedralDataflowIR.git
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=pdfl)
endif()

set(PDFL_INC
    pdfl/src
    pdfl/lib/iegenlib/src
    pdfl/lib/isl
    pdfl/lib/isl/imath
    pdfl/lib/isl/include
    pdfl/lib/isl/bld
    pdfl/lib/isl/bld/include
    pdfl/lib/chill
    pdfl/lib/chill/omega/basic/include
    pdfl/lib/chill/omega/omega_lib/include
    pdfl/lib/chill/omega/code_gen/include
    pdfl/lib/chill/omega/parser
    pdfl/lib/gmp
    pdfl/test
)

set(GTEST_INC pdfl/lib/gtest pdfl/lib/gtest/include)

include_directories(${PROTO_INC} ${PDFL_INC} ${GTEST_INC})

set(PROTO_SRC
        examples/Euler/src/EulerOp.cpp
        examples/Euler/src/EulerOp.H
        include/Proto.H
        include/Proto_Box.H
        include/Proto_BoxData.H
        include/Proto_BoxDataImplem.H
        include/Proto_BoxImplem.H
        include/Proto_ClockTicks.H
        include/Proto_DebugHooks.H
        include/Proto_DebugHooksImplem.H
        include/Proto_DisjointBoxLayout.H
        include/Proto_DisjointBoxLayoutI.H
        include/Proto_LevelData.H
        include/Proto_LevelDataI.H
        include/Proto_PAssert.H
        include/Proto_Point.H
        include/Proto_PointImplem.H
        include/Proto_PowerItoI.H
        include/Proto_RK4.H
        include/Proto_Stack.H
        include/Proto_StackImplem.H
        include/Proto_Stencil.H
        include/Proto_StencilDefs.H
        include/Proto_StencilImplem.H
        include/Proto_Timer.H
        include/Proto_TimerImplem.H
        include/Proto_VisitWriter.H
        include/Proto_VisitWriterImplem.H
        include/Proto_WriteBoxData.H
        include/Proto_WriteBoxDataImplem.H
        include/Proto_DataFlow.H)

set(IEGEN_PATH pdfl/lib/iegenlib/src)
set(IEGEN_SRC  ${IEGEN_PATH}/set_relation/environment.cc
        ${IEGEN_PATH}/set_relation/expression.cc
        ${IEGEN_PATH}/set_relation/isl_str_manipulation.cc
        ${IEGEN_PATH}/set_relation/set_relation.cc
        ${IEGEN_PATH}/set_relation/SubMap.cc
        ${IEGEN_PATH}/set_relation/TupleDecl.cc
        ${IEGEN_PATH}/set_relation/UFCallMap.cc
        ${IEGEN_PATH}/set_relation/UninterpFunc.cc
        ${IEGEN_PATH}/util/util.cc
        ${IEGEN_PATH}/parser/gen_scanner.cc
        ${IEGEN_PATH}/parser/gen_parser.cc
        ${IEGEN_PATH}/parser/parser.cc
)

set(ISL_PATH pdfl/lib/isl)
set(ISL_SRC ${ISL_PATH}/imath/gmp_compat.c
        ${ISL_PATH}/imath/imath.c
        ${ISL_PATH}/imath/imrat.c
        ${ISL_PATH}/isl_map.c
        ${ISL_PATH}/isl_input.c
        ${ISL_PATH}/isl_output.c
        ${ISL_PATH}/isl_printer.c
        ${ISL_PATH}/isl_coalesce.c
        ${ISL_PATH}/isl_union_map.c
        ${ISL_PATH}/isl_map_simplify.c
        ${ISL_PATH}/isl_map_subtract.c
        ${ISL_PATH}/isl_ctx.c
        ${ISL_PATH}/isl_space.c
        ${ISL_PATH}/isl_mat.c
        ${ISL_PATH}/isl_local_space.c
        ${ISL_PATH}/isl_dim_map.c
        ${ISL_PATH}/isl_map_list.c
        ${ISL_PATH}/isl_map_simplify.c
        ${ISL_PATH}/isl_val.c
        ${ISL_PATH}/isl_int_sioimath.h
        ${ISL_PATH}/isl_int_sioimath.c
        ${ISL_PATH}/isl_val_sioimath.c
        ${ISL_PATH}/isl_imath.c
        ${ISL_PATH}/isl_obj.c
        ${ISL_PATH}/isl_affine_hull.c
        ${ISL_PATH}/isl_convex_hull.c
        ${ISL_PATH}/isl_bernstein.c
        ${ISL_PATH}/isl_blk.c
        ${ISL_PATH}/isl_bound.c
        ${ISL_PATH}/isl_id.c
        ${ISL_PATH}/isl_id_to_ast_expr.c
        ${ISL_PATH}/isl_seq.c
        ${ISL_PATH}/isl_equalities.c
        ${ISL_PATH}/isl_factorization.c
        ${ISL_PATH}/isl_hash.c
        ${ISL_PATH}/isl_tab.c
        ${ISL_PATH}/isl_tab_pip.c
        ${ISL_PATH}/isl_transitive_closure.c
        ${ISL_PATH}/isl_vec.c
        ${ISL_PATH}/isl_aff.c
        ${ISL_PATH}/isl_reordering.c
        ${ISL_PATH}/isl_constraint.c
        ${ISL_PATH}/isl_lp.c
        ${ISL_PATH}/isl_polynomial.c
        ${ISL_PATH}/isl_sort.c
        ${ISL_PATH}/isl_farkas.c
        ${ISL_PATH}/isl_tarjan.c
        ${ISL_PATH}/isl_fold.c
        ${ISL_PATH}/isl_morph.c
        ${ISL_PATH}/isl_arg.c
        ${ISL_PATH}/isl_options.c
        ${ISL_PATH}/isl_point.c
        ${ISL_PATH}/isl_range.c
        ${ISL_PATH}/isl_sample.c
        ${ISL_PATH}/isl_scan.c
        ${ISL_PATH}/isl_stream.c
        ${ISL_PATH}/isl_version.c
        ${ISL_PATH}/isl_vertices.c
        ${ISL_PATH}/basis_reduction_tab.c
        ${ISL_PATH}/isl_set_list.c
        ${ISL_PATH}/isl_ast.c
        ${ISL_PATH}/isl_ast_build.c
        ${ISL_PATH}/isl_ast_build_expr.c
        ${ISL_PATH}/isl_ast_codegen.c
        ${ISL_PATH}/isl_ast_graft.c
        ${ISL_PATH}/isl_band.c
        ${ISL_PATH}/isl_ilp.c
        ${ISL_PATH}/isl_schedule.c
        ${ISL_PATH}/isl_schedule_node.c
        ${ISL_PATH}/isl_schedule_tree.c
        ${ISL_PATH}/isl_schedule_band.c
        ${ISL_PATH}/isl_schedule_read.c
        ${ISL_PATH}/print.c
)

set(OMEGA_PATH pdfl/lib/chill/omega)
set(OMEGA_SRC ${OMEGA_PATH}/basic/src/ConstString.cc
        ${OMEGA_PATH}/parser/parser.tab.cc
        ${OMEGA_PATH}/parser/lex.yy.cc
        ${OMEGA_PATH}/parser/AST.cc
        ${OMEGA_PATH}/omega_lib/src/closure.cc
        ${OMEGA_PATH}/omega_lib/src/farkas.cc
        ${OMEGA_PATH}/omega_lib/src/hull.cc
        ${OMEGA_PATH}/omega_lib/src/pres_beaut.cc
        ${OMEGA_PATH}/omega_lib/src/pres_col.cc
        ${OMEGA_PATH}/omega_lib/src/pres_dnf.cc
        ${OMEGA_PATH}/omega_lib/src/pres_cnstr.cc
        ${OMEGA_PATH}/omega_lib/src/pres_decl.cc
        ${OMEGA_PATH}/omega_lib/src/pres_gen.cc
        ${OMEGA_PATH}/omega_lib/src/pres_logic.cc
        ${OMEGA_PATH}/omega_lib/src/pres_form.cc
        ${OMEGA_PATH}/omega_lib/src/pres_conj.cc
        ${OMEGA_PATH}/omega_lib/src/pres_print.cc
        ${OMEGA_PATH}/omega_lib/src/pres_quant.cc
        ${OMEGA_PATH}/omega_lib/src/pres_rear.cc
        ${OMEGA_PATH}/omega_lib/src/pres_var.cc
        ${OMEGA_PATH}/omega_lib/src/reach.cc
        ${OMEGA_PATH}/omega_lib/src/hull_simple.cc
        ${OMEGA_PATH}/omega_lib/src/Relation.cc
        ${OMEGA_PATH}/omega_lib/src/RelBody.cc
        ${OMEGA_PATH}/omega_lib/src/Relations.cc
        ${OMEGA_PATH}/omega_lib/src/RelVar.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_eq.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_exp_kill.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_global.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_print.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_problems.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_query.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_quick_kill.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_simple.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_solve.cc
        ${OMEGA_PATH}/omega_lib/src/omega_core/oc_util.cc
        ${OMEGA_PATH}/code_gen/src/codegen.cc
        ${OMEGA_PATH}/code_gen/src/CG.cc
        ${OMEGA_PATH}/code_gen/src/CG_utils.cc
        ${OMEGA_PATH}/code_gen/src/CG_stringBuilder.cc
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fPIC -DPIC -O3 -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-deprecated-declarations -Wno-format-overflow")

#unset(TESTING)
set(TESTING "true")

if (TESTING)
    set(GTEST_PATH pdfl/lib/gtest/src)
    set(GTEST_SRC ${GTEST_PATH}/gtest.cc
                  ${GTEST_PATH}/gtest-filepath.cc
                  ${GTEST_PATH}/gtest-death-test.cc
                  ${GTEST_PATH}/gtest-test-part.cc
                  ${GTEST_PATH}/gtest-printers.cc
                  ${GTEST_PATH}/gtest-port.cc)

    set(UNIT_TESTS test/EulerTest.cpp test/driver.cpp)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -O3")

    enable_testing()
    add_executable(eulerTest ${UNIT_TESTS} ${GTEST_SRC} ${PROTO_SRC} ${PDFL_SRC} ${IEGEN_SRC} ${ISL_SRC} ${OMEGA_SRC})
    target_link_libraries(eulerTest m omp)
    target_compile_definitions(eulerTest PRIVATE DIM=3 NUMCELLS=64 DATAFLOW_ON=0)
    add_test(eulerTest driver EulerTest)
else()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
    add_executable(euler examples/Euler/exec/Euler.cpp
            examples/Euler/src/EulerRK4.cpp
            examples/Euler/src/EulerRK4.H
            ${PROTO_SRC} ${IEGEN_SRC} ${ISL_SRC} ${OMEGA_SRC})
    target_compile_definitions(euler PRIVATE DIM=3 NUMCELLS=64 DATAFLOW_ON=1)
endif()
