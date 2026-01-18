# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(orderbook LANGUAGES CXX)

# ---- Options (match your script knobs) ----
option(ENABLE_SANITIZERS "Enable ASan+UBSan (Debug only; best with Clang)" OFF)
option(ENABLE_TIDY       "Enable clang-tidy (if available)" ON)
option(TIDY_FIX          "Apply clang-tidy fixes (may rewrite files)" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Produce compile_commands.json for clangd/clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---- Target ----
add_executable(main
  src/main.cpp
)

# ---- Common compile options ----
# Equivalent to: -fcolor-diagnostics (Clang) and strong warnings + Werror
target_compile_options(main PRIVATE
  $<$<CXX_COMPILER_ID:Clang>:-fcolor-diagnostics>

  -Wall
  -Wextra
  -Wpedantic
  -Werror
  -Wshadow
  -Wconversion
  -Wsign-conversion
  -Wfloat-conversion
  -Wdouble-promotion
  -Wformat=2
  -Wimplicit-fallthrough
  -Wextra-semi
  -Woverloaded-virtual
  -Wnon-virtual-dtor
  -Wold-style-cast
  -Wcast-qual
  -Wnull-dereference
  -Wundef
  -Wunreachable-code
  $<$<CXX_COMPILER_ID:Clang>:-Wrange-loop-analysis>
  $<$<CXX_COMPILER_ID:Clang>:-Wloop-analysis>
)

# ---- Debug vs Release flags ----
# CMake already sets -O0/-g for Debug and -O3/-DNDEBUG for Release typically,
# but we force the ones you used, including -fno-omit-frame-pointer in Debug.
target_compile_options(main PRIVATE
  $<$<CONFIG:Debug>:-O0 -g -fno-omit-frame-pointer>
  $<$<CONFIG:Release>:-O3>
)
target_compile_definitions(main PRIVATE
  $<$<CONFIG:Release>:NDEBUG>
)

# ---- Sanitizers (Debug-only) ----
# Matches: -fsanitize=address,undefined -fno-sanitize-recover=all
if(ENABLE_SANITIZERS)
  target_compile_options(main PRIVATE
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:Clang>>:-fsanitize=address,undefined -fno-sanitize-recover=all>
  )
  target_link_options(main PRIVATE
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:Clang>>:-fsanitize=address,undefined -fno-sanitize-recover=all>
  )
endif()

# ---- clang-tidy integration ----
# Notes:
# - This runs clang-tidy via CMake's CXX_CLANG_TIDY for the target (per-file).
# - "TIDY_FIX" maps to -fix -fix-errors.
# - header filter is best-effort; paths are relative to source dir.
if(ENABLE_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    set(_tidy_checks
      "misc-const-correctness,readability-make-member-function-const,readability-const-return-type"
    )

    set(_tidy_args
      "${CLANG_TIDY_EXE}"
      "-checks=${_tidy_checks}"
      "-warnings-as-errors=*"
      "-header-filter=^(.*/)?(src|include)/"
      "-quiet"
    )

    if(TIDY_FIX)
      list(APPEND _tidy_args "-fix" "-fix-errors")
    endif()

    set_target_properties(main PROPERTIES CXX_CLANG_TIDY "${_tidy_args}")
  else()
    message(STATUS "clang-tidy not found; configure with -DENABLE_TIDY=OFF to silence this message.")
  endif()
endif()

# ---- Convenience custom targets ----
# Mirrors your script usage patterns.
add_custom_target(run
  COMMAND main
  DEPENDS main
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# ---- Usage (examples) ----
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
#   cmake --build build -j
#   ./build/main
#
#   cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-rel -j
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON
#   cmake --build build -j
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_TIDY=ON -DTIDY_FIX=ON
#   cmake --build build -j