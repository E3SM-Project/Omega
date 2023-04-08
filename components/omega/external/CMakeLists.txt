# Add external packages
#add_subdirectory(external_package1)

# Add the spdlog header-only library
add_library(spdlog INTERFACE)
target_include_directories(spdlog INTERFACE spdlog/include)
