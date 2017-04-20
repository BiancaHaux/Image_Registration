project(MeanSquareMetric)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
add_executable(MeanSquareMetric main.cpp)
 
target_link_libraries(MeanSquareMetric ${ITK_LIBRARIES})
