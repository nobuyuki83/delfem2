cmake_minimum_required(VERSION 3.12)

project(examples_cuda)

add_subdirectory(dfm2_cuda_static)

add_subdirectory(000_Simple)
add_subdirectory(001_ThrustTest)
add_subdirectory(002_MatMat)
add_subdirectory(003_VecDot)
add_subdirectory(004_MinMaxPoints3D)
add_subdirectory(005_ImageGrayscaleSobel)
