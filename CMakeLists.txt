cmake_minimum_required(VERSION 3.6)
set(PROJECT_NAME  engine)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        easy_image.cc
        easy_image.hh
        engine.cc
        ini_configuration.cc
        ini_configuration.hh
        IntroParser.cc
        IntroParser.hh
        l_parser.cc
        l_parser.hh
        Line2D.cc
        Line2D.hh
        LSystems.cc
        LSystems.hh
        UsefulFunctions.cc
        UsefulFunctions.hh
        vector.cc
        vector.hh Figure3D.cc Figure3D.hh WireFrameParser.cc WireFrameParser.hh Light.cc Light.hh)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})