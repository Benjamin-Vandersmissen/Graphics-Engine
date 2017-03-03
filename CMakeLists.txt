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
        Intro.cc
        Intro.hh
        l_parser.cc
        l_parser.hh
        Line2D.cc
        Line2D.hh
        LSystems.cc
        LSystems.hh
        Point2D.cc
        Point2D.hh
        UsefulFunctions.cc
        UsefulFunctions.hh
        vector.cc
        vector.hh Figure3D.cc Figure3D.hh WireFrameParser.cc WireFrameParser.hh)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})