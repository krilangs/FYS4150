TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        funcs.cpp \
        main.cpp

HEADERS += \
    funcs.h

QMAKE_CXXFLAGS += -Ofast
//QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

INCLUDEPATH += C:\armadillo-9.600.6\include
DEPENDPATH += C:\armadillo-9.600.6\include

LIBS += \
    -LC:\armadillo-9.600.6\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT

