TEMPLATE = app
CONFIG += console
CONFIG += c++1z
CONFIG -= app_bundle
CONFIG -= qt

#QMAKE_CXXFLAGS += -Wall -O3 -fopenmp

CONFIG(release) { #// , debug|release
	CONFIG += optimize_full
	QMAKE_CXXFLAGS += -fopenmp
	LIBS += -fopenmp
}




INCLUDEPATH=/usr/include/eigen3

SOURCES += \
        test.cpp

HEADERS += \
        matrix.h \
		timer.h
