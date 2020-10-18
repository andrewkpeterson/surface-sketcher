QT += core gui xml

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

EIGEN_PATH = # path to eigen3 directory inside of the include directory of eigen installation (i.e. /usr/local/Cellar/eigen/3.3.7/include/eigen3)
CGAL_PATH = # path to CGAL installation (e.g. /usr/local/Cellar/CGAL/5.1)
BOOST_PATH = # path to boost installation (e.g. /usr/local/Cellar/boost/1.73.0)
CERES_PATH = # path to ceres-solver installation (e.g. /usr/local/Cellar/ceres-solver/1.14.0_12)
GLOG_PATH = # path to glog installation (e.g. /usr/local/Cellar/glog/0.4.0)
GFLAGS_PATH = # path to gflags installation (e.g. /usr/local/Cellar/gflags/2.2.2)


SOURCES += \
    src/drawing/ScribbleArea.cpp \
    main.cpp \
    mainwindow.cpp \
    src/pipeline/DirectionFieldInitializer.cpp \
    src/pipeline/DirectionFieldOptimizer.cpp \
    src/pipeline/HeightFieldSolver.cpp \
    src/pipeline/Mesh.cpp \
    src/pipeline/Sketch.cpp \
    src/pipeline/Triangulate.cpp \

HEADERS += \
    src/drawing/ScribbleArea.h \
    mainwindow.h \
    src/pipeline/DirectionFieldInitializer.h \
    src/pipeline/DirectionFieldOptimizer.h \
    src/pipeline/HeightFieldSolver.h \
    src/pipeline/Mesh.h \
    src/pipeline/Pipeline.h \
    src/pipeline/Sketch.h \
    src/pipeline/Triangulate.h \
    src/utils/OBJWriter.h \

INCLUDEPATH = $${EIGEN_PATH} \
               $${CGAL_PATH}/include \
               $${BOOST_PATH}/include \

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

LIBS += -L$$CERES_PATH/lib/ -lceres
INCLUDEPATH += $$CERES_PATH/include
DEPENDPATH += $$CERES_PATH/include

LIBS += -L$$GLOG_PATH/lib/ -lglog
INCLUDEPATH += $$GLOG_PATH/include
DEPENDPATH += $$GLOG_PATH/include

LIBS += -L$$GFLAGS_PATH/lib/ -lgflags
INCLUDEPATH += $$GFLAGS_PATH/include
DEPENDPATH += $$GFLAGS_PATH/include
