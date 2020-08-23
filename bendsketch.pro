QT       += core gui xml

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

SOURCES += \
    ScribbleArea.cpp \
    main.cpp \
    mainwindow.cpp \
    src/pipeline/DirectionFieldInitializer.cpp \
    src/pipeline/DirectionFieldOptimizer.cpp \
    src/pipeline/HeightFieldSolver.cpp \
    src/pipeline/Mesh.cpp \
    src/pipeline/Sketch.cpp \
    src/pipeline/Triangulate.cpp \
    $$files(alglib/*.cpp, true)

HEADERS += \
    ScribbleArea.h \
    mainwindow.h \
    src/pipeline/DirectionFieldInitializer.h \
    src/pipeline/DirectionFieldOptimizer.h \
    src/pipeline/HeightFieldSolver.h \
    src/pipeline/Mesh.h \
    src/pipeline/Pipeline.h \
    src/pipeline/Sketch.h \
    src/pipeline/Triangulate.h \
    src/utils/OBJWriter.h \
    nanosvg/src/nanosvg.h \
    $$files(alglib/*.h, true)

INCLUDEPATH += /usr/local/include/eigen3 \
               /usr/local/include/CGAL \
               /usr/local/include/boost \

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/ceres-solver/1.14.0_12/lib/release/ -lceres.1.14.0
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/ceres-solver/1.14.0_12/lib/debug/ -lceres.1.14.0
else:unix: LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/ceres-solver/1.14.0_12/lib/ -lceres.1.14.0

INCLUDEPATH += $$PWD/../../../../../../../usr/local/Cellar/ceres-solver/1.14.0_12/include
DEPENDPATH += $$PWD/../../../../../../../usr/local/Cellar/ceres-solver/1.14.0_12/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/glog/0.4.0/lib/release/ -lglog.0.4.0
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/glog/0.4.0/lib/debug/ -lglog.0.4.0
else:unix: LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/glog/0.4.0/lib/ -lglog.0.4.0

INCLUDEPATH += $$PWD/../../../../../../../usr/local/Cellar/glog/0.4.0/include
DEPENDPATH += $$PWD/../../../../../../../usr/local/Cellar/glog/0.4.0/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/gflags/2.2.2/lib/release/ -lgflags.2.2.2
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/gflags/2.2.2/lib/debug/ -lgflags.2.2.2
else:unix: LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/gflags/2.2.2/lib/ -lgflags.2.2.2

INCLUDEPATH += $$PWD/../../../../../../../usr/local/Cellar/gflags/2.2.2/include
DEPENDPATH += $$PWD/../../../../../../../usr/local/Cellar/gflags/2.2.2/include
