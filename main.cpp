#include "src/pipeline/Pipeline.h"


/*int main(int argc, char *argv[])
{
    Pipeline::runPipeline(argv[1], argv[2]);
} */

#include <QApplication>

#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    MainWindow window;
#if defined(Q_OS_SYMBIAN)
    window.showMaximized();
#else
    window.show();
#endif
    return app.exec();
}
