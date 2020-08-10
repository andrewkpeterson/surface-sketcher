#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QList>
#include <QMainWindow>
#include <QToolBar>
#include "src/pipeline/Pipeline.h"

class ScribbleArea;

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void open();
    void save();
    void penColor();
    void penWidth();
    void about();
    void runBendsketch();

private:
    void createActions();
    void createMenus();
    bool maybeSave();
    bool saveFile(const QByteArray &fileFormat);

    ScribbleArea *scribbleArea;

    QMenu *saveAsMenu;
    QMenu *fileMenu;
    QMenu *optionMenu;
    QMenu *helpMenu;
    QMenu *strokeMenu;

    QAction *openAct;
    QList<QAction *> saveAsActs;
    QAction *exitAct;
    QAction *penColorAct;
    QAction *penWidthAct;
    QAction *printAct;
    QAction *clearScreenAct;
    QAction *aboutAct;
    QAction *aboutQtAct;

    QAction *setToBoundary;
    QAction *setToConcave;
    QAction *setToConvex;
    QAction *runBendsketchAction;

    Ui::MainWindow *ui;
};

#endif
