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
    void runSurfaceSketcher();

    void on_boundary_stroke_clicked();

    void on_contour_stroke_clicked();

    void on_start_height_valueChanged(double arg1);

    void on_end_height_valueChanged(double arg1);

    void on_concave_stroke_2_clicked();

    void on_concave_stroke_clicked();

    void on_radius_valueChanged(double arg1);

    void on_run_surface_sketcher_clicked();

private:
    void createActions();
    void createMenus();
    bool maybeSave();
    bool saveFile(const QByteArray &fileFormat);

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
    QAction *runAction;

    Ui::MainWindow *ui;
};

#endif
