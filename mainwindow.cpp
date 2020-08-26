#include "ui_mainwindow.h"

#include <QtGui>
#include <QFileDialog>
#include <QColorDialog>
#include <QMessageBox>
#include <QInputDialog>

#include "mainwindow.h"
#include "ScribbleArea.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    scribbleArea = new ScribbleArea;
    ui->window->layout()->addWidget(scribbleArea);
    //setCentralWidget(scribbleArea);

    createActions();
    createMenus();

    setWindowTitle(tr("bendsketch-lite"));
    resize(500, 500);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::runBendsketch() {
    Pipeline::runPipelineScribble(scribbleArea->getSketchData());
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    //if (maybeSave()) {
        event->accept();
    //} else {
        //event->ignore();
    //}
}

void MainWindow::open()
{
    if (maybeSave()) {
        QString fileName = QFileDialog::getOpenFileName(this,
                                   tr("Open File"), QDir::currentPath());
        if (!fileName.isEmpty())
            scribbleArea->openImage(fileName);
    }
}

void MainWindow::save()
{
    QAction *action = qobject_cast<QAction *>(sender());
    QByteArray fileFormat = action->data().toByteArray();
    saveFile(fileFormat);
}

void MainWindow::penColor()
{
    QColor newColor = QColorDialog::getColor(scribbleArea->penColor());
    if (newColor.isValid())
        scribbleArea->setPenColor(newColor);
}

void MainWindow::penWidth()
{
    bool ok;
    int newWidth = QInputDialog::getInt(this, tr("Scribble"),
                                            tr("Select pen width:"),
                                            scribbleArea->penWidth(),
                                            1, 50, 1, &ok);
    if (ok)
        scribbleArea->setPenWidth(newWidth);
}

void MainWindow::about()
{

}

void MainWindow::createActions()
{
    openAct = new QAction(tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

    foreach (QByteArray format, QImageWriter::supportedImageFormats()) {
        QString text = tr("%1...").arg(QString(format).toUpper());

        QAction *action = new QAction(text, this);
        action->setData(format);
        connect(action, SIGNAL(triggered()), this, SLOT(save()));
        saveAsActs.append(action);
    }

    printAct = new QAction(tr("&Print..."), this);
    connect(printAct, SIGNAL(triggered()), scribbleArea, SLOT(print()));

    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcuts(QKeySequence::Quit);
    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

    penColorAct = new QAction(tr("&Pen Color..."), this);
    connect(penColorAct, SIGNAL(triggered()), this, SLOT(penColor()));

    penWidthAct = new QAction(tr("Pen &Width..."), this);
    connect(penWidthAct, SIGNAL(triggered()), this, SLOT(penWidth()));

    clearScreenAct = new QAction(tr("&Clear Screen"), this);
    clearScreenAct->setShortcut(tr("Ctrl+L"));
    connect(clearScreenAct, SIGNAL(triggered()),
            scribbleArea, SLOT(clearImage()));

    aboutAct = new QAction(tr("&About"), this);
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));

    aboutQtAct = new QAction(tr("About &Qt"), this);
    connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));

    setToBoundary = new QAction(tr("&Boundary Stroke"), this);
    connect(setToBoundary, SIGNAL(triggered()), scribbleArea, SLOT(setStrokeTypeToBoundary()));

    setToConvex = new QAction(tr("&Convex Bending Stroke"), this);
    connect(setToConvex, SIGNAL(triggered()), scribbleArea, SLOT(setStrokeTypeToConvexBend()));

    setToConcave = new QAction(tr("&Concave Bending Stroke"), this);
    connect(setToConcave, SIGNAL(triggered()), scribbleArea, SLOT(setStrokeTypeToConcaveBend()));

    runBendsketchAction = new QAction(tr("&Run Bendsketch"), this);
    connect(runBendsketchAction, SIGNAL(triggered()), this, SLOT(runBendsketch()));
}

void MainWindow::createMenus()
{
    saveAsMenu = new QMenu(tr("&Save As"), this);
    foreach (QAction *action, saveAsActs)
        saveAsMenu->addAction(action);

    fileMenu = new QMenu(tr("&File"), this);
    fileMenu->addAction(openAct);
    fileMenu->addMenu(saveAsMenu);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    optionMenu = new QMenu(tr("&Options"), this);
    optionMenu->addAction(clearScreenAct);

    helpMenu = new QMenu(tr("&Help"), this);
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);

    strokeMenu = new QMenu(tr("&Bendsketch Actions"), this);
    strokeMenu->addAction(setToBoundary);
    strokeMenu->addAction(setToConvex);
    strokeMenu->addAction(setToConcave);
    strokeMenu->addAction(runBendsketchAction);

    menuBar()->addMenu(fileMenu);
    menuBar()->addMenu(optionMenu);
    menuBar()->addMenu(helpMenu);
    menuBar()->addMenu(strokeMenu);


}

bool MainWindow::maybeSave()
{
    if (scribbleArea->isModified()) {
       QMessageBox::StandardButton ret;
       ret = QMessageBox::warning(this, tr("Scribble"),
                          tr("The image has been modified.\n"
                             "Do you want to save your changes?"),
                          QMessageBox::Save | QMessageBox::Discard
                          | QMessageBox::Cancel);
        if (ret == QMessageBox::Save) {
            return saveFile("png");
        } else if (ret == QMessageBox::Cancel) {
            return false;
        }
    }
    return true;
}

bool MainWindow::saveFile(const QByteArray &fileFormat)
{
    QString initialPath = QDir::currentPath() + "/untitled." + fileFormat;

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"),
                               initialPath,
                               tr("%1 Files (*.%2);;All Files (*)")
                               .arg(QString(fileFormat.toUpper()))
                               .arg(QString(fileFormat)));
    if (fileName.isEmpty()) {
        return false;
    } else {
        return scribbleArea->saveImage(fileName, fileFormat);
    }
}



void MainWindow::on_boundary_stroke_clicked()
{
    scribbleArea->setStrokeTypeToBoundary();
}

void MainWindow::on_contour_stroke_clicked()
{
    scribbleArea->setStrokeTypeToContour();
}

void MainWindow::on_start_height_valueChanged(double arg1)
{
    scribbleArea->setStartHeight(arg1);
}

void MainWindow::on_end_height_valueChanged(double arg1)
{
    scribbleArea->setEndHeight(arg1);
}

void MainWindow::on_concave_stroke_2_clicked()
{
    scribbleArea->setStrokeTypeToConvexBend();
}

void MainWindow::on_concave_stroke_clicked()
{
    scribbleArea->setStrokeTypeToConcaveBend();
}

void MainWindow::on_run_bendsketch_clicked()
{
    Pipeline::runPipelineScribble(scribbleArea->getSketchData());
}
