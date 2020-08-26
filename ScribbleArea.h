#ifndef SCRIBBLEAREA_H
#define SCRIBBLEAREA_H

#include "ui_mainwindow.h"
#include <QColor>
#include <QImage>
#include <QPoint>
#include <QWidget>
#include "Eigen/Dense"

namespace Ui {
    class ScribbleArea;
}

using Stroke = std::vector<Eigen::Vector2f>;

struct BoundaryStroke {
    std::vector<Eigen::Vector2f> points;
    std::vector<float> heights;
    bool contour = false;
};

struct SketchData {
    std::vector<BoundaryStroke> boundary;
    std::vector<BoundaryStroke> contour;

    std::vector<Stroke> convex;
    std::vector<Stroke> concave;

    float width;
    float height;
};

class ScribbleArea : public QWidget
{
    Q_OBJECT

public:
    enum class StrokeType {
        BOUNDARY, CONTOUR, CONVEX_BEND, CONCAVE_BEND
    };

    ScribbleArea(QWidget *parent = 0);

    bool openImage(const QString &fileName);
    bool saveImage(const QString &fileName, const char *fileFormat);
    void setPenColor(const QColor &newColor);
    void setPenWidth(int newWidth);
    void setStartHeight(float height) { start_height = height; }
    void setEndHeight(float height) { end_height = height; }

    bool isModified() const { return modified; }
    QColor penColor() const { return myPenColor; }
    int penWidth() const { return myPenWidth; }

    SketchData& getSketchData() { return data; }

public slots:
    void clearImage();
    void print();
    void setStrokeTypeToBoundary() {
        m_type = StrokeType::BOUNDARY;
        myPenColor = QColor(0,0,0);
    }
    void setStrokeTypeToContour() {
        m_type = StrokeType::CONTOUR;
        myPenColor = QColor(0,0,0);
    }
    void setStrokeTypeToConvexBend() {
        m_type = StrokeType::CONVEX_BEND;
        myPenColor = QColor(255,0,0);
    }
    void setStrokeTypeToConcaveBend() {
        m_type = StrokeType::CONCAVE_BEND;
        myPenColor = QColor(0,0,255);
    }

protected:
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void paintEvent(QPaintEvent *event);
    void resizeEvent(QResizeEvent *event);

private:
    void drawLineTo(const QPoint &endPoint);
    void resizeImage(QImage *image, const QSize &newSize);
    void addStrokeToData();
    BoundaryStroke createBoundaryStroke();

    bool modified;
    bool scribbling;
    int myPenWidth;
    QColor myPenColor;
    QImage image;
    QPoint lastPoint;
    StrokeType m_type;

    SketchData data;
    std::vector<Eigen::Vector2f> currentStroke;
    int strokeCount = 0;

    float start_height = 0;
    float end_height = 0;

    float DISTANCE_BETWEEN_POINTS = 2;

    Ui::ScribbleArea *ui;
};

#endif
