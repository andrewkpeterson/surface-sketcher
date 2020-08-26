#include <QtGui>

#include "ScribbleArea.h"
#include <iostream>
#include "src/pipeline/Sketch.h"

ScribbleArea::ScribbleArea(QWidget *parent)
    : QWidget(parent)
{
    setAttribute(Qt::WA_StaticContents);
    modified = false;
    scribbling = false;
    myPenWidth = 2;
    myPenColor = Qt::black;
    m_type = StrokeType::BOUNDARY;
    data.width = width();
    data.height = height();
}

bool ScribbleArea::openImage(const QString &fileName)
{
    QImage loadedImage;
    if (!loadedImage.load(fileName))
        return false;

    int idx = fileName.lastIndexOf(".");
    QString text = fileName.mid(0, idx);
    QFile file(text + ".txt");
    if(!file.open(QIODevice::ReadOnly)) {
        return false;
    }

    data.convex.clear();
    data.boundary.clear();
    data.concave.clear();
    QTextStream in(&file);
    QString line = in.readLine();
    BoundaryStroke curr_boundary;
    while (line != "contour") {
        if (QString::compare(line,QString("boundary")) != 0 && QString::compare(line, QString("")) != 0) {
            float x;
            float y;
            float height;
            std::sscanf(line.toUtf8().constData(), "%f %f %f", &x, &y, &height);
            curr_boundary.points.push_back(Eigen::Vector2f(x,y));
            curr_boundary.heights.push_back(height);
        } else if (curr_boundary.points.size() > 0) {
            curr_boundary.contour = false;
            data.boundary.push_back(curr_boundary);
            curr_boundary = BoundaryStroke();
        }
        line = in.readLine();
    }
    if (curr_boundary.points.size() > 0) {
        curr_boundary.contour = false;
        data.boundary.push_back(curr_boundary);
        curr_boundary = BoundaryStroke();
    }

    while (line != "convex") {
        if (QString::compare(line,QString("contour")) != 0 && QString::compare(line, QString("")) != 0) {
            float x;
            float y;
            float height;
            std::sscanf(line.toUtf8().constData(), "%f %f %f", &x, &y, &height);
            curr_boundary.points.push_back(Eigen::Vector2f(x,y));
            curr_boundary.heights.push_back(height);
        } else if (curr_boundary.points.size() > 0) {
            curr_boundary.contour = true;
            data.contour.push_back(curr_boundary);
            curr_boundary = BoundaryStroke();
        }
        line = in.readLine();
    }
    if (curr_boundary.points.size() > 0) {
        curr_boundary.contour = true;
        data.contour.push_back(curr_boundary);
        curr_boundary = BoundaryStroke();
    }

    std::vector<Eigen::Vector2f> curr_vec = std::vector<Eigen::Vector2f>();
    while (line != "concave") {
        if (line.compare("convex") != 0 && line.compare("") != 0) {
            float x;
            float y;
            std::sscanf(line.toUtf8().constData(), "%f %f", &x, &y);
            curr_vec.push_back(Eigen::Vector2f(x,y));
        } else if (curr_vec.size() > 0) {
            data.convex.push_back(curr_vec);
            curr_vec = std::vector<Eigen::Vector2f>();
        }
        line = in.readLine();
    }
    if (curr_vec.size() > 0) {
        data.convex.push_back(curr_vec);
        curr_vec = std::vector<Eigen::Vector2f>();
    }

    curr_vec = std::vector<Eigen::Vector2f>();
    while (!in.atEnd()) {
        if (line.compare("concave") != 0 && line.compare("") != 0) {
            float x;
            float y;
            std::sscanf(line.toUtf8().constData(), "%f %f", &x, &y);
            curr_vec.push_back(Eigen::Vector2f(x,y));
        } else if (curr_vec.size() > 0) {
            data.concave.push_back(curr_vec);
            curr_vec = std::vector<Eigen::Vector2f>();
        }
        line = in.readLine();
    }
    if (curr_vec.size() > 0) {
        data.concave.push_back(curr_vec);
    }

    file.close();


    QSize newSize = loadedImage.size().expandedTo(size());
    resizeImage(&loadedImage, newSize);
    image = loadedImage;
    modified = false;
    update();
    return true;
}

bool ScribbleArea::saveImage(const QString &fileName, const char *fileFormat)
{
    QImage visibleImage = image;
    resizeImage(&visibleImage, size());

    if (visibleImage.save(fileName, fileFormat)) {
        modified = false;
    }
    int idx = fileName.lastIndexOf(".");
    QString text = fileName.mid(0, idx);
    QFile file(text + ".txt");
    if (file.open(QIODevice::ReadWrite | QFile::Truncate)) {
        QTextStream stream(&file);
        char buf[512];
        stream << "boundary" << endl;
        for (int i = 0; i < data.boundary.size(); i++) {
            for (int j = 0; j < data.boundary[i].points.size(); j++) {
                std::sprintf(buf, "%f %f %f", data.boundary[i].points[j].x(), data.boundary[i].points[j].y(), data.boundary[i].heights[j]);
                stream << buf << endl;
            }
            stream << endl;
        }
        stream << "contour" << endl;
        for (int i = 0; i < data.contour.size(); i++) {
            for (int j = 0; j < data.contour[i].points.size(); j++) {
                std::sprintf(buf, "%f %f %f", data.contour[i].points[j].x(), data.contour[i].points[j].y(), data.contour[i].heights[j]);
                stream << buf << endl;
            }
            stream << endl;
        }
        stream << "convex" << endl;
        for (int i = 0; i < data.convex.size(); i++) {
            for (int j = 0; j < data.convex[i].size(); j++) {
                std::sprintf(buf, "%f %f", data.convex[i][j].x(), data.convex[i][j].y());
                stream << buf << endl;
            }
            stream << endl;
        }
        stream << "concave" << endl;
        for (int i = 0; i < data.concave.size(); i++) {
            for (int j = 0; j < data.concave[i].size(); j++) {
                std::sprintf(buf, "%f %f", data.concave[i][j].x(), data.concave[i][j].y());
                stream << buf << endl;
            }
            stream << endl;
        }
    }

    return true;
}

void ScribbleArea::setPenColor(const QColor &newColor)
{
    myPenColor = newColor;
}

void ScribbleArea::setPenWidth(int newWidth)
{
    myPenWidth = newWidth;
}

void ScribbleArea::clearImage()
{
    image.fill(qRgb(255, 255, 255));
    modified = true;
    update();
}

void ScribbleArea::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        lastPoint = event->pos();
        scribbling = true;
    }
}

void ScribbleArea::mouseMoveEvent(QMouseEvent *event)
{
    if ((event->buttons() & Qt::LeftButton) && scribbling)
        drawLineTo(event->pos());
}

void ScribbleArea::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton && scribbling) {
        drawLineTo(event->pos());
        addStrokeToData();
        scribbling = false;
    }
}

BoundaryStroke ScribbleArea::createBoundaryStroke() {
    BoundaryStroke stroke;
    stroke.points = currentStroke;
    if (m_type == StrokeType::CONTOUR) {
        stroke.contour = true;
    }
    float L = end_height - start_height;
    float stroke_length = 0;
    float k = 1;
    for (int i = 1; i < currentStroke.size(); i++) {
        stroke_length += (currentStroke[i-1] - currentStroke[i]).norm();
    }

    float accumulated_length = 0;
    for (int i = 0; i < currentStroke.size() - 1; i++) {
        float t = 12 * (.5 * accumulated_length - stroke_length);
        float height = L / (1 + std::exp(-k * t));
        stroke.heights.push_back(height);
        if (i < currentStroke.size() - 1) { accumulated_length += (currentStroke[i] - currentStroke[i+1]).norm(); }
    }
    return stroke;
}

void ScribbleArea::addStrokeToData() {
    if (m_type == StrokeType::BOUNDARY) {
        auto boundary_stroke = createBoundaryStroke();
        data.boundary.push_back(std::move(boundary_stroke));
    } else if (m_type == StrokeType::CONTOUR) {
        auto boundary_stroke = createBoundaryStroke();
        data.contour.push_back(std::move(boundary_stroke));
    } else if (m_type == StrokeType::CONVEX_BEND) {
        data.convex.push_back(std::move(currentStroke));
    } else if (m_type == StrokeType::CONCAVE_BEND) {
        data.concave.push_back(std::move(currentStroke));
    }
    currentStroke = Stroke();
}

void ScribbleArea::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    QRect dirtyRect = event->rect();
    painter.drawImage(dirtyRect, image, dirtyRect);
}

void ScribbleArea::resizeEvent(QResizeEvent *event)
{
    if (width() > image.width() || height() > image.height()) {
        int newWidth = qMax(width() + 128, image.width());
        int newHeight = qMax(height() + 128, image.height());
        resizeImage(&image, QSize(newWidth, newHeight));
        update();
    }
    QWidget::resizeEvent(event);
}

void ScribbleArea::drawLineTo(const QPoint &endPoint)
{
    QPainter painter(&image);
    painter.setPen(QPen(myPenColor, myPenWidth, Qt::SolidLine, Qt::RoundCap,
                        Qt::RoundJoin));
    painter.drawLine(lastPoint, endPoint);
    modified = true;

    Eigen::Vector2f start = Eigen::Vector2f(lastPoint.x(), lastPoint.y());
    Eigen::Vector2f end = Eigen::Vector2f(endPoint.x(), endPoint.y());
    float length = (end - start).norm();
    if (length > 0) {
        float t = 0;
        while (t < 1) {
            t += DISTANCE_BETWEEN_POINTS / length;
            Eigen::Vector2f point = start + (end - start) * t;
            currentStroke.push_back(point / Sketch::SKETCH_SCALE);
        }
        currentStroke.push_back(end / Sketch::SKETCH_SCALE);
    }
    strokeCount = 0;
    strokeCount++;

    int rad = (myPenWidth / 2) + 2;
    update(QRect(lastPoint, endPoint).normalized()
                                     .adjusted(-rad, -rad, +rad, +rad));
    lastPoint = endPoint;
}

void ScribbleArea::resizeImage(QImage *image, const QSize &newSize)
{
    if (image->size() == newSize)
        return;

    data.width = width();
    data.height = height();
    QImage newImage(newSize, QImage::Format_RGB32);
    newImage.fill(qRgb(255, 255, 255));
    QPainter painter(&newImage);
    painter.drawImage(QPoint(0, 0), *image);
    *image = newImage;
}

void ScribbleArea::print()
{
#ifndef QT_NO_PRINTER

#endif // QT_NO_PRINTER
}
