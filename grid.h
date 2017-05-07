#ifndef GRID_H
#define GRID_H

#include <QWidget>
#include <memory>

#include "settings.hpp"

#include <eigen3/Eigen/Dense>


namespace Ui
{
    class Grid;
}

class Grid : public QWidget
{
    Q_OBJECT

public:
    explicit Grid(int m, int n, QWidget *parent = 0);
    ~Grid();
    Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& getObstacles();
    Eigen::Array<Eigen::Array<unsigned char,3,1>,Eigen::Dynamic,Eigen::Dynamic>& getColor();
public slots:
    void clearGrid();
    void updateColor();
signals:
    void compute(unsigned int nIter);
protected:
    const int windowWidth, windowHeight;
    int mWidth, mHeight;
    Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> mGrid;
    Eigen::Array<Eigen::Array<unsigned char,3,1>,Eigen::Dynamic,Eigen::Dynamic> mColorGrid;
    Eigen::Array<QRect,Eigen::Dynamic,Eigen::Dynamic> mGridRect;

    int mCursorI, mCursorJ;
    bool mIsLeftClick;
    bool mIsRightClick;
    bool mIsRunning;


    void setGridRect();


    void paintEvent(QPaintEvent *e);
    void mouseMoveEvent(QMouseEvent *mouseEvent);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void drawGrid(QPainter *painter);
    void keyPressEvent(QKeyEvent *event);
    void draw(int i, int j, int color);
private:
    Ui::Grid *ui;
};

#endif // GRID_H
