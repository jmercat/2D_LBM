#ifndef GRID_H
#define GRID_H

#include <QWidget>
#include <memory>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
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
    std::shared_ptr<boost::numeric::ublas::matrix<int> > getObstacles();
public slots:
    void clearGrid();
signals:
    void compute(unsigned int nIter);
protected:
    static int windowWidth, windowHeight;
    int mWidth, mHeight;
    std::shared_ptr<boost::numeric::ublas::matrix<int> > mGrid;
    std::unique_ptr<boost::numeric::ublas::matrix<QRect> > mGridRect;

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
