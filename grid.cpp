#include "grid.h"
#include "ui_grid.h"
#include <QGenericMatrix>
#include <QtWidgets>
#include <QtGui>
#include <QtCore>

using namespace boost::numeric::ublas;
int Grid::windowWidth = 600;
int Grid::windowHeight = 600;

Grid::Grid(int m, int n, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Grid)
{
    mWidth = m;
    mHeight = n;
    mGrid = std::shared_ptr<matrix<int> >(new matrix<int>(n+2,m+2));
    mGridRect = std::unique_ptr<matrix<QRect> >(new matrix<QRect>(n,m));
    this->setFixedSize(600,600);
    mCursorI = 0;
    mCursorJ = 0;
    mIsLeftClick = false;
    mIsRightClick = false;
    mIsRunning = false;
    setGridRect();
    this->setMouseTracking(true);
    ui->setupUi(this);
}

Grid::~Grid()
{
    delete ui;
}

std::shared_ptr<boost::numeric::ublas::matrix<int> > Grid::getObstacles()
{
    return mGrid;
}

void Grid::clearGrid()
{
    mGrid->clear();
    this->update();
}

void Grid::setGridRect()
{
    int rectWidth = windowWidth/mWidth;
    int rectHeight = windowHeight/mHeight;

    mGrid->clear();
    for(int i = 0; i<mWidth; i++)
    {
        for (int j= 0; j<mHeight; j++)
        {
            (*mGrid)(i+1,j+1) = 0;
            (*mGridRect)(i,j).setCoords(i*rectWidth+1,j*rectHeight+1,(i+1)*rectWidth-1,(j+1)*rectHeight-1);
        }
    }
    this->update();
}

void Grid::drawGrid(QPainter* painter)
{
    QPen linePen(Qt::gray);
    linePen.setWidth(2);
    QBrush brushWhite(Qt::white,Qt::SolidPattern);
    QBrush brushGray(Qt::gray,Qt::SolidPattern);
    QBrush brushBlack(Qt::black,Qt::SolidPattern);

    for(int i = 0; i<mWidth; i++)
    {
        for (int j= 0; j<mHeight; j++)
        {
            painter->setPen(linePen);
            painter->drawRect((*mGridRect)(i,j));
            if ((*mGrid)(i+1,j+1)==1) // divide by 2
            {
                painter->fillRect((*mGridRect)(i,j),brushBlack);
            }
            else
                painter->fillRect((*mGridRect)(i,j),brushWhite);
        }
    }
    painter->fillRect((*mGridRect)(mCursorI,mCursorJ),brushGray);
}

void Grid::draw(int i, int j, int color)
{
    (*mGrid)(i+1,j+1)=color;
}

void Grid::paintEvent(QPaintEvent *e)
{
    QPainter painter(this);
    drawGrid(&painter);
}

void Grid::mouseMoveEvent(QMouseEvent *mouseEvent)
{
    int i,j;

    int rectWidth = windowWidth/mWidth;
    int rectHeight = windowHeight/mHeight;

      i = qMin(qMax(mouseEvent->pos().x()/rectWidth,0),mWidth-1);
      j = qMin(qMax(mouseEvent->pos().y()/rectHeight,0),mHeight-1);
      if (mCursorI != i || mCursorJ != j)
      {
//          drawCursor(mGrid,mCursorI,mCursorJ,-1);
//          drawCursor(mGrid,i,j,1);

          if (mIsLeftClick)
              this->draw(i,j,1);
          else if (mIsRightClick)
              this->draw(i,j,0);

          mCursorI = i;
          mCursorJ = j;
      }
      this->update();
}

void Grid::mousePressEvent(QMouseEvent *event)
{
    int i,j;

    int rectWidth = windowWidth/mWidth;
    int rectHeight = windowHeight/mHeight;

    i = qMin(qMax(event->pos().x()/rectWidth,0),mWidth-1);
    j = qMin(qMax(event->pos().y()/rectHeight,0),mHeight-1);
    if (event->button()==Qt::LeftButton)
    {
        draw(i,j,1);
        mIsLeftClick = 1;
    }else if (event->button()==Qt::RightButton)
    {
        draw(i,j,0);
        mIsRightClick = 1;
    }
}


    void Grid::mouseReleaseEvent(QMouseEvent *event)
    {
        mIsLeftClick = 0;
        mIsRightClick = 0;
    }


    void Grid::keyPressEvent(QKeyEvent *event)
    {
        if (event->key()==Qt::Key_Enter || event->key()==Qt::Key_Return)
        {   
            std::cout << "Enter key event" << std::endl;
            emit compute(1000);

        }else if(event->key()==Qt::Key_0)
        {
            mGrid->clear();
            this->update();
        }
    }

