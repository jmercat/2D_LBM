#include "grid.h"
#include "ui_grid.h"
#include <QGenericMatrix>
#include <QtWidgets>
#include <QtGui>
#include <QtCore>

#include <iostream>

Grid::Grid(int m, int n, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Grid),
    windowWidth(n*600/m),
    windowHeight(600)
{
    mWidth = n;
    mHeight = m;
    mGrid.resize(n+2,m+2);
    mColorGrid.resize(n+2,m+2);
    mGridRect.resize(n,m);
    this->setFixedSize(windowWidth,windowHeight);
    mCursorI = 0;
    mCursorJ = 0;
    mIsLeftClick = false;
    mIsRightClick = false;
    mIsRunning = false;
    setGridRect();
    this->setMouseTracking(true);
    ui->setupUi(this);
//    for(int i = 0; i<mWidth; i++)
//    {
//        for (int j= 0; j<mHeight; j++)
//        {
//            mGrid(i+1,j+1) = 1;
//        }
//    }
    mGrid.setOnes();
    for(int i = 0; i<mWidth+2; i++)
    {
        for (int j= 0; j<mHeight+2; j++)
        {
            mColorGrid(i,j)(0)=0;
            mColorGrid(i,j)(1)=0;
            mColorGrid(i,j)(2)=0;
        }
    }
    this->update();
}

Grid::~Grid()
{
    delete ui;
}

Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> &Grid::getObstacles()
{
    return mGrid;
}

Eigen::Array<Eigen::Array<unsigned char,3,1>, Eigen::Dynamic, Eigen::Dynamic> &Grid::getColor()
{
    return mColorGrid;
}

void Grid::clearGrid()
{
    mGrid.setOnes();
    mColorGrid.setOnes();
    for(int i = 0; i<mWidth; i++)
    {
        for (int j= 0; j<mHeight; j++)
        {
            mGrid(i+1,j+1) = 1;
        }
    }
    this->update();
}

void Grid::updateColor()
{
    this->repaint();
}

void Grid::setGridRect()
{
    int rectWidth = windowWidth/mWidth;
    int rectHeight = windowHeight/mHeight;

    mGrid.setZero();
    for(int i = 0; i<mWidth; i++)
    {
        for (int j= 0; j<mHeight; j++)
        {
            mGrid(i+1,j+1) = 0;
            mGridRect(i,j).setCoords(i*rectWidth,j*rectHeight,(i+1)*rectWidth,(j+1)*rectHeight);
//            mGridRect(i,j).setCoords(i*rectWidth+1,j*rectHeight+1,(i+1)*rectWidth-1,(j+1)*rectHeight-1); // show grid
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
//            painter->drawRect(mGridRect(i,j));
            if (mGrid(i+1,j+1)==1)
            {
                QColor color(mColorGrid(i,j)(0),mColorGrid(i,j)(1),mColorGrid(i,j)(2));
                QBrush brushColor(color,Qt::SolidPattern);
                painter->fillRect(mGridRect(i,j),brushColor);
            }
            else
            {
                painter->fillRect(mGridRect(i,j),brushWhite);
            }
        }
    }
    painter->fillRect(mGridRect(mCursorI,mCursorJ),brushGray);
}

void Grid::draw(int i, int j, int color)
{
    mGrid(i+1,j+1)=color;
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
          {
              this->draw(i,j,0);
          }
          else if (mIsRightClick)
          {
              this->draw(i,j,1);
          }
          this->update(mGridRect(mCursorI,mCursorJ));
          this->update(mGridRect(i,j));
          mCursorI = i;
          mCursorJ = j;
      }
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
        draw(i,j,0);
        mIsLeftClick = 1;
    }else if (event->button()==Qt::RightButton)
    {
        draw(i,j,1);
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
            emit compute(10000);

        }else if(event->key()==Qt::Key_0)
        {
            mGrid.setOnes();
            this->update();
        }
    }

