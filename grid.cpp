#include "grid.h"
#include "ui_grid.h"

#include <math.h>
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
    mResultGrid.resize(n+2,m+2);
    mGridRect.resize(n,m);
    this->setFixedSize(windowWidth,windowHeight);
    mCursorI = 0;
    mCursorJ = 0;
    mIsLeftClick = false;
    mIsRightClick = false;
    mIsRunning = false;
    mShowSpeed = true;
    mBrush.resize(1,1);
    mBrush(0,0) = true;
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
    mRescaler(0) = 1;
    mRescaler(1) = 0;
    mPRescaler(0) = 1;
    mPRescaler(1) = 0;
    for(int i = 0; i<mWidth+2; i++)
    {
        for (int j= 0; j<mHeight+2; j++)
        {
            mResultGrid(i,j)(0)=0;
            mResultGrid(i,j)(1)=0;
            mResultGrid(i,j)(2)=0;
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

Eigen::Array<Eigen::Array<float,3,1>, Eigen::Dynamic, Eigen::Dynamic> &Grid::getResults()
{
    return mResultGrid;
}

void Grid::clearGrid()
{
    mGrid.setOnes();
    mResultGrid.setOnes();
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

void Grid::setShowSpeed(bool showSpeed)
{
    mShowSpeed = showSpeed;
    this->updateColor();
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

void Grid::setColor(QColor* colorToSet, float norm)
{
    float currentRescaler;
    if(mShowSpeed)
        currentRescaler = mRescaler(0);
    else
        currentRescaler = mPRescaler(0);
    if (norm*currentRescaler<0.5)
    {
        colorToSet->setRgb(255*currentRescaler*norm*2,255*currentRescaler*norm*2,255);
    }else if(norm*currentRescaler<1)
    {
        colorToSet->setRgb(255,255*(1-(norm*currentRescaler-0.5)*2),255*(1-(norm*currentRescaler-0.5)*2));
    }else
    {
        colorToSet->setRgb(200,0,0);
    }
}

void Grid::setBrush(const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& brush)
{
    mBrush = brush;
}

void Grid::drawGrid(QPainter* painter)
{
    QPen linePen(Qt::gray);
    linePen.setWidth(2);
    QBrush brushWhite(Qt::white,Qt::SolidPattern);
    QBrush brushGray(Qt::gray,Qt::SolidPattern);
    QBrush brushBlack(Qt::black,Qt::SolidPattern);
    QColor* color = new QColor;
    QColor* colorR = new QColor;

    float normMax = 0;
    float normMin = 1;
    float pMax = 0;
    float pMin = 1000;
    Eigen::Vector2f speed;
    float pressure;


    for(int i = 0; i<mWidth; i++)
    {
        for (int j= 0; j<mHeight; j++)
        {
            painter->setPen(linePen);
//            painter->drawRect(mGridRect(i,j));
            if (mGrid(i+1,j+1)==1)
            {
                speed(0) = mResultGrid(i,j)(0);
                speed(1) = mResultGrid(i,j)(1);
                pressure = mResultGrid(i,j)(2);
                float norm = speed.norm();
                QPointF p1, p2, pv;

                pv.setX(speed(0)/norm);
                pv.setY(speed(1)/norm);

                normMax = std::max(normMax,norm);
                normMin = std::min(normMin,norm);

                pMax = std::max(pMax,pressure);
                pMin = std::min(pMin,pressure);

                norm -= mRescaler(1);
                pressure -= mPRescaler(1);

                if(mShowSpeed)
                    this->setColor(color,norm);
                else
                    this->setColor(color,pressure);

                QBrush brushColor(*color,Qt::SolidPattern);
                painter->fillRect(mGridRect(i,j),brushColor);

                if(mShowSpeed)
                {
                    if(norm>1.e-5)
                    {
                        int r,g,b;
                        r = std::max(color->red()-directionContrast,0);
                        g = std::max(color->green()-directionContrast,0);
                        b = std::max(color->blue()-directionContrast,0);
                        QColor dirColor(r,g,b);
                        QPen dirPen(dirColor);
                        dirPen.setWidth(1);
                        painter->setPen(dirPen);
                        pv.setX(mResultGrid(i,j)(0)/norm);
                        pv.setY(mResultGrid(i,j)(1)/norm);

                        p1 = mGridRect(i,j).height()*pv/2+mGridRect(i,j).center();
                        p2 = -mGridRect(i,j).height()*pv/2+mGridRect(i,j).center();

                    }else
                    {
                        p1 = mGridRect(i,j).center();
                        p2 = mGridRect(i,j).center();
                    }
                    painter->drawLine(p1,p2);
                }
            }
            else
            {
                painter->fillRect(mGridRect(i,j),brushWhite);
            }
        }
    }
    mRescaler(1) = 0.2*normMin+0.8*mRescaler(1);
    mPRescaler(1) = 0.2*pMin+0.8*mPRescaler(1);
    if(normMax != normMin)
    {
        mRescaler(0) = 0.2/(normMax-normMin)+0.8*mRescaler(0);
    }
    if(pMax != pMin)
    {
        mPRescaler(0) = 0.2/(pMax-pMin)+0.8*mPRescaler(0);
    }
    painter->fillRect(mGridRect(mCursorI,mCursorJ),brushGray);
    delete color;
    delete colorR;
}

void Grid::drawPoint(int i, int j, int color)
{
    mGrid(i+1,j+1)=color;
}

void Grid::draw(int i, int j, int color)
{
    for (unsigned int jj = std::max(j-(mBrush.rows()+1)/2,long(0)), jjj = 0;
         jj<std::min(j+mBrush.rows()/2,long(gridSizeX)); jj++, jjj++)
    {
        for (unsigned int ii = std::max(i-(mBrush.cols()+1)/2,long(0)),
             iii = 0; ii<std::min(i+mBrush.cols()/2,long(gridSizeY)); ii++, iii++)
        {
            if(mBrush(iii,jjj))
                this->drawPoint(ii,jj,color);
        }
    }
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
            emit compute(iterPerCall);

        }else if(event->key()==Qt::Key_0)
        {
            mGrid.setOnes();
            this->update();
        }
    }

