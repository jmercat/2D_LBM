#ifndef LBMCONTROLLER_HPP
#define LBMCONTROLLER_HPP


#include "lbm.hpp"

#include <QThread>


class LBMWrapper : public QObject
{
    Q_OBJECT
    LBM<gridSizeX,gridSizeY>* lbm;

public:
    LBMWrapper(Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& obstacles, Eigen::Array<Eigen::Array<unsigned char,3,1>,Eigen::Dynamic,Eigen::Dynamic>& colors)
    {
        lbm = new LBM<gridSizeX,gridSizeY>(obstacles,colors);

    }
public slots:
    void compute()
    {
        lbm->compute(iterPerCall);
        emit end();
    }
signals:
    void end();
};

class LBMController : public QObject
{
    Q_OBJECT
    QThread LBMThread;
    Grid* w;
public:
    LBMController()
    {
        w = new Grid(gridSizeX,gridSizeY);
        LBMWrapper* lbm = new LBMWrapper(w->getObstacles(),w->getColor());
        lbm->moveToThread(&LBMThread);
        connect(&LBMThread,&QThread::finished, lbm, &QObject::deleteLater);
        connect(this, &LBMController::operate, lbm, &LBMWrapper::compute);
        connect(lbm, &LBMWrapper::end, this, &LBMController::handleResults);
        LBMThread.start();
        w->show();
    }

    ~LBMController()
    {
        LBMThread.quit();
        LBMThread.wait();
    }

public slots:
    void handleResults()
    {
        w->updateColor();
        emit operate();
    }

signals:
    void operate();
};

#endif // LBMCONTROLLER_HPP
