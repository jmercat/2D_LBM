#ifndef LBMCONTROLLER_HPP
#define LBMCONTROLLER_HPP


#include "lbm.hpp"

#include <QThread>
#include <QKeyEvent>
#include <QLayout>
#include <QPushButton>
#include <QCoreApplication>


class LBMWrapper : public QObject
{
    Q_OBJECT
    LBM<gridSizeX,gridSizeY>* lbm;

public:
    LBMWrapper(Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic>& obstacles, Eigen::Array<Eigen::Array<float,3,1>,Eigen::Dynamic,Eigen::Dynamic>& results)
    {
        lbm = new LBM<gridSizeX,gridSizeY>(obstacles,results);

    }
public slots:
    void compute()
    {
        lbm->compute(iterPerCall);
        emit end();
    }
    void start(bool isOn)
    {
        if(isOn)
        {
            lbm->compute(iterPerCall);
            emit end();
        }
    }

signals:
    void end();
};

class LBMController : public QObject
{
    Q_OBJECT
    QThread LBMThread;
    Grid* w;
    QWidget* mainWidget;
    QHBoxLayout* HLayout;
    QVBoxLayout* VLayout;
    QPushButton* startButton;
    QPushButton* showSpeedButton;
public:
    LBMController()
    {
        mainWidget = new QWidget;
        w = new Grid(gridSizeX,gridSizeY,mainWidget);
        HLayout = new QHBoxLayout;
        VLayout = new QVBoxLayout;
        startButton = new QPushButton(mainWidget);
        startButton->setCheckable(true);

        showSpeedButton = new QPushButton(mainWidget);
        showSpeedButton->setCheckable(true);
        this->showSpeedButtonName(true);
        LBMWrapper* lbm = new LBMWrapper(w->getObstacles(),w->getResults());
        lbm->moveToThread(&LBMThread);

        connect(&LBMThread,&QThread::finished, lbm, &QObject::deleteLater);
        connect(this, &LBMController::operate, lbm, &LBMWrapper::compute);
        connect(lbm, &LBMWrapper::end, this, &LBMController::handleResults);
        LBMThread.start();
        HLayout->addWidget(w);
        VLayout->addWidget(startButton);
        VLayout->addWidget(showSpeedButton);
        HLayout->addLayout(VLayout);
        mainWidget->setLayout(HLayout);
        connect(startButton, &QPushButton::clicked, lbm, &LBMWrapper::start);
        connect(startButton, &QPushButton::clicked, this, &LBMController::startButtonName);
//        connect(showSpeedButton, &QPushButton::clicked, lbm, &LBMWrapper::start);
        connect(showSpeedButton, &QPushButton::clicked, this, &LBMController::showSpeedButtonName);
        mainWidget->show();
        startButton->click();
        showSpeedButton->click();
    }

    void keyPressEvent(QKeyEvent *event)
    {

        QCoreApplication::postEvent(w,event);
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
        if(startButton->isChecked())
            emit operate();
    }
    void startButtonName(bool isOn)
    {
        if(isOn)
            startButton->setText("Pause");
        else
            startButton->setText("Continue");
    }
    void showSpeedButtonName(bool isOn)
    {
        w->setShowSpeed(isOn);
        if(isOn)
            showSpeedButton->setText("Show pressure");
        else
            showSpeedButton->setText("Show speed");
    }

signals:
    void operate();
};

#endif // LBMCONTROLLER_HPP
