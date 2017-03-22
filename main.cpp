#include "lbm.hpp"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    constexpr unsigned int x = 100, y = 200;
    Grid w(x,y);
//    std::shared_ptr<boost::numeric::ublas::matrix<int> > grid();
    LBM<x,y> lbm(w.getObstacles(),w.getColor());
    QObject::connect(&w,SIGNAL(compute(unsigned int)),&lbm,SLOT(compute(unsigned int)));
    QObject::connect(&lbm,SIGNAL(colorUpdated()),&w,SLOT(updateColor()));
    w.show();

    return a.exec();
}
