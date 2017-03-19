#include "lbm.hpp"
#include <QApplication>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Grid w(10,10);
//    std::shared_ptr<boost::numeric::ublas::matrix<int> > grid();
    LBM<10,10> lbm(w.getObstacles());
    QObject::connect(&w,SIGNAL(compute(unsigned int)),&lbm,SLOT(compute(unsigned int)));
    w.show();

    return a.exec();
}
