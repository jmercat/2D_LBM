#include "LBMController.hpp"
#include <QApplication>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    LBMController controller;
    controller.operate();
    return a.exec();
}
