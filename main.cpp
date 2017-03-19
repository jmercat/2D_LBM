#include "grid.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Grid w;
    w.show();

    return a.exec();
}
