#ifndef GRID_H
#define GRID_H

#include <QWidget>

namespace Ui {
class Grid;
}

class Grid : public QWidget
{
    Q_OBJECT

public:
    explicit Grid(QWidget *parent = 0);
    ~Grid();

private:
    Ui::Grid *ui;
};

#endif // GRID_H
