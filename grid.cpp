#include "grid.h"
#include "ui_grid.h"

Grid::Grid(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Grid)
{
    ui->setupUi(this);
}

Grid::~Grid()
{
    delete ui;
}
