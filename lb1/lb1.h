#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_lb1.h"

class lb1 : public QMainWindow
{
    Q_OBJECT

public:
    lb1(QWidget *parent = nullptr);
    ~lb1();

private:
    Ui::lb1Class ui;
};
