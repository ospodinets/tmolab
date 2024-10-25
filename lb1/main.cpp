#include "stdafx.h"
#include "lb1.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    lb1 w;
    w.show();
    return a.exec();
}
