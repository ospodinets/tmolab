#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_lb1.h"

#include <QVector>

class lb1 : public QMainWindow
{
    Q_OBJECT

public:
    lb1(QWidget *parent = nullptr);
    ~lb1();

private:
    void disableUI();

private slots:
    bool loadIntervas(const QString& path);
    void browse();
    void evaluate();
    bool proc1();
    bool proc2();
    bool proc3();
    bool proc4();

private:
    Ui::lb1Class ui;

    QVector<int> m_intervals;
    bool m_flowIsTrivial;
};
