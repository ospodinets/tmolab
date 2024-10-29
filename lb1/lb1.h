#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_lb1.h"

#include <QVector>
#include <functional>

class lb1 : public QMainWindow
{
    Q_OBJECT

private:
    struct LambdaS
    {
        int t0;
        int t1;
        double lambda_s;
        double lambda_s1;
        double lambda_s2;
    };  
    struct ApproxFunction
    {
        std::function<double(double, double, double, double)> f;
        QVector<double> params;
        QString name;
    };

public:
    lb1(QWidget *parent = nullptr);
    ~lb1();

private:
    void disableUI();
    void initApproximationFunctions();

private slots:
    bool loadIntervas(const QString& path);
    void browse();
    void evaluate();
    bool proc1();
    bool proc2();
    bool proc3();
    bool proc4();
    bool proc5();

private:
    Ui::lb1Class ui;

    QVector<int> m_intervals;
    QVector<LambdaS> m_lambdaS;
    QVector<ApproxFunction> m_approxFunctions;
    bool m_flowIsTrivial;
};
