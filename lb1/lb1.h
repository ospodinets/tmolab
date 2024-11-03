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
        int ns;
        int t0;
        int t1;
        double lambda_s;
        double lambda_s1;
        double lambda_s2;
        double dist;

        LambdaS(int ns = 0, int t0 = 0, int t1 = 0, double lambda_s = 0.0, double lambda_s1 = 0.0, double lambda_s2 = 0.0, double dist = 0.0)
            : ns(ns), t0(t0), t1(t1), lambda_s(lambda_s), lambda_s1(lambda_s1), lambda_s2(lambda_s2), dist(dist)
        {
        }
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
    QVector<int> loadIntervals(const QString& path);

private slots:
    bool showIntervals(const QString& path);
    void browse();
    void browse1();
    void browse2();
    void evaluate();
    bool proc1();
    bool proc2();
    bool proc3();
    bool proc4();
    bool proc56();
    bool proc6();
    void evaluateCompareFlows();

private:
    Ui::lb1Class ui;

    QVector<int> m_intervals;
    QVector<LambdaS> m_lambdaS;
    QVector<ApproxFunction> m_approxFunctions;
    bool m_flowIsTrivial;
    bool m_proc5isOk;
    QVector<LambdaS> m_lambdaSDiscrete;

    QVector<int> m_intervals1;
    bool m_flow1IsTrivial;
    QVector<int> m_intervals2;
    bool m_flow2IsTrivial;
};
