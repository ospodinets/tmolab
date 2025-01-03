﻿#include "stdafx.h"
#include "lb1.h"
#include "MyTableItemDelegate.h"
#include "CustomTicker.h"


#include <math.h>

namespace
{
    static const double s_alpha = 0.05;
    double inverseNormalCDF(double p) 
    {
        // Перевірка на допустимість значення p
        if (p <= 0.0 || p >= 1.0) {
            throw std::invalid_argument("p має бути в межах (0, 1)");
        }

        // Коефіцієнти для апроксимації
        static const double a1 = -39.69683028665376;
        static const double a2 = 220.9460984245205;
        static const double a3 = -275.9285104469687;
        static const double a4 = 138.3577518672690;
        static const double a5 = -30.66479806614716;
        static const double a6 = 2.506628277459239;

        static const double b1 = -54.47609879822406;
        static const double b2 = 161.5858368580409;
        static const double b3 = -155.6989798598866;
        static const double b4 = 66.80131188771972;
        static const double b5 = -13.28068155288572;

        static const double c1 = -0.007784894002430293;
        static const double c2 = -0.3223964580411365;
        static const double c3 = -2.400758277161838;
        static const double c4 = -2.549732539343734;
        static const double c5 = 4.374664141464968;
        static const double c6 = 2.938163982698783;

        static const double d1 = 0.007784695709041462;
        static const double d2 = 0.3224671290700398;
        static const double d3 = 2.445134137142996;
        static const double d4 = 3.754408661907416;

        // Обчислення квантиля за апроксимаційною формулою
        double q, r;
        if (p < 0.02425) {
            // Нижня частина інтервалу
            q = std::sqrt(-2 * std::log(p));
            return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
        }
        else if (p > 1 - 0.02425) {
            // Верхня частина інтервалу
            q = std::sqrt(-2 * std::log(1 - p));
            return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
        }
        else {
            // Центральна частина інтервалу
            q = p - 0.5;
            r = q * q;
            return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
                (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
        }
    }

    // Функція гамма (наближення за Стірлінгом)
    double gamma_stirling(double x) {
        return sqrt(2 * M_PI / x) * pow(x / M_E, x);
    }

    // Функція для обчислення неповної бета-функції
    double incompleteBeta(double x, double a, double b) {
        double sum = 0.0;
        double term = 1.0;
        for (int k = 0; k < 1000; ++k) {
            if (k > 0) {
                term *= (a + k - 1) * x / (k * (b - k));
            }
            sum += term;
            if (std::fabs(term) < 1e-10) break;  // Похибка
        }
        return sum * std::pow(x, a) * std::pow(1 - x, b) / (a * tgamma(a) * tgamma(b));
    }

    // Інтегральна функція (метод трапецій для наближеного обчислення)
    double chi_square_cdf(int k, double x) {
        int n = 10000;  // кількість підінтервалів для точності
        double step = x / n;
        double sum = 0.0;

        for (int i = 1; i < n; ++i) {
            double t = i * step;
            sum += pow(t, k / 2.0 - 1) * exp(-t / 2.0);
        }

        sum += (pow(0, k / 2.0 - 1) * exp(-0 / 2.0) + pow(x, k / 2.0 - 1) * exp(-x / 2.0)) / 2.0;  // Кінцеві точки
        return sum * step / (pow(2, k / 2.0) * tgamma(k / 2.0));
    }

    // Функція для обчислення квантиля Хі-квадрат методом бінарного пошуку
    double chisqr_q(int k, double p)
    {
        double low = 0.0;
        double high = 1000.0;
        double mid = 0.0;

        while (high - low > 1e-6) {  // Точність до 1e-6
            mid = (low + high) / 2.0;
            double cdf_value = chi_square_cdf(k, mid);

            if (cdf_value < p) {
                low = mid;
            }
            else {
                high = mid;
            }
        }

        return mid;
    }

    // Функція для обчислення визначеного інтегралу методом Сімпсона
    double integrate_simpson(const std::function<double(double)>& F, double a, double b, int n = 1000)
    {
        if (n % 2 != 0) ++n;  // Переконуємося, що n є парним
        double h = (b - a) / n;
        double integral = F(a) + F(b);

        for (int i = 1; i < n; i += 2)
        {
            integral += 4 * F(a + i * h);  // Множимо на 4 для непарних індексів
        }
        for (int i = 2; i < n - 1; i += 2)
        {
            integral += 2 * F(a + i * h);  // Множимо на 2 для парних індексів
        }

        integral *= h / 3;
        return integral;
    }

    // Функція для щільності розподілу Стьюдента
    double t_density(double x, int df)
    {
        double gamma_part = tgamma((df + 1) / 2.0) / (sqrt(df * M_PI) * tgamma(df / 2.0));
        return gamma_part * pow(1 + x * x / df, -(df + 1) / 2.0);
    }

    // Наближення інтегралом для обчислення CDF розподілу Стьюдента
    double t_cdf(double x, int df, int num_steps = 10000) {
        double step = x / num_steps;
        double sum = 0.0;
        for (int i = 0; i < num_steps; ++i) {
            double xi = i * step;
            sum += t_density(xi, df) * step;
        }
        return 0.5 + sum;  // додаємо 0.5 для симетрії (t для двостороннього тесту)
    }

    // Метод Ньютона для знаходження критичного значення
    double student_critical_value(int df, double alpha, double tol = 1e-6, int max_iter = 100)
    {
        double x = 1.0;  // початкове наближення
        for (int i = 0; i < max_iter; ++i) {
            double fx = t_cdf(x, df) - (1.0 - alpha / 2.0); // для двостороннього тесту
            double dfx = t_density(x, df);  // похідна за x

            if (fabs(fx) < tol) break;

            x -= fx / dfx;
        }
        return x;
    }

    bool isFlowTrivial( const QVector<int>& intervals_ )
    {
        QVector <int> intervals = intervals_;
        const double N = intervals.size();

        std::sort(intervals.begin(), intervals.end());

        QVector<int> d;
        for (int i = 1; i < intervals.size(); ++i)
        {
            d.push_back(intervals[i] - intervals[i - 1]);
        }

        QVector<int> dN;
        // значення нормованих величин
        for (int i = 0; i < d.size(); ++i)
        {
            dN.push_back((N - i + 1) * d[i]);
        }

        double V = 0;
        for (int i = 0; i < dN.size() - 1; ++i)
        {
            for (int j = i + 1; j < dN.size(); ++j)
            {
                if (dN[i] > dN[j])
                {
                    V++;
                }
                else if (dN[i] == dN[j])
                {
                    V += 0.5;
                }
            }
        }

        double EV = (N * (N - 1)) / 4;
        double sigmaV = sqrt((N * (N - 1) * (2 * N + 5)) / 72);
        double u = (V + 0.5 - EV) / sigmaV;


        // квантіль стандартного нормального розподілу
        double u_alpha = inverseNormalCDF(1 - s_alpha / 2);

        return abs(u) <= u_alpha;
    }

    // Кумулятивна функція розподілу Фішера
    double fisherCDF(double x, int df1, int df2) {
        double a = df1 / 2.0;
        double b = df2 / 2.0;
        double y = (df1 * x) / (df1 * x + df2);
        return incompleteBeta(y, a, b);
    }

    double fisher_critical_value(double alpha, int df1, int df2, double tol = 1e-6) {
        double left = 0.0;
        double right = 10.0;  // Початковий верхній межу можна збільшити
        double mid;

        // Бінарний пошук для знаходження квантиля
        while (right - left > tol) {
            mid = (left + right) / 2;
            if (fisherCDF(mid, df1, df2) < 1.0 - alpha) {
                left = mid;
            }
            else {
                right = mid;
            }
        }

        return mid;
    }

    
    void test_t_quantile()
    {
        QVector<double> p = { 0.05 };

        for (auto alpha : p)
        {
            for (int i = 1; i < 130; ++i)
            {
                double x = student_critical_value(i, alpha);
                printf("student_critical_value(%f, %d) = %f\n", alpha, i, x);
            }

        }

    }

    void test_chi_square_quantile()
    {
        QVector<double> p = { 0.025, 0.975 };

        for (auto alpha : p)
        {
            for( int i = 1; i < 10; ++i )
            { 
                double x = chisqr_q(i, alpha);
                printf("chi_square_quantile(%f, %d) = %f\n", alpha, i, x);
            }
            
        }
        
    }
}


lb1::lb1(QWidget *parent)
    : QMainWindow(parent)
    , m_flowIsTrivial(false)
    , m_proc5isOk(false)
{
    ui.setupUi(this);
    ui.inputTab->layout()->setAlignment(Qt::AlignTop);
    ui.intervalsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 2, ui.intervalsWidget));
    ui.intervalsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.intervalsWidget->setAlternatingRowColors(true);

    ui.check1ResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.check1ResultsWidget->setAlternatingRowColors(true);
    ui.check1ResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));

    ui.check2ResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.check2ResultsWidget->setAlternatingRowColors(true);
    ui.check2ResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));

    ui.flowParameterTab->layout()->setAlignment(Qt::AlignTop);
    ui.check3Group->layout()->setAlignment(Qt::AlignTop);
    ui.check3ResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.check3ResultsWidget->setAlternatingRowColors(true);
    ui.check3ResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));

    ui.intensityParameterTab->layout()->setAlignment(Qt::AlignTop);
    ui.check4Group->layout()->setAlignment(Qt::AlignTop);
    ui.check4ResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.check4ResultsWidget->setAlternatingRowColors(true);
    ui.check4ResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));

    ui.approximationTab->layout()->setAlignment(Qt::AlignTop);
    ui.check5Group->layout()->setAlignment(Qt::AlignTop);
    ui.check5ResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.check5ResultsWidget->setAlternatingRowColors(true);
    ui.check5ResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));

    ui.approximationTab2->layout()->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    ui.check6Group->layout()->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    ui.check6ResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.check6ResultsWidget->setAlternatingRowColors(true);
    ui.check6ResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));


    ui.flowComparisonTab->layout()->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    ui.intervals1Widget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 2, ui.intervalsWidget));
    ui.intervals1Widget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.intervals1Widget->setAlternatingRowColors(true);

    ui.intervals2Widget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 2, ui.intervalsWidget));
    ui.intervals2Widget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.intervals2Widget->setAlternatingRowColors(true);

    ui.checkCompResultsWidget->horizontalHeader()->setStyleSheet("QHeaderView::section { font-weight: bold; background-color: #d3d7cf; border: 1px; }"
        " QHeaderView::section:hover{background-color: #babdb6;} ");
    ui.checkCompResultsWidget->setAlternatingRowColors(true);
    ui.checkCompResultsWidget->setItemDelegate(new CustomDelegate(Qt::AlignCenter, 5, ui.intervalsWidget));

    connect(ui.browseBtn, SIGNAL(clicked()), this, SLOT(browse()));
    connect(ui.pathEdit, SIGNAL(textChanged(const QString&)), this, SLOT(showIntervals(const QString&)));
    connect(ui.mInputBox, SIGNAL(valueChanged(int)), this, SLOT(proc3()));
    connect(ui.mInputBox, SIGNAL(valueChanged(int)), this, SLOT(proc4()));
    connect(ui.mInputBox, SIGNAL(valueChanged(int)), this, SLOT(proc56()));
    connect(ui.approxFunctionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(proc56()));


    connect(ui.browse1Btn, SIGNAL(clicked()), this, SLOT(browse1()));
    connect(ui.browse2Btn, SIGNAL(clicked()), this, SLOT(browse2()));

    initApproximationFunctions();

    disableUI();
}

lb1::~lb1()
{}

void lb1::disableUI()
{
    ui.tabWidget->setTabEnabled(1, false);
    ui.tabWidget->setTabEnabled(2, false);
    ui.tabWidget->setTabEnabled(3, false);
    ui.tabWidget->setTabEnabled(4, false);

    ui.intervalsWidget->clear();
    ui.intervalsWidget->hide();
    ui.check1ResultsWidget->clear();
    ui.check1Group->hide();
    ui.check2ResultsWidget->clear();
    ui.check2Group->hide(); 
    ui.check3Group->hide();
    ui.check3ResultsWidget->clear();
    ui.check4Group->hide();
    ui.check4ResultsWidget->clear();
    ui.check5Group->hide();
    ui.check5ResultsWidget->clear();
    ui.check6Group->hide();
    ui.check6ResultsWidget->clear();

    ui.inputGraphWidget->clearGraphs();
    ui.inputGraphWidget->replot();
    ui.inputGraphWidget->hide();

    ui.flowParameterGraphWidget->clearGraphs();
    ui.flowParameterGraphWidget->replot();
    ui.flowParameterGraphWidget->hide();

    ui.intensityParameterGraphWidget->clearGraphs();
    ui.intensityParameterGraphWidget->replot();
    ui.intensityParameterGraphWidget->hide();

    ui.approximationGraphWidget->clearGraphs();
    ui.approximationGraphWidget->replot();
    ui.approximationGraphWidget->hide();

    ui.distributionGraphWidget->clearGraphs();
    ui.distributionGraphWidget->replot();
    ui.distributionGraphWidget->hide();

    ui.approximation2GraphWidget->clearGraphs();
    ui.approximation2GraphWidget->replot();
    ui.approximation2GraphWidget->hide();

    ui.distribution2GraphWidget->clearGraphs();
    ui.distribution2GraphWidget->replot();
    ui.distribution2GraphWidget->hide();

    ui.intervals1Widget->clear();
    ui.intervals1Widget->hide();
    ui.intervals2Widget->clear();
    ui.intervals2Widget->hide();
    ui.checkCompResultsWidget->clear();
    ui.checkCompGroup->hide();

    ui.flow1IsNotTrivial->setText("");
    ui.flow2IsNotTrivial->setText("");
}

void lb1::initApproximationFunctions()
{
    ApproxFunction f;
    f.f = [](double t, double a, double b, double = 0) { return a + b * t; };
    f.name = "λ(t) = a + bt";
    f.params = { 0.0, 0.0 };
    m_approxFunctions.push_back(f);

    f.f = [](double t, double a, double b, double c) { return a + b * t + c * t * t; };
    f.params = { 0.0, 0.0, 0.0 };
    f.name = "λ(t) = a + bt + ct^2";
    m_approxFunctions.push_back(f);

    f.f = [](double t, double a, double b, double = 0) { return a + pow(t, b); };
    f.params = { 0.0, 0.0001 };
    f.name = "λ(t) = a + t^b";
    m_approxFunctions.push_back(f);

    f.f = [](double t, double a, double b, double = 0) { return a + pow(b, t); };
    f.params = { 0.0, 0.0 };
    f.name = "λ(t) = a + b^t";
    m_approxFunctions.push_back(f);

    f.f = [](double t, double a, double b, double = 0) { return a + exp( a + b*t ); };
    f.params = { 0.0, 0.0 };
    f.name = "λ(t) = exp(a + bt)";
    m_approxFunctions.push_back(f);

    ui.approxFunctionBox->blockSignals(true);
    ui.approxFunctionBox->clear();
    for (const auto& f : m_approxFunctions)
    {
        ui.approxFunctionBox->addItem(f.name);
    }    
    ui.approxFunctionBox->setCurrentIndex(-1);
    ui.approxFunctionBox->blockSignals(false);
}

QVector<int> lb1::loadIntervals(const QString& path)
{
    QVector<int> intervals;
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
    {
        QMessageBox::warning(this, "Warning", "Cannot open file");
        return intervals;
    }

    QTextStream in(&file);
    // load parse numbers separated by spaces
    QString line = in.readAll();
    QStringList numbers = line.split(" ");

    for (const QString& number : numbers)
    {
        intervals.push_back(number.toInt());
    }
    file.close();
    return intervals;
}

bool lb1::showIntervals(const QString & path)
{
    m_intervals = loadIntervals(path);

    ui.intervalsWidget->clear();
    ui.intervalsWidget->show();

    ui.intervalsWidget->setRowCount(m_intervals.size());
    ui.intervalsWidget->setColumnCount(2);

    QVector<double> indices;
    QVector<double> values;


    int max = 0;

    for (int i = 0; i < m_intervals.size(); ++i)
    {
        QTableWidgetItem* indexItem = new QTableWidgetItem();
        indexItem->setData(Qt::DisplayRole, i);
        ui.intervalsWidget->setItem(i, 0, indexItem);

        QTableWidgetItem* intervalItem = new QTableWidgetItem();
        intervalItem->setData(Qt::DisplayRole, m_intervals[i]);
        ui.intervalsWidget->setItem(i, 1, intervalItem);

        indices.push_back(i);
        values.push_back(m_intervals[i]);

        if (m_intervals[i] > max)
        {
            max = m_intervals[i];
        }
    }

    // update header text
    QStringList headerLabels;
    headerLabels << "№" << "Інтервал";
    ui.intervalsWidget->setHorizontalHeaderLabels(headerLabels);


    // update alignment 
    ui.intervalsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.intervalsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.intervalsWidget->horizontalHeader()->setMinimumSectionSize(75);

    // change comparator for items
    ui.intervalsWidget->setSortingEnabled(true);

    // graph
    ui.inputGraphWidget->clearGraphs();
    ui.inputGraphWidget->show();
    ui.inputGraphWidget->addGraph();
    ui.inputGraphWidget->graph(0)->setData(indices, values);
    ui.inputGraphWidget->xAxis->setLabel("номер");
    ui.inputGraphWidget->yAxis->setLabel("інтервал");
    ui.inputGraphWidget->xAxis->setRange(0, m_intervals.size());
    ui.inputGraphWidget->yAxis->setRange(0, max);
    
    ui.inputGraphWidget->graph(0)->setLineStyle(QCPGraph::lsImpulse);


    ui.mInputBox->blockSignals(true);
    ui.mInputBox->setMinimum( 2 );
    ui.mInputBox->setMaximum(sqrt(m_intervals.size()) + 10);
    ui.mInputBox->setValue(sqrt(m_intervals.size()));
    ui.mInputBox->blockSignals(false);

    ui.inputGraphWidget->replot();

    evaluate();

    return true;
}

void lb1::evaluate()
{
    m_flowIsTrivial = proc1();
    proc2();
    proc3();
    proc4();
    proc56();
}

bool lb1::proc1()
{
    QVector <int> intervals = m_intervals;
    const double N = intervals.size();

    std::sort(intervals.begin(), intervals.end());

    QVector<int> d;
    for (int i = 1; i < intervals.size(); ++i)
    {
        d.push_back(intervals[i] - intervals[i - 1]);
    }

    QVector<int> dN;
    // значення нормованих величин
    for (int i = 0; i < d.size(); ++i)
    {
        dN.push_back( (N - i + 1) * d[i] );
    }

    double V = 0;
    for (int i = 0; i < dN.size() - 1; ++i)
    {
        for (int j = i + 1; j < dN.size(); ++j)
        {
            if (dN[i] > dN[j])
            {
                V++;
            }
            else if (dN[i] == dN[j])
            {
                V += 0.5;
            }
        }
    }

    double EV = (N * (N - 1)) / 4;
    double sigmaV = sqrt((N * (N - 1) * (2 * N + 5)) / 72);
    double u = (V + 0.5 - EV) / sigmaV;


    // квантіль стандартного нормального розподілу
    double u_alpha = inverseNormalCDF(1 - s_alpha / 2);


    ui.check1ResultsWidget->clear();
    ui.check1ResultsWidget->setRowCount(1);
    ui.check1ResultsWidget->setColumnCount(4);


    QStringList headerLabels;
    QTableWidgetItem* Vitem = new QTableWidgetItem();
    Vitem->setData(Qt::DisplayRole, V);
    ui.check1ResultsWidget->setItem(0, 0, Vitem);
    headerLabels << "V";

    QTableWidgetItem* uItem = new QTableWidgetItem();
    uItem->setData(Qt::DisplayRole, u);
    ui.check1ResultsWidget->setItem(0, 1, uItem);
    headerLabels << "u";

    QTableWidgetItem* u_alphaItem = new QTableWidgetItem();
    u_alphaItem->setData(Qt::DisplayRole, u_alpha);
    ui.check1ResultsWidget->setItem(0, 2, u_alphaItem);
    headerLabels << "u alpha/2";

    QTableWidgetItem* resultItem = new QTableWidgetItem();

    bool outcome = abs(u) <= u_alpha;

    if (outcome)
    {
        resultItem->setText("Гіпотеза про експоненціальний закон розподілу приймається");
    }
    else
    {
        resultItem->setText("Гіпотеза про експоненціальний закон розподілу відхиляється");
    }
    ui.check1ResultsWidget->setItem(0, 3, resultItem);
    headerLabels << "Висновок";

    ui.check1ResultsWidget->setHorizontalHeaderLabels(headerLabels);

    // update alignment 
    ui.check1ResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.check1ResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.check1ResultsWidget->horizontalHeader()->setMinimumSectionSize(100);

    ui.check1Group->show();

    return outcome;
}

bool lb1::proc2()
{
    if (m_flowIsTrivial)
        return false;

    double W = 0;
    for (int i = 0; i < m_intervals.size() - 1; ++i)
    {
        for (int j = i + 1; j < m_intervals.size(); ++j)
        {
            if (m_intervals[i] < m_intervals[j])
            {
                W++;
            }
            else if (m_intervals[i] == m_intervals[j])
            {
                W += 0.5;
            }
        }
    }

    double EW = (m_intervals.size() * (m_intervals.size() - 1)) / 4;
    double sigmaW = sqrt((m_intervals.size() * (m_intervals.size() - 1) * (2 * m_intervals.size() + 5)) / 72);
    // статистика Манна
    double u = (W + 0.5 - EW) / sigmaW;

    // рівень значущості

    // квантіль стандартного нормального розподілу
    double u_alpha = inverseNormalCDF(1 - s_alpha / 2);

    ui.check2ResultsWidget->show();

    ui.check2ResultsWidget->clear();
    ui.check2ResultsWidget->setRowCount(1);
    ui.check2ResultsWidget->setColumnCount(4);


    QStringList headerLabels;
    QTableWidgetItem* Vitem = new QTableWidgetItem();
    Vitem->setData(Qt::DisplayRole, W);
    ui.check2ResultsWidget->setItem(0, 0, Vitem);
    headerLabels << "W";

    QTableWidgetItem* uItem = new QTableWidgetItem();
    uItem->setData(Qt::DisplayRole, u);
    ui.check2ResultsWidget->setItem(0, 1, uItem);
    headerLabels << "u";

    QTableWidgetItem* u_alphaItem = new QTableWidgetItem();
    u_alphaItem->setData(Qt::DisplayRole, u_alpha);
    ui.check2ResultsWidget->setItem(0, 2, u_alphaItem);
    headerLabels << "u alpha/2";

    QTableWidgetItem* resultItem = new QTableWidgetItem();

    if (abs(u) <= u_alpha)
    {
        resultItem->setText("Істотної тенденції до зміни інтервалів немає");
    }
    else if (u > u_alpha)
    {
        resultItem->setText("Існує тендеція до збільшення інтервалів");
    }
    else
    {
        resultItem->setText("Існує тендеція до зменшення інтервалів");
    }

    
    ui.check2ResultsWidget->setItem(0, 3, resultItem);
    headerLabels << "Висновок";

    ui.check2ResultsWidget->setHorizontalHeaderLabels(headerLabels);

    // update alignment 
    ui.check2ResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.check2ResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.check2ResultsWidget->horizontalHeader()->setMinimumSectionSize(100);

    ui.check2Group->show();

    return true;
}

bool lb1::proc3()
{

    QVector<double> t;
    t.push_back(0.0);

    const int N = m_intervals.size();

    for (int i = 0; i < N; ++i)
    {
        t.push_back(t.back() + m_intervals[i]);
    }

    // можливість введеня кількості класівю користувачем
    int m = ui.mInputBox->value();
    double dt = t.back() / m;

    const double u_alpha = inverseNormalCDF(1 - s_alpha / 2);
    const double u_alpha_2 = u_alpha * u_alpha;

    struct ClassDefinition
    {
        double dt;
        int mi;
        int ns;
        double mju_s;
        double mju_s1;
        double mju_s2;
        double variation;
    };

    QVector<ClassDefinition> cds;
    cds.reserve(m);

    for (int i = 0; i < m; ++i)
    {
        ClassDefinition cd;
        cd.dt = dt;
        cd.mi = i;
        cd.ns = 0;
        int min = cd.mi * cd.dt;
        int max = (cd.mi + 1) * cd.dt;

        for (auto ti : t)
        {
            if (ti >= min && ti < max)
            {
                ++cd.ns;
            }
        }

        cd.mju_s = cd.ns / (N * cd.dt);
        

        double a = 1/(cd.ns * cd.dt);
        double b = u_alpha_2 / 2;
        double c = u_alpha * sqrt(cd.ns + u_alpha_2 / 4);

        cd.mju_s1 = cd.mju_s + a * ( b - c );
        cd.mju_s2 = cd.mju_s + a * (b + c);

        double D = (cd.mju_s2 - cd.mju_s1) / 2 ;
        cd.variation = pow( D / u_alpha_2, 2 );

        cds.push_back(cd);
    }

    ui.tabWidget->setTabEnabled(1, true);
    ui.check3Group->show();
    ui.check3ResultsWidget->clear();

    ui.check3ResultsWidget->setRowCount(m);
    ui.check3ResultsWidget->setColumnCount(4);

    QStringList headerLabels;
    headerLabels << "Нижня границя" << "Параметр потоку" << "Верхня границя" << "Дисперсія";
    ui.check3ResultsWidget->setHorizontalHeaderLabels(headerLabels);

    // graph
    QVector<double> x_mju_s, y_mju_s;
    QVector<double> y_mju_s1;
    QVector<double> y_mju_s2;

    double ymax = DBL_MIN;
    double ymin = DBL_MAX;

    for (const auto& cd : cds)
    {
        QTableWidgetItem* item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, cd.mju_s1 );
        ui.check3ResultsWidget->setItem(cd.mi, 0, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, cd.mju_s);
        ui.check3ResultsWidget->setItem(cd.mi, 1, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, cd.mju_s2);
        ui.check3ResultsWidget->setItem(cd.mi, 2, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, cd.variation);
        ui.check3ResultsWidget->setItem(cd.mi, 3, item);

        x_mju_s.push_back(cd.mi * cd.dt);
        y_mju_s.push_back(cd.mju_s);
        y_mju_s1.push_back(cd.mju_s1);
        y_mju_s2.push_back(cd.mju_s2);


        x_mju_s.push_back((cd.mi + 1) * cd.dt);
        y_mju_s.push_back(cd.mju_s);
        y_mju_s1.push_back(cd.mju_s1);
        y_mju_s2.push_back(cd.mju_s2);

        if (cd.mju_s > ymax)
            ymax = cd.mju_s;

        if (cd.mju_s1 > ymax)
            ymax = cd.mju_s1;

        if (cd.mju_s2 > ymax)
            ymax = cd.mju_s2;

        if (cd.mju_s < ymin)
            ymin = cd.mju_s;

        if (cd.mju_s1 < ymin)
            ymin = cd.mju_s1;

        if (cd.mju_s2 < ymin)
            ymin = cd.mju_s2;
    }

    // update alignment 
    ui.check3ResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.check3ResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.check3ResultsWidget->horizontalHeader()->setMinimumSectionSize(200);

    ui.flowParameterGraphWidget->clearItems();
    ui.flowParameterGraphWidget->clearGraphs();
    ui.flowParameterGraphWidget->show();
    ui.flowParameterGraphWidget->addGraph();
    ui.flowParameterGraphWidget->graph(0)->setData(x_mju_s, y_mju_s);
    ui.flowParameterGraphWidget->addGraph();
    ui.flowParameterGraphWidget->graph(1)->setData(x_mju_s, y_mju_s1);
    ui.flowParameterGraphWidget->graph(1)->setPen(QPen(Qt::red));
    ui.flowParameterGraphWidget->addGraph();
    ui.flowParameterGraphWidget->graph(2)->setData(x_mju_s, y_mju_s2);
    ui.flowParameterGraphWidget->graph(2)->setPen(QPen(Qt::red));

    ui.flowParameterGraphWidget->xAxis->setLabel("t");
    ui.flowParameterGraphWidget->yAxis->setLabel("mju");
    ui.flowParameterGraphWidget->xAxis->setRange(x_mju_s.front(), x_mju_s.back());
    double yrange = ymax - ymin;
    ui.flowParameterGraphWidget->yAxis->setRange(ymin - 0.1 * yrange, ymax + 0.1 * yrange);

    ui.flowParameterGraphWidget->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui.flowParameterGraphWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.flowParameterGraphWidget->graph(0)->setLineStyle(QCPGraph::lsStepLeft);

    ui.flowParameterGraphWidget->graph(1)->setLineStyle(QCPGraph::lsLine);
    ui.flowParameterGraphWidget->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.flowParameterGraphWidget->graph(1)->setLineStyle(QCPGraph::lsStepLeft);

    ui.flowParameterGraphWidget->graph(2)->setLineStyle(QCPGraph::lsLine);
    ui.flowParameterGraphWidget->graph(2)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.flowParameterGraphWidget->graph(2)->setLineStyle(QCPGraph::lsStepLeft);

    for (const auto& x : x_mju_s)
    {
        QCPItemStraightLine* infLine = new QCPItemStraightLine(ui.flowParameterGraphWidget);
        infLine->point1->setCoords(x, 0);  // location of point 1 in plot coordinate
        infLine->point2->setCoords(x, 1);  // location of point 2 in plot coordinate
        infLine->setPen(QPen(Qt::gray));
    }

    QSharedPointer <CustomTicker> ticker{ new CustomTicker(x_mju_s, 0) };
    ui.flowParameterGraphWidget->xAxis->setTicker(ticker);

    ui.flowParameterGraphWidget->replot();

    return true;
}

// визначається кусково-стала функція інтенсивності
// тільки якщо потік не найпростіший
bool lb1::proc4()
{
    m_lambdaS.clear();
    ui.tabWidget->setTabEnabled(2, true);
    // update alignment 
    ui.check4ResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.check4ResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.check4ResultsWidget->horizontalHeader()->setMinimumSectionSize(200);

    ui.check4Group->show();

    int tmax = INT_MIN;
    int tmin = INT_MAX;
    int T = 0;
    for (auto& t : m_intervals)
    {
        if (t > tmax)
            tmax = t;

        if (t < tmin)
            tmin = t;

        T += t;
    }

    if (m_flowIsTrivial)
    {
        ui.check4ResultsWidget->hide();

        LambdaS ls;
        ls.t0 = 0;
        ls.t1 = T;
        
        ls.lambda_s1 = ls.lambda_s2 = ls.lambda_s = (double)m_intervals.size() / (double)T;
        m_lambdaS.push_back(ls);
    }
    else
    {
        // можливість введеня кількості класівю користувачем
        int m = ui.mInputBox->value();
        double dt = (tmax - tmin) / m;


        // get rid of it
        struct ClassDefinition
        {
            double dt;
            int mi;
            int ns;
            double lambda_s;
            double lambda_s1;
            double lambda_s2;
            double variation;
        };

        QVector<ClassDefinition> cds;
        cds.reserve(m);

        const double u_alpha = inverseNormalCDF(1 - s_alpha / 2);
        double u_alpha_2 = u_alpha * u_alpha;

        int Ns = 0;

        for (int i = 0; i < m - 1; ++i)
        {
            ClassDefinition cd;
            cd.dt = dt;
            cd.mi = i;

            cd.ns = 0;
            int min = tmin + cd.mi * cd.dt;
            int max = tmin + (cd.mi + 1) * cd.dt;

            for (auto ti : m_intervals)
            {
                if (ti >= min && ti < max)
                {
                    ++cd.ns;
                }
            }

            Ns += cd.ns;

            if (cd.ns != 0)
            {
                cd.lambda_s = cd.ns / ((m_intervals.size() - Ns) * dt);
                double Chi2_0 = chisqr_q(2 * cd.ns, 1 - s_alpha / 2);
                double Chi2_1 = chisqr_q(2 * cd.ns, s_alpha / 2);

                cd.lambda_s1 = cd.lambda_s * Chi2_0 / (2 * cd.ns);
                cd.lambda_s2 = cd.lambda_s * Chi2_1 / (2 * cd.ns);

                double D = pow(cd.lambda_s2 - cd.lambda_s1, 2);
                cd.variation = D / (2 * u_alpha_2);
            }
            else
            {
                cd.lambda_s = 0.0;
                cd.lambda_s1 = 0.0;
                cd.lambda_s2 = 0.0;
            }
            cds.push_back(cd);


            LambdaS ls;
            ls.ns = cd.ns;
            ls.t0 = min;
            ls.t1 = max;
            ls.lambda_s = cd.lambda_s;
            ls.lambda_s1 = cd.lambda_s1;
            ls.lambda_s2 = cd.lambda_s2;
            m_lambdaS.push_back(ls);
        }

        ui.check4ResultsWidget->show();
        ui.check4ResultsWidget->clear();

        ui.check4ResultsWidget->setRowCount(m - 1);
        ui.check4ResultsWidget->setColumnCount(4);

        QStringList headerLabels;
        headerLabels << "Нижня границя" << "Інтенсивність потоку" << "Верхня границя" << "Дисперсія";
        ui.check4ResultsWidget->setHorizontalHeaderLabels(headerLabels);

        for (const auto& cd : cds)
        {
            QTableWidgetItem* item = new QTableWidgetItem();
            item->setData(Qt::DisplayRole, cd.lambda_s1);
            ui.check4ResultsWidget->setItem(cd.mi, 0, item);

            item = new QTableWidgetItem();
            item->setData(Qt::DisplayRole, cd.lambda_s);
            ui.check4ResultsWidget->setItem(cd.mi, 1, item);

            item = new QTableWidgetItem();
            item->setData(Qt::DisplayRole, cd.lambda_s2);
            ui.check4ResultsWidget->setItem(cd.mi, 2, item);

            item = new QTableWidgetItem();
            item->setData(Qt::DisplayRole, cd.variation);
            ui.check4ResultsWidget->setItem(cd.mi, 3, item);
        }
    }

    // graph
    double ymax = DBL_MIN;
    double ymin = DBL_MAX;
    QVector<double> x_lambda_s, y_lambda_s;
    QVector<double> y_lambda_s1;
    QVector<double> y_lambda_s2;

    for( auto& ls : m_lambdaS )
    {
        x_lambda_s.push_back(ls.t0);
        y_lambda_s.push_back(ls.lambda_s);
        y_lambda_s1.push_back(ls.lambda_s1);
        y_lambda_s2.push_back(ls.lambda_s2);


        x_lambda_s.push_back(ls.t1);
        y_lambda_s.push_back(ls.lambda_s);
        y_lambda_s1.push_back(ls.lambda_s1);
        y_lambda_s2.push_back(ls.lambda_s2);

        if (ls.lambda_s > ymax)
            ymax = ls.lambda_s;

        if (ls.lambda_s1 > ymax)
            ymax = ls.lambda_s1;

        if (ls.lambda_s2 > ymax)
            ymax = ls.lambda_s2;

        if (ls.lambda_s < ymin)
            ymin = ls.lambda_s;

        if (ls.lambda_s1 < ymin)
            ymin = ls.lambda_s1;

        if (ls.lambda_s2 < ymin)
            ymin = ls.lambda_s2;
    }

    double yrange = ymax - ymin;
    if( yrange < 1e-06 )
        yrange = 2 * ymax;

    ui.intensityParameterGraphWidget->clearItems();
    ui.intensityParameterGraphWidget->clearGraphs();
    ui.intensityParameterGraphWidget->show();
    ui.intensityParameterGraphWidget->addGraph();
    ui.intensityParameterGraphWidget->graph(0)->setData(x_lambda_s, y_lambda_s);

    if (!m_flowIsTrivial)
    {
        ui.intensityParameterGraphWidget->addGraph();
        ui.intensityParameterGraphWidget->graph(1)->setData(x_lambda_s, y_lambda_s1);
        ui.intensityParameterGraphWidget->graph(1)->setPen(QPen(Qt::red));
        ui.intensityParameterGraphWidget->addGraph();
        ui.intensityParameterGraphWidget->graph(2)->setData(x_lambda_s, y_lambda_s2);
        ui.intensityParameterGraphWidget->graph(2)->setPen(QPen(Qt::red));

        ui.intensityParameterGraphWidget->graph(1)->setLineStyle(QCPGraph::lsLine);
        ui.intensityParameterGraphWidget->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);
        ui.intensityParameterGraphWidget->graph(1)->setLineStyle(QCPGraph::lsStepLeft);

        ui.intensityParameterGraphWidget->graph(2)->setLineStyle(QCPGraph::lsLine);
        ui.intensityParameterGraphWidget->graph(2)->setScatterStyle(QCPScatterStyle::ssNone);
        ui.intensityParameterGraphWidget->graph(2)->setLineStyle(QCPGraph::lsStepLeft);

        for (const auto& x : x_lambda_s)
        {
            QCPItemStraightLine* infLine = new QCPItemStraightLine(ui.intensityParameterGraphWidget);
            infLine->point1->setCoords(x, 0);  // location of point 1 in plot coordinate
            infLine->point2->setCoords(x, 1);  // location of point 2 in plot coordinate
            infLine->setPen(QPen(Qt::gray));
        }

        QSharedPointer <CustomTicker> ticker{ new CustomTicker(x_lambda_s, 0) };
        ui.intensityParameterGraphWidget->xAxis->setTicker(ticker);

        ui.intensityParameterGraphWidget->graph(0)->setLineStyle(QCPGraph::lsLine);
        ui.intensityParameterGraphWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
        ui.intensityParameterGraphWidget->graph(0)->setLineStyle(QCPGraph::lsStepLeft);
    }
    else
    {
        // add label with lambda value
        QCPItemText* textLabel = new QCPItemText(ui.intensityParameterGraphWidget);
        textLabel->setPositionAlignment(Qt::AlignHCenter | Qt::AlignTop);
        textLabel->position->setType(QCPItemPosition::ptPlotCoords);
        textLabel->position->setCoords( (x_lambda_s[0] + x_lambda_s[1]) / 2, ymax);
        textLabel->setText("λ = " + QString::number(y_lambda_s[0]));
        textLabel->setFont(QFont(font().family(), 10));
        textLabel->setPen(QPen(Qt::black));
    }

    ui.intensityParameterGraphWidget->xAxis->setLabel("t");
    ui.intensityParameterGraphWidget->yAxis->setLabel("lambda");
    ui.intensityParameterGraphWidget->xAxis->setRange(x_lambda_s.first(), x_lambda_s.last());    
    ui.intensityParameterGraphWidget->yAxis->setRange(ymin - 0.1*yrange, ymax+ 0.1*yrange);


    

    ui.intensityParameterGraphWidget->replot();
    return true;
}

// аппроксимація функції інтенсивності
// заданої функцією
bool lb1::proc56()
{
    m_proc5isOk = false;
    if( m_flowIsTrivial )
        return false;

    ui.tabWidget->setTabEnabled(3, true);

    ui.check5Group->hide();
    qApp->processEvents();

    if( ui.approxFunctionBox->currentIndex() < 0 )
       return false;

    ApproxFunction af = m_approxFunctions[ui.approxFunctionBox->currentIndex()];

    auto S_2 = [&af, lambdaValues = m_lambdaS](double a, double b, double c) {
        double sum = 0.0;
        for (auto& ls : lambdaValues)
        {
            double approx_val = af.f((ls.t0 + ls.t1) / 2, a, b, c);
            double err = ls.lambda_s - approx_val;
            //printf("av = %f ls = %f, err = %f\n", approx_val, ls.lambda_s, err);
            sum += err * err;
        }
        double S2 = sum / (lambdaValues.size() - af.params.size());
        return S2;
    };

    // розв'язуємо систему диференційних рівнянь 
    // dS_2 / da = 0, dS_2 / db = 0, dS_2 / dc = 0
    // для знаходження коефіцієнтів a, b, c
    // застосовуємо метод градієнтного спуску

    // початкові значення параметрів
    QVector<double> p = af.params;

    int max_iterations = 200000;
    double epsilon = 1e-8;
    double learningRate = 1e-3;
    double prev_S2 = DBL_MAX;

    // show busy cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    for (int iteration = 0; iteration < max_iterations; ++iteration)
    {
        // обчислення градієнтів
        QVector<double> dS2(p.size(), 0.0 );
            // Обчислення dS2/da, dS2/db, dS2/dc для кожного t
        if (p.size() == 3)
        {
            dS2[0] = (S_2(p[0] + epsilon, p[1], p[2]) - S_2(p[0] - epsilon, p[1], p[2])) / (2 * epsilon);
            dS2[1] = (S_2(p[0], p[1] + epsilon, p[2]) - S_2(p[0], p[1] - epsilon, p[2])) / (2 * epsilon);
            dS2[2] = (S_2(p[0], p[1], p[2] + epsilon) - S_2(p[0], p[1], p[2] - epsilon)) / (2 * epsilon);
        }
        else if (p.size() == 2)
        {
            dS2[0] = (S_2(p[0] + epsilon, p[1], 0) - S_2(p[0] - epsilon, p[1], 0)) / (2 * epsilon);
            dS2[1] = (S_2(p[0], p[1] + epsilon, 0) - S_2(p[0], p[1] - epsilon, 0)) / (2 * epsilon);
        }
        else if( p.size() == 1 )
        {
            dS2[0] = (S_2(p[0] + epsilon, 0, 0) - S_2(p[0] - epsilon, 0, 0)) / (2 * epsilon);
        }
        else
        {
            throw std::runtime_error("Invalid number of parameters");
        }

        QVector<double> pPrev = p;

        for (int i = 0; i < p.size(); ++i)
        {
            p[i] -= learningRate * dS2[i];
        }

        double current_S2 = S_2(p[0], p.size() > 1 ? p[1] : 0, p.size() > 2 ? p[2] : 0);
        if (current_S2 < prev_S2)
        {
            learningRate *= 1.1;
            prev_S2 = current_S2;
        }
        else
        {
            p = pPrev;
            current_S2 = prev_S2;
            learningRate *= 0.1;
        }


        // Умова завершення
        bool complete = true;
        for (int i = 0; i < dS2.size(); ++i)
        {
            if (std::abs(dS2[i]) > epsilon)
            {
                complete = false;
                break;
            }
        }

        if (complete)
            break;

    //    printf("#%d S2(%f, %f, %f) = %f, lr = %f\n", iteration, p[0], p.size() > 1 ? p[1] : 0, p.size() > 2 ? p[2] : 0,
    //        current_S2, learningRate );
    }

    // restore cursor
    QApplication::restoreOverrideCursor();
    
    // значення аппроксимуючої функції 
    // разом з кусково-стаалою функцією лямбда у довірчих інтервалах

    QVector<double> x, y;
    QVector<double> x_lambda_s, y_lambda_s;
    QVector<double> y_lambda_s1;
    QVector<double> y_lambda_s2;
    double ymin = DBL_MAX;
    double ymax = DBL_MIN;
    for (auto& ls : m_lambdaS)
    {
        x.push_back((ls.t0 + ls.t1) / 2);
        y.push_back(af.f(x.back(), p[0], p.size() > 1 ? p[1] : 0, p.size() > 2 ? p[2] : 0));

        x_lambda_s.push_back(ls.t0);
        y_lambda_s.push_back(ls.lambda_s);
        y_lambda_s1.push_back(ls.lambda_s1);
        y_lambda_s2.push_back(ls.lambda_s2);


        x_lambda_s.push_back(ls.t1);
        y_lambda_s.push_back(ls.lambda_s);
        y_lambda_s1.push_back(ls.lambda_s1);
        y_lambda_s2.push_back(ls.lambda_s2);

        if (ls.lambda_s > ymax)
            ymax = ls.lambda_s;

        if (ls.lambda_s1 > ymax)
            ymax = ls.lambda_s1;

        if (ls.lambda_s2 > ymax)
            ymax = ls.lambda_s2;

        if (ls.lambda_s < ymin)
            ymin = ls.lambda_s;

        if (ls.lambda_s1 < ymin)
            ymin = ls.lambda_s1;

        if (ls.lambda_s2 < ymin)
            ymin = ls.lambda_s2;
    }

    // заносимо коефіціенти у таблицю
    ui.check5Group->show();

    ui.check5ResultsWidget->clear();
    ui.check5ResultsWidget->setRowCount(1);
    ui.check5ResultsWidget->setColumnCount(af.params.size());

    QStringList headerLabels;
    for (int i = 0; i < af.params.size(); ++i)
    {
        QTableWidgetItem* item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, p[i]);
        ui.check5ResultsWidget->setItem(0, i, item);
        headerLabels << ( i == 0 ? "a" : i == 1 ? "b" : "c");
    }

    ui.check5ResultsWidget->setHorizontalHeaderLabels(headerLabels);
    ui.check5ResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.check5ResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.check5ResultsWidget->horizontalHeader()->setMinimumSectionSize(200);

    ui.approximationGraphWidget->clearItems();
    ui.approximationGraphWidget->clearGraphs();
    ui.approximationGraphWidget->show();
    ui.approximationGraphWidget->addGraph();
    ui.approximationGraphWidget->graph(0)->setData(x_lambda_s, y_lambda_s);

    ui.approximationGraphWidget->addGraph();
    ui.approximationGraphWidget->graph(1)->setData(x_lambda_s, y_lambda_s1);
    ui.approximationGraphWidget->graph(1)->setPen(QPen(Qt::red));
    ui.approximationGraphWidget->addGraph();
    ui.approximationGraphWidget->graph(2)->setData(x_lambda_s, y_lambda_s2);
    ui.approximationGraphWidget->graph(2)->setPen(QPen(Qt::red));

    ui.approximationGraphWidget->addGraph();
    ui.approximationGraphWidget->graph(3)->setData(x, y);

    // set bold style for approximation graph
    QPen pen(Qt::darkGreen);
    pen.setWidth(2);
    ui.approximationGraphWidget->graph(3)->setPen(pen);   


    ui.approximationGraphWidget->xAxis->setLabel("t");
    ui.approximationGraphWidget->yAxis->setLabel("lambda");
    ui.approximationGraphWidget->xAxis->setRange(x_lambda_s.first(), x_lambda_s.last());
    double yrange = ymax - ymin;
    ui.approximationGraphWidget->yAxis->setRange(ymin - 0.1 * yrange, ymax + 0.1 * yrange);

    ui.approximationGraphWidget->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui.approximationGraphWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.approximationGraphWidget->graph(0)->setLineStyle(QCPGraph::lsStepLeft);

    ui.approximationGraphWidget->graph(1)->setLineStyle(QCPGraph::lsLine);
    ui.approximationGraphWidget->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.approximationGraphWidget->graph(1)->setLineStyle(QCPGraph::lsStepLeft);

    ui.approximationGraphWidget->graph(2)->setLineStyle(QCPGraph::lsLine);
    ui.approximationGraphWidget->graph(2)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.approximationGraphWidget->graph(2)->setLineStyle(QCPGraph::lsStepLeft);

    for (const auto& x : x_lambda_s)
    {
        QCPItemStraightLine* infLine = new QCPItemStraightLine(ui.approximationGraphWidget);
        infLine->point1->setCoords(x, 0);  // location of point 1 in plot coordinate
        infLine->point2->setCoords(x, 1);  // location of point 2 in plot coordinate
        infLine->setPen(QPen(Qt::gray));
    }

    QSharedPointer <CustomTicker> ticker{ new CustomTicker(x_lambda_s, 0) };
    ui.approximationGraphWidget->xAxis->setTicker(ticker);

    ui.approximationGraphWidget->replot();

    auto lambda_t = [&af, p](double t) {
        return af.f(t, p[0], p.size() > 1 ? p[1] : 0, p.size() > 2 ? p[2] : 0);
        };

    // додаємо функцію розподілу
    QVector<double> y_dist;
    double y_dist_max = DBL_MIN;
    double y_dist_min = DBL_MAX;
    for (auto& ls : m_lambdaS)
    {
        double t = (ls.t0 + ls.t1) / 2;
        double integrV = integrate_simpson(lambda_t, 0, t);
        ls.dist = 1 - exp(-integrV);
        y_dist.push_back( ls.dist );

        if (y_dist.back() > y_dist_max)
            y_dist_max = y_dist.back();

        if (y_dist.back() < y_dist_min)
            y_dist_min = y_dist.back();

        
    }

    ui.distributionGraphWidget->clearItems();
    ui.distributionGraphWidget->clearGraphs();
    ui.distributionGraphWidget->show();
    ui.distributionGraphWidget->addGraph();
    ui.distributionGraphWidget->graph(0)->setData(x, y_dist);

    ui.distributionGraphWidget->xAxis->setLabel("t");
    ui.distributionGraphWidget->yAxis->setLabel("F(t)");
    ui.distributionGraphWidget->xAxis->setRange(x.first(), x.last());
    yrange = y_dist_max - y_dist_min;
    ui.distributionGraphWidget->yAxis->setRange(y_dist_min - 0.1 * yrange, y_dist_max + 0.1 * yrange);

    ui.distributionGraphWidget->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui.distributionGraphWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);

    // set bold style for approximation graph
    QPen pen3(Qt::red);
    pen3.setWidth(2);
    ui.distributionGraphWidget->graph(0)->setPen(pen3);

    ui.distributionGraphWidget->replot();

    m_proc5isOk = true;

    proc6();

    return true;
}

// визначення вірогідної кусково-сталої функції інтенсивності
bool lb1::proc6()
{
    ui.tabWidget->setTabEnabled(4, false);

    if (m_flowIsTrivial || !m_proc5isOk)
        return false;

    ui.tabWidget->setTabEnabled(4, true);

    double t_student = student_critical_value(m_lambdaS.size() - 1, s_alpha / 2 );

    QVector<LambdaS> remaining = m_lambdaS;
    m_lambdaSDiscrete.clear();
    m_lambdaSDiscrete.push_back(remaining.front());
    remaining.pop_front();

    m_lambdaSDiscrete.back().lambda_s1 = m_lambdaSDiscrete.back().lambda_s2 = m_lambdaSDiscrete.back().lambda_s;

    while (!remaining.empty())
    {
        // H0: lambda_s = lambda_s1
        double l0 = m_lambdaSDiscrete.back().lambda_s;
        double l1 = remaining.front().lambda_s;

        double n0 = m_lambdaSDiscrete.back().ns;
        double n1 = remaining.front().ns;

        double A = (l1 - l0) / sqrt((n0 - 1) * pow(l1, 2) + (n1 - 1) * pow(l0, 2));
        double B = sqrt((n0 * n1 * (n0 + n1 - 2)) / (n0 + n1));

        double t = A * B;

        if (abs(t) <= t_student)
        {
            // H0 правильна
            auto& last = m_lambdaSDiscrete.back();
            last.ns += remaining.front().ns;
            last.t1 = remaining.front().t1;
            last.lambda_s = (n0 * l0 + n1 * l1) / (n0 + n1);
        }
        else
        {
            // H0 відхиляється
            m_lambdaSDiscrete.push_back(remaining.front());
            m_lambdaSDiscrete.back().lambda_s1 = m_lambdaSDiscrete.back().lambda_s2 = m_lambdaSDiscrete.back().lambda_s;
        }
        remaining.pop_front();
    }

    // заносимо значення у таблицю
    ui.check6Group->show();

    ui.check6ResultsWidget->clear();
    ui.check6ResultsWidget->setRowCount(m_lambdaSDiscrete.size());
    ui.check6ResultsWidget->setColumnCount(2);

    QStringList headerLabels;
    headerLabels << "t0..t1" << "Інтенсивність потоку";
    ui.check6ResultsWidget->setHorizontalHeaderLabels(headerLabels);

    for (int i = 0; i < m_lambdaSDiscrete.size(); ++i)
    {
        QTableWidgetItem* item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, QString::number(m_lambdaSDiscrete[i].t0) + ".." + QString::number(m_lambdaSDiscrete[i].t1));
        ui.check6ResultsWidget->setItem(i, 0, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, m_lambdaSDiscrete[i].lambda_s);
        ui.check6ResultsWidget->setItem(i, 1, item);
    }

    ui.check6ResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.check6ResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.check6ResultsWidget->horizontalHeader()->setMinimumSectionSize(50);

    // графік кусково-сталої функції з попередніх процедур
    
    QVector<double> x_lambda_s, y_lambda_s;
    QVector<double> y_lambda_s1;
    QVector<double> y_lambda_s2;
    double ymin = DBL_MAX;
    double ymax = DBL_MIN;
    for (auto& ls : m_lambdaS)
    {
        x_lambda_s.push_back(ls.t0);
        y_lambda_s.push_back(ls.lambda_s);
        y_lambda_s1.push_back(ls.lambda_s1);
        y_lambda_s2.push_back(ls.lambda_s2);


        x_lambda_s.push_back(ls.t1);
        y_lambda_s.push_back(ls.lambda_s);
        y_lambda_s1.push_back(ls.lambda_s1);
        y_lambda_s2.push_back(ls.lambda_s2);

        if (ls.lambda_s > ymax)
            ymax = ls.lambda_s;

        if (ls.lambda_s1 > ymax)
            ymax = ls.lambda_s1;

        if (ls.lambda_s2 > ymax)
            ymax = ls.lambda_s2;

        if (ls.lambda_s < ymin)
            ymin = ls.lambda_s;

        if (ls.lambda_s1 < ymin)
            ymin = ls.lambda_s1;

        if (ls.lambda_s2 < ymin)
            ymin = ls.lambda_s2;
    }

    QVector<double> x_lambda, y_lambda;
    for (auto& ls : m_lambdaSDiscrete)
    {
        x_lambda.push_back(ls.t0);
        y_lambda.push_back(ls.lambda_s);

        x_lambda.push_back(ls.t1);
        y_lambda.push_back(ls.lambda_s);

        if (ls.lambda_s > ymax)
            ymax = ls.lambda_s;

        if (ls.lambda_s1 > ymax)
            ymax = ls.lambda_s1;

        if (ls.lambda_s2 > ymax)
            ymax = ls.lambda_s2;

        if (ls.lambda_s < ymin)
            ymin = ls.lambda_s;

        if (ls.lambda_s1 < ymin)
            ymin = ls.lambda_s1;

        if (ls.lambda_s2 < ymin)
            ymin = ls.lambda_s2;
    }

    ui.approximation2GraphWidget->clearItems();
    ui.approximation2GraphWidget->clearGraphs();
    ui.approximation2GraphWidget->show();
    ui.approximation2GraphWidget->addGraph();
    ui.approximation2GraphWidget->graph(0)->setData(x_lambda_s, y_lambda_s);

    ui.approximation2GraphWidget->addGraph();
    ui.approximation2GraphWidget->graph(1)->setData(x_lambda_s, y_lambda_s1);
    ui.approximation2GraphWidget->graph(1)->setPen(QPen(Qt::red));
    ui.approximation2GraphWidget->addGraph();
    ui.approximation2GraphWidget->graph(2)->setData(x_lambda_s, y_lambda_s2);
    ui.approximation2GraphWidget->graph(2)->setPen(QPen(Qt::red));

    ui.approximation2GraphWidget->xAxis->setLabel("t");
    ui.approximation2GraphWidget->yAxis->setLabel("lambda");
    ui.approximation2GraphWidget->xAxis->setRange(x_lambda_s.first(), x_lambda_s.last());
    double yrange = ymax - ymin;
    ui.approximation2GraphWidget->yAxis->setRange(ymin - 0.1 * yrange, ymax + 0.1 * yrange);

    ui.approximation2GraphWidget->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui.approximation2GraphWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.approximation2GraphWidget->graph(0)->setLineStyle(QCPGraph::lsStepLeft);

    ui.approximation2GraphWidget->graph(1)->setLineStyle(QCPGraph::lsLine);
    ui.approximation2GraphWidget->graph(1)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.approximation2GraphWidget->graph(1)->setLineStyle(QCPGraph::lsStepLeft);

    ui.approximation2GraphWidget->graph(2)->setLineStyle(QCPGraph::lsLine);
    ui.approximation2GraphWidget->graph(2)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.approximation2GraphWidget->graph(2)->setLineStyle(QCPGraph::lsStepLeft);

    ui.approximation2GraphWidget->addGraph();
    ui.approximation2GraphWidget->graph(3)->setData(x_lambda, y_lambda);
    QPen pen(Qt::darkGreen);
    pen.setWidth(2);
    ui.approximation2GraphWidget->graph(3)->setPen(pen);

    for (const auto& x : x_lambda_s)
    {
        QCPItemStraightLine* infLine = new QCPItemStraightLine(ui.approximation2GraphWidget);
        infLine->point1->setCoords(x, 0);  // location of point 1 in plot coordinate
        infLine->point2->setCoords(x, 1);  // location of point 2 in plot coordinate
        infLine->setPen(QPen(Qt::gray));
    }

    QSharedPointer <CustomTicker> ticker{ new CustomTicker(x_lambda_s, 0) };
    ui.approximation2GraphWidget->xAxis->setTicker(ticker);

    ui.approximation2GraphWidget->replot();

    // функція сплайн-експоненційного розподілу

    double y_dist_max = DBL_MIN;
    double y_dist_min = DBL_MAX;

    QVector<double> x_prevDist, y_prevDist;
    for (auto& ls : m_lambdaS)
    {
        x_prevDist.push_back( (ls.t0 + ls.t1) / 2 );
        y_prevDist.push_back( ls.dist );

        if (y_prevDist.back() > y_dist_max)
            y_dist_max = y_prevDist.back();

        if (y_prevDist.back() < y_dist_min)
            y_dist_min = y_prevDist.back();
    }

    QVector<double> x_spline, y_spline;
    QVector<double> x_dots, y_dots;

    // для кожного класу створюємо сплайн
    for (int i = 0; i < m_lambdaSDiscrete.size(); ++i)
    {
        double t0 = m_lambdaSDiscrete[i].t0;
        double t1 = m_lambdaSDiscrete[i].t1;

        auto f = [i, lambdaS = m_lambdaSDiscrete](double t) ->double {

            double sum = 0.0;
            for (int s = 0; s < i; ++s)
            {
                sum += ( lambdaS[s].lambda_s - lambdaS[s+1].lambda_s ) * lambdaS[s].t1;
            }

            return 1 - exp( -lambdaS[i].lambda_s * t  - sum );
        };

        if (i < m_lambdaSDiscrete.size() - 1)
        {
            x_dots.push_back(t1);
            y_dots.push_back(f(t1));
        }
        
        const int Nsamples = 100;
        for (int p = 0; p < Nsamples; ++p )
        {
            double t = t0 + ( t1 - t0 ) * ( (double)p  / (double)Nsamples);
            x_spline.push_back(t);
            y_spline.push_back(f(t));

            if (y_spline.back() > y_dist_max)
                y_dist_max = y_spline.back();

            if (y_spline.back() < y_dist_min)
                y_dist_min = y_spline.back();
        }
    }

    ui.distribution2GraphWidget->clearItems();
    ui.distribution2GraphWidget->clearGraphs();
    ui.distribution2GraphWidget->show();
    ui.distribution2GraphWidget->addGraph();
    ui.distribution2GraphWidget->graph(0)->setData(x_prevDist, y_prevDist);

    ui.distribution2GraphWidget->xAxis->setLabel("t");
    ui.distribution2GraphWidget->yAxis->setLabel("F(t)");
    ui.distribution2GraphWidget->xAxis->setRange(x_spline.first(), x_spline.last());
    yrange = y_dist_max - y_dist_min;
    ui.distribution2GraphWidget->yAxis->setRange(y_dist_min - 0.1 * yrange, y_dist_max + 0.1 * yrange);

    ui.distribution2GraphWidget->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui.distribution2GraphWidget->graph(0)->setScatterStyle(QCPScatterStyle::ssNone);
    ui.distribution2GraphWidget->graph(0)->setPen(QPen(Qt::red));

    ui.distribution2GraphWidget->addGraph();
    QPen pen2(Qt::darkGreen);
    pen2.setWidth(2);
    ui.distribution2GraphWidget->graph(1)->setPen(pen2);
    ui.distribution2GraphWidget->graph(1)->setData(x_spline, y_spline);

    ui.distribution2GraphWidget->addGraph();
    ui.distribution2GraphWidget->graph(2)->setData(x_dots, y_dots);
    ui.distribution2GraphWidget->graph(2)->setLineStyle(QCPGraph::lsNone);
    ui.distribution2GraphWidget->graph(2)->setScatterStyle(QCPScatterStyle::ssDisc);
    ui.distribution2GraphWidget->graph(2)->setPen(pen2);


    for (const auto& x : x_dots)
    {
        QCPItemStraightLine* infLine = new QCPItemStraightLine(ui.distribution2GraphWidget);
        infLine->point1->setCoords(x, 0);  // location of point 1 in plot coordinate
        infLine->point2->setCoords(x, 1);  // location of point 2 in plot coordinate
        infLine->setPen(QPen(Qt::gray));
    }

    QSharedPointer<CustomTicker> ticker2{ new CustomTicker(x_dots, 0) };
    ui.distribution2GraphWidget->xAxis->setTicker(ticker);


    ui.distribution2GraphWidget->replot();


    return true;
}

void lb1::evaluateCompareFlows()
{
    ui.intervals1Widget->hide();
    ui.intervals2Widget->hide();

    ui.intervals1Widget->clear();
    ui.intervals1Widget->show();

    ui.intervals1Widget->setRowCount(m_intervals1.size());
    ui.intervals1Widget->setColumnCount(2);

    for (int i = 0; i < m_intervals1.size(); ++i)
    {
        QTableWidgetItem* indexItem = new QTableWidgetItem();
        indexItem->setData(Qt::DisplayRole, i);
        ui.intervals1Widget->setItem(i, 0, indexItem);

        QTableWidgetItem* intervalItem = new QTableWidgetItem();
        intervalItem->setData(Qt::DisplayRole, m_intervals1[i]);
        ui.intervals1Widget->setItem(i, 1, intervalItem);
    }

    // update header text
    QStringList headerLabels;
    headerLabels << "№" << "Інтервал";
    ui.intervals1Widget->setHorizontalHeaderLabels(headerLabels);


    // update alignment 
    ui.intervals1Widget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.intervals1Widget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.intervals1Widget->horizontalHeader()->setMinimumSectionSize(75);

    // change comparator for items
    ui.intervals1Widget->setSortingEnabled(true);
    ui.intervals1Widget->show();
    
    ui.intervals2Widget->clear();

    ui.intervals2Widget->setRowCount(m_intervals2.size());
    ui.intervals2Widget->setColumnCount(2);

    for (int i = 0; i < m_intervals2.size(); ++i)
    {
        QTableWidgetItem* indexItem = new QTableWidgetItem();
        indexItem->setData(Qt::DisplayRole, i);
        ui.intervals2Widget->setItem(i, 0, indexItem);

        QTableWidgetItem* intervalItem = new QTableWidgetItem();
        intervalItem->setData(Qt::DisplayRole, m_intervals2[i]);
        ui.intervals2Widget->setItem(i, 1, intervalItem);
    }

    // update header text
    ui.intervals2Widget->setHorizontalHeaderLabels(headerLabels);


    // update alignment 
    ui.intervals2Widget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.intervals2Widget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.intervals2Widget->horizontalHeader()->setMinimumSectionSize(75);

    // change comparator for items
    ui.intervals2Widget->setSortingEnabled(true);
    ui.intervals2Widget->show();

    ui.checkCompGroup->hide();

    double T1 = 0.0;
    double n1 = m_intervals1.size();
    for (int i = 0; i < m_intervals1.size(); ++i)
    {
        T1 += m_intervals1[i];
    }

    double T2 = 0.0;
    double n2 = m_intervals2.size();
    for (int i = 0; i < m_intervals2.size(); ++i)
    {
        T2 += m_intervals2[i];
    }

    bool outcome = false;

    double R = (n1 * T1) / (n2 * T2);
    if (n1 < 30 && n2 < 30)
    {
        double F = fisher_critical_value( 2 * n1, 2 * n2, s_alpha);
        outcome = R < F;

        ui.checkCompResultsWidget->clear();
        ui.checkCompResultsWidget->setRowCount(1);
        ui.checkCompResultsWidget->setColumnCount(3);

        QStringList headerLabels;
        headerLabels << "R" << "F" << "Висновок";
        ui.checkCompResultsWidget->setHorizontalHeaderLabels(headerLabels);

        QTableWidgetItem* item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, R);
        ui.checkCompResultsWidget->setItem(0, 0, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, F);
        ui.checkCompResultsWidget->setItem(0, 1, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, outcome ? "Потоки збігаються" : "Потоки не збігаються");
        ui.checkCompResultsWidget->setItem(0, 2, item);
        ui.checkCompGroup->show();
    }
    else
    {
        double z = log(R);
        double Ez = 0.5 * (1 / n1 - 1 / n2);
        double sigmaz = sqrt(1 / n1 + 1 / n2);

        double U = (z - Ez) / sigmaz;

        double u_alpha = inverseNormalCDF(1 - s_alpha / 2);
        outcome = abs(U) < u_alpha;

        ui.checkCompResultsWidget->clear();
        ui.checkCompResultsWidget->setRowCount(1);
        ui.checkCompResultsWidget->setColumnCount(3);

        QStringList headerLabels;
        headerLabels << "U" << "U_alpha" << "Висновок";
        ui.checkCompResultsWidget->setHorizontalHeaderLabels(headerLabels);

        QTableWidgetItem* item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, U);
        ui.checkCompResultsWidget->setItem(0, 0, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, u_alpha);
        ui.checkCompResultsWidget->setItem(0, 1, item);

        item = new QTableWidgetItem();
        item->setData(Qt::DisplayRole, outcome ? "Потоки збігаються" : "Потоки не збігаються");
        ui.checkCompResultsWidget->setItem(0, 2, item);
        ui.checkCompGroup->show();
    }   


    ui.checkCompResultsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.checkCompResultsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.checkCompResultsWidget->horizontalHeader()->setMinimumSectionSize(200);
}

void lb1::browse()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open File");
    if (!filename.isEmpty())
    {
        ui.pathEdit->blockSignals(true);
        ui.pathEdit->setText(filename);
        ui.pathEdit->blockSignals(false);
        disableUI();

        showIntervals( filename );
    }
}

void lb1::browse1()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open File");
    if (!filename.isEmpty())
    {
        ui.path1Edit->blockSignals(true);
        ui.path1Edit->setText(filename);
        ui.path1Edit->blockSignals(false);

        ui.flow1IsNotTrivial->setText("");

        ui.intervals1Widget->clear();
        ui.intervals1Widget->hide();
        ui.checkCompResultsWidget->clear();
        ui.checkCompGroup->hide();

        m_intervals1 = loadIntervals(filename);
        m_flow1IsTrivial = isFlowTrivial(m_intervals1);
        ui.flow1IsNotTrivial->setText(m_flow1IsTrivial ? "" : "Потік не є найпростішим!");

        if (!m_intervals1.isEmpty() && !m_intervals2.isEmpty() && m_flow1IsTrivial && m_flow2IsTrivial )
            evaluateCompareFlows();
    }
}

void lb1::browse2()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open File");
    if (!filename.isEmpty())
    {
        ui.path2Edit->blockSignals(true);
        ui.path2Edit->setText(filename);
        ui.path2Edit->blockSignals(false);

        ui.intervals2Widget->clear();
        ui.intervals2Widget->hide();
        ui.checkCompResultsWidget->clear();
        ui.checkCompGroup->hide();

        m_intervals2 = loadIntervals(filename);
        m_flow2IsTrivial = isFlowTrivial(m_intervals2);
        ui.flow2IsNotTrivial->setText(m_flow2IsTrivial ? "" : "Потік не є найпростішим!");
        if (!m_intervals1.isEmpty() && !m_intervals2.isEmpty() && m_flow1IsTrivial && m_flow2IsTrivial )
            evaluateCompareFlows();
    }
}

