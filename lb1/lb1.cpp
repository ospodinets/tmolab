#include "stdafx.h"
#include "lb1.h"


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

    double gamma_lanczos(double x) {
        static const double coefficients[] = {
            676.5203681218851,
           -1259.1392167224028,
            771.32342877765313,
           -176.61502916214059,
            12.507343278686905,
           -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };
        int g = 7;  // Константа для зміщення x

        if (x < 0.5) {
            // Використовуємо відношення гамма-функції для від’ємних x
            return M_PI / (sin(M_PI * x) * gamma_lanczos(1 - x));
        }

        x -= 1;
        double a = 0.99999999999980993;  // Початкове значення для сумування
        for (int i = 0; i < sizeof(coefficients) / sizeof(double); ++i) {
            a += coefficients[i] / (x + i + 1);
        }

        double t = x + g + 0.5;
        return sqrt(2 * M_PI) * pow(t, x + 0.5) * exp(-t) * a;
    }

    // Функція гамма (наближення за Стірлінгом)
    double gamma_stirling(double x) {
        return sqrt(2 * M_PI / x) * pow(x / M_E, x);
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
        return sum * step / (pow(2, k / 2.0) * gamma_lanczos(k / 2.0));
    }

    // Функція для обчислення квантиля Хі-квадрат методом бінарного пошуку
    double chisqr_q(int k, double p) {
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
    double integrate_simpson(const std::function<double(double)>& F, double a, double b, int n = 1000) {
        if (n % 2 != 0) ++n;  // Переконуємося, що n є парним
        double h = (b - a) / n;
        double integral = F(a) + F(b);

        for (int i = 1; i < n; i += 2) {
            integral += 4 * F(a + i * h);  // Множимо на 4 для непарних індексів
        }
        for (int i = 2; i < n - 1; i += 2) {
            integral += 2 * F(a + i * h);  // Множимо на 2 для парних індексів
        }

        integral *= h / 3;
        return integral;
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
{
    ui.setupUi(this);
    ui.inputTab->layout()->setAlignment(Qt::AlignTop);
    ui.flowParameterTab->layout()->setAlignment(Qt::AlignTop);
    ui.check3Group->layout()->setAlignment(Qt::AlignTop);
    ui.intensityParameterTab->layout()->setAlignment(Qt::AlignTop);
    ui.check4Group->layout()->setAlignment(Qt::AlignTop);
    ui.approximationTab->layout()->setAlignment(Qt::AlignTop);
    ui.check5Group->layout()->setAlignment(Qt::AlignTop);
    connect(ui.browseBtn, SIGNAL(clicked()), this, SLOT(browse()));
    connect(ui.pathEdit, SIGNAL(textChanged(const QString&)), this, SLOT(loadIntervas(const QString&)));
    connect(ui.mInputBox, SIGNAL(valueChanged(int)), this, SLOT(proc3()));
    connect(ui.mInputBox, SIGNAL(valueChanged(int)), this, SLOT(proc4()));
    connect(ui.mInputBox, SIGNAL(valueChanged(int)), this, SLOT(proc5()));
    connect(ui.approxFunctionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(proc5()));

    initApproximationFunctions();

    disableUI();

    ui.pathEdit->setText("C:/Users/ospodynets/Documents/EDU/TMO/opts/opt5.txt");
    loadIntervas(ui.pathEdit->text());
    evaluate();
}

lb1::~lb1()
{}

void lb1::disableUI()
{
    ui.tabWidget->setTabEnabled(1, false);
    ui.tabWidget->setTabEnabled(2, false);
    ui.tabWidget->setTabEnabled(3, false);
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

bool lb1::loadIntervas(const QString & path)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
    {
        QMessageBox::warning(this, "Warning", "Cannot open file");
        return false;
    }

    QTextStream in(&file);
    // load parse numbers separated by spaces
    QString line = in.readAll();
    QStringList numbers = line.split(" ");
    m_intervals.clear();
    for (const QString& number : numbers)
    {
        m_intervals.push_back(number.toInt());
    }
    file.close();

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
    headerLabels << "№" << "Інтервали між вимогами";
    ui.intervalsWidget->setHorizontalHeaderLabels(headerLabels);

    // update alignment 
    ui.intervalsWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignCenter);
    // set minimuw width for columns to fit header text
    ui.intervalsWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    // set minimum width for columns
    ui.intervalsWidget->horizontalHeader()->setMinimumSectionSize(100);

    // change comparator for items
    ui.intervalsWidget->setSortingEnabled(true);
    ui.intervalsWidget->sortByColumn(1, Qt::AscendingOrder);

    // graph
    ui.inputGraphWidget->clearGraphs();
    ui.inputGraphWidget->show();
    ui.inputGraphWidget->addGraph();
    ui.inputGraphWidget->graph(0)->setData(indices, values);
    ui.inputGraphWidget->xAxis->setLabel("номер");
    ui.inputGraphWidget->yAxis->setLabel("інтервал");
    ui.inputGraphWidget->xAxis->setRange(0, m_intervals.size());
    ui.inputGraphWidget->yAxis->setRange(0, max);

    ui.mInputBox->blockSignals(true);
    ui.mInputBox->setMinimum( 2 );
    ui.mInputBox->setMaximum(sqrt(m_intervals.size()) + 10);
    ui.mInputBox->setValue(sqrt(m_intervals.size()));
    ui.mInputBox->blockSignals(false);

    ui.inputGraphWidget->replot();

    return true;
}

void lb1::evaluate()
{
    m_flowIsTrivial = proc1();
    if (!m_flowIsTrivial)
        proc2();
    proc3();
    proc4();
    proc5();
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
        textLabel->position->setCoords( (x_lambda_s[0] + x_lambda_s[1]) / 2, ymax - 0.1 * yrange);
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

bool lb1::proc5()
{
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
        y_dist.push_back(1 - exp(-integrV));

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

    return true;
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
        if (loadIntervas(filename))
            evaluate();
    }
}

