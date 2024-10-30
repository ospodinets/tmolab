#ifndef CUSTOMTICKER_H
#define CUSTOMTICKER_H

#include "qcustomplot.h"
#include <QVector>

class CustomTicker : public QCPAxisTicker {
public:
    explicit CustomTicker(const QVector<double>& customValues, int precision = 2);

    void setCustomValues(const QVector<double>& values);

protected:
    virtual QVector<double> createTickVector(double tickStep, const QCPRange& range);
    // Overrides the default label generator
    virtual QString getTickLabel(double tick, const QLocale &locale, QChar formatChar, int precision) override;

private:
    QVector<double> customValues;
    int precision;
};

#endif // CUSTOMTICKER_H
