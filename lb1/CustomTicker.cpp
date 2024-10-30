#include "stdafx.h"
#include "CustomTicker.h"

CustomTicker::CustomTicker(const QVector<double>& customValues, int precision)
    : QCPAxisTicker(), customValues(customValues), precision(precision) {}

void CustomTicker::setCustomValues(const QVector<double>& values) {
    customValues = values;
}

QVector<double> CustomTicker::createTickVector(double tickStep, const QCPRange& range)
{
    return customValues;
}

QString CustomTicker::getTickLabel(double tick, const QLocale& locale, QChar formatChar, int precision) {
    return QString::number(tick, 'f', this->precision);
}
