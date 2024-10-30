#include "stdafx.h"
#include "MyTableItemDelegate.h"

CustomDelegate::CustomDelegate(Qt::Alignment alignment, int precision, QObject* parent)
    : QStyledItemDelegate(parent), alignment(alignment), precision(precision) {}

void CustomDelegate::setAlignment(Qt::Alignment align) {
    alignment = align;
}

void CustomDelegate::setPrecision(int prec) {
    precision = prec;
}

void CustomDelegate::paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index) const {
    painter->save();

    // Retrieve and format the data
    QVariant value = index.data(Qt::DisplayRole);
    QString text;
    if (value.typeId() == QMetaType::Double) {
        text = QString::number(value.toDouble(), 'f', precision);
        // if there more zeroes - enable scientific notation
        if (text.count("0") > precision - 1 ) {
            text = QString::number(value.toDouble(), 'e', precision);
        }
    }
    else {
        text = value.toString();
    }

    // Set up the alignment and draw the text
    painter->drawText(option.rect, alignment, text);
    painter->restore();
}

void CustomDelegate::setEditorData(QWidget* editor, const QModelIndex& index) const {
    if (QDoubleSpinBox* spinBox = qobject_cast<QDoubleSpinBox*>(editor)) {
        double value = index.data(Qt::EditRole).toDouble();
        spinBox->setValue(value);
    }
    else {
        QStyledItemDelegate::setEditorData(editor, index);
    }
}

QWidget* CustomDelegate::createEditor(QWidget* parent, const QStyleOptionViewItem&, const QModelIndex& index) const {
    return nullptr;
    if (index.data().typeId() == QMetaType::Double) {
        QDoubleSpinBox* spinBox = new QDoubleSpinBox(parent);
        spinBox->setDecimals(precision);
        spinBox->setRange(-1e6, 1e6); // Set range as needed
        return spinBox;
    }
    return QStyledItemDelegate::createEditor(parent, QStyleOptionViewItem(), index);
}

void CustomDelegate::setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const {
    if (QDoubleSpinBox* spinBox = qobject_cast<QDoubleSpinBox*>(editor)) {
        spinBox->interpretText();
        model->setData(index, spinBox->value(), Qt::EditRole);
    }
    else {
        QStyledItemDelegate::setModelData(editor, model, index);
    }
}
