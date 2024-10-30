#ifndef CUSTOMDELEGATE_H
#define CUSTOMDELEGATE_H

#include <QStyledItemDelegate>
#include <QDoubleSpinBox>
#include <QPainter>

class CustomDelegate : public QStyledItemDelegate
{
    Q_OBJECT

public:
    explicit CustomDelegate(Qt::Alignment alignment = Qt::AlignLeft, int precision = 2, QObject* parent = nullptr);

    void setAlignment(Qt::Alignment align);
    void setPrecision(int prec);

    void paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index) const override;
    void setEditorData(QWidget* editor, const QModelIndex& index) const override;
    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const override;
    void setModelData(QWidget* editor, QAbstractItemModel* model, const QModelIndex& index) const override;

private:
    Qt::Alignment alignment;
    int precision;
};

#endif // CUSTOMDELEGATE_H
