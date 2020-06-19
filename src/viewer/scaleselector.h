#ifndef SCALESELECTOR_H
#define SCALESELECTOR_H

#include <QWidget>
#include <QComboBox>

class QDoubleValidator;

class ScaleSelector : public QComboBox
{
    Q_OBJECT

public:

    ScaleSelector(QWidget *parent = 0);
    ~ScaleSelector();

    double getScale() { return this->currentData().toDouble(); }

private:

    QDoubleValidator *doubleValidator;

private slots:

    void emitScaleChanged();

signals:

    void scaleChanged(double scale);
};

#endif // SCALESELECTOR_H
