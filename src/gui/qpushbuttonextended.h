#ifndef QPUSHBUTTONEXTENDED_H
#define QPUSHBUTTONEXTENDED_H

#include <QPushButton>
#include <QWidget>
#include <QIcon>

class QPushButtonExtended : public QPushButton
{
    Q_OBJECT

public:

    explicit QPushButtonExtended(QWidget *parent = 0);

signals:

public slots:

    void changeIcon(const QIcon &anIcon);
};

#endif // QPUSHBUTTONEXTENDED_H
