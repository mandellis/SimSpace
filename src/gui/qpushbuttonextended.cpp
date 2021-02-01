#include "qpushbuttonextended.h"

QPushButtonExtended::QPushButtonExtended(QWidget *parent) : QPushButton(parent)
{
    ;
}

void QPushButtonExtended::changeIcon(const QIcon &anIcon)
{
    this->setIcon(anIcon);
}
