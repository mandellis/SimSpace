#ifndef COLORSELECTOR_H
#define COLORSELECTOR_H

#include <QMetaType>
#include <QComboBox>

class ColorSelector : public QComboBox
{
public:

    ColorSelector(const ColorSelector&);
    ColorSelector(QWidget *parent=0);
    ~ColorSelector() { ; }

public:

    QColor color() const;
    void setColor(QColor c);

private:

    void populateList();
};

Q_DECLARE_METATYPE(ColorSelector)

#endif // COLORSELECTOR_H
