#ifndef VIEWOPTIONS_H
#define VIEWOPTIONS_H

#include <QMetaType>

enum viewOptions
{
    typeOfGradient_uniform,
    typeOfGradient_horizontal,
    typeOfGradient_vertical
};

Q_DECLARE_METATYPE(viewOptions)

#endif // VIEWOPTIONS_H
