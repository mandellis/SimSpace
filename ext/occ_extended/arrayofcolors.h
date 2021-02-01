#ifndef ArrayOfColors_H
#define ArrayOfColors_H

#include <Quantity_Color.hxx>
#include <Quantity_Array1OfColor.hxx>
#include <NCollection_Array1.hxx>

class ArrayOfColors
{
public:

    //! constructor
    ArrayOfColors();

    //! destructor
    ~ArrayOfColors();

    //! get a color
    Quantity_Color getColor(int index);

private:

    //! array of colors
    NCollection_Array1<Quantity_Color> *myArrayOfColors;
};

#endif // ArrayOfColors_H
