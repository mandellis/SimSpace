#include "ArrayOfColors.h"
#include <Quantity_NameOfColor.hxx>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
ArrayOfColors::ArrayOfColors()
{
    myArrayOfColors = new NCollection_Array1<Quantity_Color>(0,16);

    myArrayOfColors->SetValue(0,Quantity_NOC_AZURE);
    myArrayOfColors->SetValue(1,Quantity_NOC_AZURE2);
    myArrayOfColors->SetValue(2,Quantity_NOC_AZURE3);
    myArrayOfColors->SetValue(3,Quantity_NOC_AZURE4);
    myArrayOfColors->SetValue(4,Quantity_NOC_BISQUE);
    myArrayOfColors->SetValue(5,Quantity_NOC_BISQUE2);
    myArrayOfColors->SetValue(6,Quantity_NOC_BISQUE3);
    myArrayOfColors->SetValue(7,Quantity_NOC_BISQUE4);

    myArrayOfColors->SetValue(8,Quantity_NOC_BURLYWOOD);
    myArrayOfColors->SetValue(9,Quantity_NOC_BURLYWOOD1);
    myArrayOfColors->SetValue(10,Quantity_NOC_BURLYWOOD2);
    myArrayOfColors->SetValue(11,Quantity_NOC_BURLYWOOD3);
    myArrayOfColors->SetValue(12,Quantity_NOC_BURLYWOOD4);

    myArrayOfColors->SetValue(13,Quantity_NOC_LEMONCHIFFON1);
    myArrayOfColors->SetValue(14,Quantity_NOC_LEMONCHIFFON2);
    myArrayOfColors->SetValue(15,Quantity_NOC_LEMONCHIFFON3);
    myArrayOfColors->SetValue(16,Quantity_NOC_LEMONCHIFFON4);
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
ArrayOfColors::~ArrayOfColors()
{
    delete(myArrayOfColors);
}

//! --------------------------
//! function: returns a color
//! details:
//! --------------------------
Quantity_Color ArrayOfColors::getColor(int index)
{
    return myArrayOfColors->Value(index%(myArrayOfColors->Size()));
}
