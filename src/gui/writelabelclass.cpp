#include "writelabelclass.h"

#include <Graphic3d_TransformPers.hxx>
#include <TCollection_ExtendedString.hxx>
#include <iostream>

using namespace std;

writeLabelClass::~writeLabelClass()
{
    //!cout<<"writeLabelClass()->____destructor called____"<<endl;
}

writeLabelClass::writeLabelClass(QObject *parent):QObject(parent)
{
    myText = new AIS_TextLabel();
}

writeLabelClass::writeLabelClass(const occHandle(AIS_InteractiveContext) &aCTX, QObject *parent):QObject(parent)
{
    myText = new AIS_TextLabel();
    myCTX = aCTX;
}

void writeLabelClass::setContext(const occHandle(AIS_InteractiveContext) &aCTX)
{
    if(!aCTX.IsNull()) myCTX = aCTX;
}


void writeLabelClass::setText(const QString &aText)
{
    cout<<"writeLabelClass::setText()->____function called____"<<endl;
    occHandle(Graphic3d_TransformPers) trs = new Graphic3d_TransformPers(Graphic3d_TMF_2d);
    trs->SetCorner2d(Aspect_TOTP_RIGHT_UPPER);
    Graphic3d_Vec2i offset(150,20);
    trs->SetOffset2d(offset);
    myText->SetZLayer(Graphic3d_ZLayerId_TopOSD);
    myText->SetTransformPersistence(trs);

    //! set the color
    myText->SetColor(Quantity_NOC_BLACK);
    myText->SetFontAspect(Font_FA_Regular);
    myText->SetFont("Arial");
    myText->SetHeight(14);

    //! the text content
    Standard_CString theCString = aText.toUtf8();
    TCollection_ExtendedString theExtString(theCString);
    myText->SetText(theExtString);

    if(!myCTX.IsNull())
    {
        cout<<"writeLabelClass::setText()->____displaying text____"<<endl;
        myCTX->Display(myText,true);
    }
    else
    {
        cout<<"writeLabelClass::setText()->____cannot display text: context NULL____"<<endl;
    }
}
