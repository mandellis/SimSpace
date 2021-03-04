#ifndef SHAPECOMPARISON_H
#define SHAPECOMPARISON_H

#include <TopoDS_Shape.hxx>
#include <TopLoc_Datum3D.hxx>
#include "occhandle.h"
#include <Standard_NoSuchObject.hxx>
#include <iostream>
using namespace std;

class shapeComparison
{
public:

    shapeComparison(const TopoDS_Shape &aShape):myShape(aShape){;}
    void getLast(occHandle(TopLoc_Datum3D) &coordSys)
    {
        if(coordSys.IsNull()) coordSys = new TopLoc_Datum3D();
        cout<<"____get last____"<<endl;
        coordSys = myShape.Location().FirstDatum();
        cout<<"____"<<coordSys->Transformation().TranslationPart().X()<<", "<<coordSys->Transformation().TranslationPart().Y()
           <<", "<<coordSys->Transformation().TranslationPart().Z()<<"____"<<endl;
        try
        {
            cout<<"____tag02____"<<endl;
            myShape.Location().NextLocation();
            cout<<"____tag03____"<<endl;
        }
        catch(Standard_NoSuchObject)
        {
            return;
        }
        getLast(coordSys);
    }

private:

    TopoDS_Shape myShape;
};

#endif // SHAPECOMPARISON_H
