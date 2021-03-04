//! ----------------
//! custom includes
//! ----------------
#include "occface.h"

//! ----
//! OCC
//! ----
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <TopAbs_State.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopLoc_Location.hxx>

//! ----
//! C++
//! ----
#include <algorithm>

//! ---------------------
//! function: deflection
//! details:
//! ---------------------
double OCCFace::deflection()
{
    return myMaxDeflection;
}

//! -------------------------------
//! function: computeMaxDefletcion
//! details:
//! -------------------------------
void OCCFace::computeDeflection()
{
    double maxDeflection = 0.0;
    std::vector<double> deflections;
    for(std::vector<TopoDS_Face>::const_iterator it = myFaces.cbegin(); it!=myFaces.cend(); it++)
    {
        const TopoDS_Face &curFace = *it;
        TopLoc_Location loc;
        double deflection = BRep_Tool::Triangulation(curFace,loc)->Deflection();
        if(deflection>=maxDeflection) maxDeflection = deflection;
        deflections.push_back(deflection);
    }
    myMaxDeflection = maxDeflection;
    cout<<"OCCFace::OCCFace()->____max deflection: "<<maxDeflection<<"____"<<endl;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
OCCFace::OCCFace(const TopoDS_Face &aFace)
{
    cout<<"OCCFace::OCCFace()->____constructor called____"<<endl;
    if(aFace.IsNull()) return;
    myFaces = std::vector<TopoDS_Face>();
    myFaces.push_back(aFace);
    this->computeDeflection();
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
OCCFace::OCCFace(const std::vector<TopoDS_Face> &listOfFaces)
{
    cout<<"OCCFace::OCCFace()->____constructor called____"<<endl;
    if(listOfFaces.size()==0) return;
    for(std::vector<TopoDS_Face>::const_iterator it=listOfFaces.cbegin(); it!=listOfFaces.cend(); ++it)
        myFaces.push_back(*it);
    this->computeDeflection();
}

//! ----------------------
//! function: setGeometry
//! details:
//! ----------------------
void OCCFace::setGeometry(const TopoDS_Face &aFace)
{
    myFaces.clear();
    myFaces.push_back(aFace);
    this->computeDeflection();
}

//! ----------------------
//! function: setGeometry
//! details:
//! ----------------------
void OCCFace::setGeometry(const std::vector<TopoDS_Face> &listOfFaces)
{
    myFaces = listOfFaces;
    this->computeDeflection();
}

//! ---------------
//! function: area
//! details:
//! ---------------
double OCCFace::area()
{
    double area = 0;
    for(std::vector<TopoDS_Face>::const_iterator it = myFaces.cbegin(); it!=myFaces.cend(); it++)
    {
        GProp_GProps prop;
        const TopoDS_Face &curFace = *it;
        BRepGProp::SurfaceProperties(curFace,prop);
        area = area + prop.Mass();
    }
    return area;
}

//! --------------------------
//! function: pointProjection
//! details:
//! --------------------------
bool OCCFace::pointProjection(double *aPoint, double *theProjectedPoint, double eps)
{
    //cout<<"OCCFace::pointProjection()->____function called for ("<<aPoint[0]<<", "<<aPoint[1]<<", "<<aPoint[2]<<")____"<<endl;
    std::vector<std::vector<double>> vecProjections;
    std::vector<double> vecDistances;

    BRepClass_FaceClassifier aFaceClassifier;
    double minDistance = 1e80;
    int k = -1;
    int k_min = -1;
    for(std::vector<TopoDS_Face>::const_iterator it = myFaces.cbegin(); it!=myFaces.cend(); it++)
    {
        const TopoDS_Face &curFace = *it;
        GeomAPI_ProjectPointOnSurf aProjector;
        BRepAdaptor_Surface adaptor(curFace,true);
        const GeomAdaptor_Surface &s_adaptor = adaptor.Surface();
        opencascade::handle<Geom_Surface> aGeomSurface =  s_adaptor.Surface();

        aProjector.Init(gp_Pnt(aPoint[0],aPoint[1],aPoint[2]),aGeomSurface);
        if(!aProjector.IsDone()) continue;

        const gp_Pnt &nst = aProjector.NearestPoint();
        double faceTolerance = 1.0*BRep_Tool::Tolerance(curFace);
        aFaceClassifier.Perform(curFace,nst,faceTolerance);
        TopAbs_State pointStatus = aFaceClassifier.State();
        if(pointStatus == TopAbs_IN)
        {
            k++;
            std::vector<double> aProjectedPoint{nst.X(),nst.Y(),nst.Z()};
            double curDistance = sqrt(pow(aPoint[0]-nst.X(),2)+pow(aPoint[1]-nst.Y(),2)+pow(aPoint[2]-nst.Z(),2));
            vecProjections.push_back(aProjectedPoint);
            vecDistances.push_back(curDistance);
            if(minDistance>=curDistance)
            {
                minDistance = curDistance;
                k_min=k;
            }
        }
    }
    if(k_min==-1)
    {
        cout<<"OCCFace::pointProjection()->____failure in projection: no valid projection found____"<<endl;
        return false;
    }
    if(minDistance >=5.0*myMaxDeflection) return false;
    for(int i=0; i<3; i++) theProjectedPoint[i] = vecProjections[k_min][i];
    return true;
}
