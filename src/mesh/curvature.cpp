//! custom includes
#include "curvature.h"

//! OCC
#include <BRepLProp_SLProps.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

#define SAFE_CURVATURE_VALUE 0.1

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
curvature::curvature(const TopoDS_Face &aFace, double *P):
    myFace(aFace)
{
    myPoint[0] = P[0];
    myPoint[1] = P[1];
    myPoint[2] = P[2];
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
curvature::curvature()
{
    myPoint[0] = myPoint[1] = myPoint[2] = 0.0;
    myFace = TopoDS_Face();
}

//! --------------------------------------------------------------
//! function: parameters
//! details:  found the (u,v) parameters of the point on the face
//! --------------------------------------------------------------
bool curvature::parameters(double *uv)
{
    double u,v;
    double tolerance = BRep_Tool::Tolerance(myFace);

    const occHandle(Geom_Surface) &aSurface = BRep_Tool::Surface(myFace);

    gp_Pnt P(myPoint[0], myPoint[1], myPoint[2]);
    bool isFound;
    isFound = GeomLib_Tool::Parameters(aSurface,P,tolerance,u,v);
    if(isFound)
    {
        uv[0] = u;
        uv[1] = v;
        return true;
    }
    else
    {
        //! --------------------------------------
        //! point not found withing the tolerance
        //! trying projection
        //! --------------------------------------
        GeomAPI_ProjectPointOnSurf projector(P,aSurface);
        if(projector.NbPoints()==0)
        {
            //! the projection has failed
            return false;
        }
        else
        {
            //! point projection success
            if(projector.LowerDistance()>tolerance)
            {
                //! if the point is too far from the face
                //! do not accept it
                return false;
            }
            else
            {
                projector.LowerDistanceParameters(u,v);
                uv[0] = u;
                uv[1] = v;
                return true;
            }
        }
    }
}

//! -----------------------
//! function: getCurvature
//! details:
//! -----------------------
bool curvature::getCurvature(double &curvature, typeOfCurvature type)
{
    double uv[2];
    bool isDone = this->parameters(uv);
    if(isDone)
    {
        int NbDerivatives = 1;
        BRepAdaptor_Surface aSurface(myFace,true);
        double tolerance = aSurface.Tolerance();
        double u = uv[0];
        double v = uv[1];
        BRepLProp_SLProps calc(aSurface,u,v,NbDerivatives,tolerance);

        if(calc.IsCurvatureDefined())
        {
            //! curvature defined
            switch(type)
            {
            case typeOfCurvature_gaussian: curvature = calc.GaussianCurvature(); break;
            case typeOfCurvature_min: curvature = calc.MinCurvature(); break;
            case typeOfCurvature_max: curvature = calc.MaxCurvature(); break;
            case typeOfCurvature_average: curvature = calc.MeanCurvature(); break;
            }
            return true;
        }
        else
        {
            //! curvature not defined
            return false;
        }
    }
    else
    {
        //! (u,v) parameters not found for some reason
        return false;
    }
}

//! --------------------------------------------------------------------------
//! function: getCurvature - static
//! details:  return false if:
//!           1) the (u,v) parameters of P on the face cannot be
//!  calculated,
//!              also using an orthogonal projection (the projection point P'
//!              does not exists, or it is too far from P)
//!           2) the curvature in P(P') is not defined
//!           The function calculates also the normal; if the normal cannot
//!           be calculated normal == NULL
//! --------------------------------------------------------------------------
bool curvature::getCurvature_static(const TopoDS_Face &aFace,
                                    double *aPoint,
                                    double &curvature,
                                    typeOfCurvature type,
                                    double *normal)
{
    if(!aFace.IsNull())
    {
        //! ----------------------------------------------------------
        //! calculate the (u,V) parameters of the point P on the face
        //! ----------------------------------------------------------
        double u,v;
        double tolerance = BRep_Tool::Tolerance(aFace);

        const occHandle(Geom_Surface) &aGeomSurface = BRep_Tool::Surface(aFace);

        gp_Pnt P(aPoint[0], aPoint[1], aPoint[2]);
        bool isFound;
        isFound = GeomLib_Tool::Parameters(aGeomSurface,P,tolerance,u,v);
        if(!isFound)
        {
            //! --------------------------------------
            //! point not found withing the tolerance
            //! trying projection
            //! --------------------------------------
            GeomAPI_ProjectPointOnSurf projector(P,aGeomSurface);
            if(projector.NbPoints()>0)
            {
                if(projector.LowerDistance()<tolerance)
                {
                    projector.LowerDistanceParameters(u,v);                    
                }
                else
                {
                    //! ----------------------------------
                    //! P' too far from P: P' is rejected
                    //! ----------------------------------
                    cout<<"____projection point too far from selected point: (u,v) coords not found: return safe curvature___"<<endl;
                    curvature = SAFE_CURVATURE_VALUE;
                    return false;
                }
            }
            else
            {
                //! ---------------------
                //! projection not found
                //! ---------------------
                cout<<"____projection point not found: (u,v) coords not found: returning safe curvature____"<<endl;
                curvature = SAFE_CURVATURE_VALUE;
                return false;
            }
        }

        //! ---------------------------------------
        //! (u,v) found: begin computing curvature
        //! ---------------------------------------
        int NbDerivatives = 2;
        BRepAdaptor_Surface aSurface(aFace,true);
        double tol = aSurface.Tolerance();

        BRepLProp_SLProps calc(aSurface,u,v,NbDerivatives,tol);

        if(calc.IsCurvatureDefined())
        {
            //! ----------------------------------
            //! the curvature in (u,v) is defined
            //! ----------------------------------
            switch(type)
            {
            case typeOfCurvature_gaussian: curvature = calc.GaussianCurvature(); break;
            case typeOfCurvature_max: curvature = calc.MaxCurvature(); break;
            case typeOfCurvature_min: curvature = calc.MinCurvature(); break;
            case typeOfCurvature_average: curvature = calc.MaxCurvature(); break;
            }
            //! --------------------------
            //! try to compute the normal
            //! --------------------------
            if(calc.IsNormalDefined())
            {
                normal[0] = calc.Normal().X();
                normal[1] = calc.Normal().Y();
                normal[2] = calc.Normal().Z();
            }
            else
            {
                normal = NULL;
            }
            return true;
        }
        else
        {
            //! --------------------------------------
            //! the curvature in (u,v) is not defined
            //! try to compute at least the normal
            //! --------------------------------------
            curvature = 0.0;
            if(calc.IsNormalDefined())
            {
                normal[0] = calc.Normal().X();
                normal[1] = calc.Normal().Y();
                normal[2] = calc.Normal().Z();
            }
            else
            {
                normal = NULL;
            }
            return false;
        }
    }
    return false; //! the face is null
}
