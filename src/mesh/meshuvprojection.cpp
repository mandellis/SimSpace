#include "meshuvprojection.h"
#include "occhandle.h"

#include <BRepAdaptor_Surface.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

//! ------------------------------------------------------------
//! function: projectUV
//! details:  project the vertices of a triangular mesh defined
//!           on the face aFace onto the reciprocal space (u,v)
//! ------------------------------------------------------------
bool projectUV(const TopoDS_Face &aFace, const Eigen::MatrixXd &V, Eigen::MatrixXd &VProj)
{
    if(aFace.IsNull()) return false;
    if(V.rows()==0) return false;

    VProj.resize(V.rows(),3);

    GeomAPI_ProjectPointOnSurf aProjector;
    BRepAdaptor_Surface adaptor(aFace,true);
    const GeomAdaptor_Surface &s_adaptor = adaptor.Surface();
    occHandle(Geom_Surface) aGeomSurface =  s_adaptor.Surface();
    for(int i=0; i<V.rows();i++)
    {
        gp_Pnt P(V(i,0),V(i,1),V(i,2));  // the point to project
        aProjector.Init(P, aGeomSurface);   // init the projection
        //! --------------------------------------------------------------------
        //! Returns the parameters (u,v) on the surface of the nearest computed
        //! orthogonal projection of the point. Exceptions StdFail_NotDone if
        //! projection fails (from OCC documentation)
        //! --------------------------------------------------------------------
        double u,v;
        aProjector.LowerDistanceParameters(u,v);
        VProj(i,0) = u;
        VProj(i,1) = v;
        VProj(i,2) = 0.0;
    }
    return true;
}
