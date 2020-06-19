//! custom includes
#include "contactparameters.h"
#include "geomtoolsclass.h"
#include "polygon.h"

//! OCC
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>

contactParameters::contactParameters()
{
    ;
}


//! -----------------------
//! function: trim_contact
//! details:
//! -----------------------


//! ----------------------------------------------------------------------------
//! function: calc_K
//! details:  provides a rough estimate of the contact stiffness (experimental)
//! ----------------------------------------------------------------------------
double contactParameters::calc_K(QList<occHandle(Ng_MeshVS_DataSourceFace)> master,QList<occHandle(Ng_MeshVS_DataSourceFace)> slave /*const opencascade::handle<Ng_MeshVS_DataSourceFace> &faceMesh*/)
{
    //! get elastic coefficients
    //! search for E and nu for each body material
    //! to do ... here static for TWO materials

    double E_1 = 200.0e5;
    double E_2 = 200.0e5;

    double nu_1 = 0.3;
    double nu_2 = 0.3;

    double E1_star = E_1/(1-pow(nu_1,2));
    double E2_star = E_2/(1-pow(nu_2,2));

    double E_star = sqrt(E1_star*E2_star);

    double A = 0;
    for(QList<occHandle(Ng_MeshVS_DataSourceFace)>::iterator it = master.begin(); it!= master.end(); ++it)
    {
        const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh = *it;
        A += calc_aveDiscretizedArea(faceMesh);
    }

    //! characteristic size
    double r = sqrt(A/3.141592654);

    //! contact stiffness
    double K = 2*E_star/r;

    return K;
}


//! ---------------------------------------------------------
//! function: calc_discretizedArea
//! details:  calculate the total area of a discretized face
//! ---------------------------------------------------------
double contactParameters::calc_discretizedArea(const opencascade::handle<Ng_MeshVS_DataSourceFace> &faceMesh)
{
    cout<<"contactParameters::calc_discretizedArea()->____function called____"<<endl;
    double A = 0.0;
    int NbNodes;
    double bufn[24];
    TColStd_Array1OfReal coords(*bufn,1,24);
    MeshVS_EntityType type;
    TColStd_PackedMapOfInteger eMap = faceMesh->GetAllElements();
    TColStd_MapIteratorOfPackedMapOfInteger eIter;

    for(eIter.Initialize(eMap);eIter.More();eIter.Next())
    {
        double elementArea;
        int globalElementID = eIter.Key();
        bool isOK = faceMesh->GetGeom(globalElementID,true,coords,NbNodes,type);
        if(isOK)
        {
            QList<QList<double>> points;
            for(int i=0; i<NbNodes; i++)
            {
                QList<double> P;
                int s = 3*i;
                P<<coords(s+1)<<coords(s+2)<<coords(s+3);
                points<<P;
            }
            elementArea = GeomToolsClass::polygon_area(points);
            A += elementArea;
            cout<<"contactParameters::calc_discretizedArea()->____area_=  "<<A<<endl;
        }
    }
    return A;
}

//! -------------------------------------------------------------------
//! function: calc_AverageDiscretizedArea
//! details:  calculate the average element area of a discretized face
//! -------------------------------------------------------------------
double contactParameters::calc_aveDiscretizedArea(const opencascade::handle<Ng_MeshVS_DataSourceFace> &faceMesh)
{
    double A = calc_discretizedArea(faceMesh);

    //! total number of elements
    int NE = faceMesh->GetAllElements().Extent();

    //! average area of an element
    double area = A/NE;
    return area;
}
