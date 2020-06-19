#define TESTS_DIR "D:/Work/Qt/pro_25.1_OCC7.3.0_MSVC2015/tests"

//! custom includes
#include "testtools.h"

//! OCC
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <MeshVS_EntityType.hxx>

//! C++
#include <random>

testTools::testTools()
{
    ;
}


//! -----------------------------------------------------------------------
//! function: fakeScalarData
//! details:  generate fake scalar data withing the BB of the mesh
//! -----------------------------------------------------------------------
void testTools::fakeScalarData(const occHandle(MeshVS_DataSource) &aDS, const QString &fileName, int Npoints)
{
    QString path = QString(TESTS_DIR)+"\\"+fileName+"_fake.txt";
    ofstream file;
    file.open(qPrintable(path));

    TColStd_PackedMapOfInteger nodes = aDS->GetAllNodes();
    TColStd_MapIteratorOfPackedMapOfInteger nodeMap;
    double x_min = 1e99;
    double y_min = 1e99;
    double z_min = 1e99;
    double x_max = -1e99;
    double y_max = -1e99;
    double z_max = -1e99;
    for(nodeMap.Initialize(nodes);nodeMap.More();nodeMap.Next())
    {
        int nodeID = nodeMap.Key();
        MeshVS_EntityType type;
        int NbNodes;
        Standard_Real aCoordsBuf[3];
        TColStd_Array1OfReal Coords(*aCoordsBuf,1,3);
        aDS->GetGeom(nodeID,Standard_False,Coords,NbNodes,type);

        if (Coords(1)>=x_max) x_max = Coords(1);
        if (Coords(2)>=y_max) y_max = Coords(2);
        if (Coords(3)>=z_max) z_max = Coords(3);

        if (Coords(1)<=x_min) x_min = Coords(1);
        if (Coords(2)<=y_min) y_min = Coords(2);
        if (Coords(3)<=z_min) z_min = Coords(3);
    }

    double val;
    double Lx = x_max-x_min;
    cout<<"Lx = "<<Lx<<endl;

    for(int i=0; i<Npoints; i++)
    {
        double M = x_min; double N = x_max;
        double xr = M+rand()/(RAND_MAX/(N-M+1)+1);
        M = y_min; N = y_max;
        double yr = M+rand()/(RAND_MAX/(N-M+1)+1);
        M = z_min; N = z_max;
        double zr = M+rand()/(RAND_MAX/(N-M+1)+1);

        // x right
        if(xr>=x_min && xr<=x_min+Lx/2) val = 10;
        // x left
        if(xr>=x_min+Lx/2 && xr<=x_max) val = 100;

        file<<xr<<"\t"<<yr<<"\t"<<zr<<"\t"<<val<<endl;
    }
}


void testTools::fakeScalarData1(const occHandle(MeshVS_DataSource) &aDS, const QString &fileName)
{
    QString path = QString(TESTS_DIR)+"\\"+fileName+"_fake.txt";
    ofstream file;
    file.open(qPrintable(path));

    TColStd_PackedMapOfInteger nodes = aDS->GetAllNodes();
    TColStd_MapIteratorOfPackedMapOfInteger nodeMap;
    double x_min = 1e99;
    double y_min = 1e99;
    double z_min = 1e99;
    double x_max = -1e99;
    double y_max = -1e99;
    double z_max = -1e99;
    for(nodeMap.Initialize(nodes);nodeMap.More();nodeMap.Next())
    {
        int nodeID = nodeMap.Key();
        MeshVS_EntityType type;
        int NbNodes;
        Standard_Real aCoordsBuf[3];
        TColStd_Array1OfReal Coords(*aCoordsBuf,1,3);
        aDS->GetGeom(nodeID,Standard_False,Coords,NbNodes,type);

        if (Coords(1)>=x_max) x_max = Coords(1);
        if (Coords(2)>=y_max) y_max = Coords(2);
        if (Coords(3)>=z_max) z_max = Coords(3);

        if (Coords(1)<=x_min) x_min = Coords(1);
        if (Coords(2)<=y_min) y_min = Coords(2);
        if (Coords(3)<=z_min) z_min = Coords(3);
    }

    double val;
    double Lx = x_max-x_min;
    double Ly = y_max-y_min;
    double Lz = z_max-z_min;
    cout<<"Lx = "<<Lx<<" Ly= "<<Ly<<" Lz= "<<Lz<<endl;

    for(nodeMap.Initialize(nodes);nodeMap.More();nodeMap.Next())
    {
        int nodeID = nodeMap.Key();
        MeshVS_EntityType type;
        int NbNodes;
        Standard_Real aCoordsBuf[3];
        TColStd_Array1OfReal Coords(*aCoordsBuf,1,3);
        aDS->GetGeom(nodeID,Standard_False,Coords,NbNodes,type);

        // x left
        if(Coords(1)>=x_min && Coords(1)<=x_min+Lx/2)
        {
            // y left
            if(Coords(2)>=y_min && Coords(2)<=y_min+Ly/2)
            {
                //z left
                if(Coords(3)>=z_min && Coords(3)<=z_min+Lz/2)
                {
                    val=10;
                }
                // z right
                if(Coords(3)>=z_min+Lz/2 && Coords(3)<=z_max)
                {
                    val=20;
                }
            }
            // y right
            if(Coords(2)>=y_min+Ly/2 && Coords(2)<=y_max)
            {
                //z left
                if(Coords(3)>=z_min && Coords(3)<=z_min+Lz/2)
                {
                    val=30;
                }
                // z right
                if(Coords(3)>=z_min+Lz/2 && Coords(3)<=z_max)
                {
                    val=40;
                }
            }
        }

        // x right
        if(Coords(1)>=x_min+Lx/2 && Coords(1)<=x_max)
        {
            // y left
            if(Coords(2)>=y_min && Coords(2)<=y_min+Ly/2)
            {
                //z left
                if(Coords(3)>=z_min && Coords(3)<=z_min+Lz/2)
                {
                    val=50;
                }
                // z right
                if(Coords(3)>=z_min+Lz/2 && Coords(3)<=z_max)
                {
                    val=60;
                }
            }
            // y right
            if(Coords(2)>=y_min+Ly/2 && Coords(2)<=y_max)
            {
                //z left
                if(Coords(3)>=z_min && Coords(3)<=z_min+Lz/2)
                {
                    val=70;
                }
                // z right
                if(Coords(3)>=z_min+Lz/2 && Coords(3)<=z_max)
                {
                    val=80;
                }
            }
        }

        file<<Coords(1)<<"\t"<<Coords(2)<<"\t"<<Coords(3)<<"\t"<<val<<endl;
    }
}
