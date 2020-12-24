//#define WRITE_DIAG_FILE

//! ----------------
//! custom includes
//! ----------------
#include "mshConvert.h"

//! -------
//! libigl
//! -------
#include <libigl/include/igl/readMSH.h>

//! -------
//! pymesh
//! -------
#include <MshLoader.h>

//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>
using namespace std;

//! ---------------------------------------------------------------------------------
//! function: constructor
//! details:  simply a container for libraries able to read .msh files (text or bin)
//!           here libigl is used
//! ---------------------------------------------------------------------------------
mshConvert::mshConvert(const QString &absoluteFilePath):myMshInputFile(absoluteFilePath)
{
    ;
}

//! ------------------------------------
//! function: toEigen
//! details:  convert to Eigen matrices
//! ------------------------------------
bool mshConvert::toEigen(Eigen::MatrixXd &V, Eigen::MatrixXi &T)
{
    try
    {
        //! -------------------
        //! nodes and elements
        //! -------------------
        igl::readMSH(myMshInputFile.toStdString(),V,T);
    }
    catch(...)
    {
        cerr<<"____error in converting the .msh file____"<<endl;
        return false;
    }

#ifdef WRITE_DIAG_FILE
    //! ----------------
    //! diagnostic file
    //! ----------------
    FILE *f = fopen("D:/eigenMesh.txt","w");
    fprintf(f,"**** nodes ****\n");

    int NbNodes = V.rows();
    for(int i=0;i<NbNodes;i++)
    {
        double x = V(i,0);
        double y = V(i,1);
        double z = V(i,2);

        cout<<"____"<<i+1<<" (x,y,z) ="<<"("<<x<<","<<y<<","<<z<<")____"<<endl;
        fprintf(f,"%d\t%lf\t%lf\t%lf\n",i+1,x,y,z);
    }

    fprintf(f,"**** elements ****\n");

    int NbElements = T.rows();
    for(int i=0;i<NbElements;i++)
    {
        int i1 = T(i,0)+1;
        int i2 = T(i,1)+1;
        int i3 = T(i,2)+1;
        int i4 = T(i,3)+1;

        cout<<"____"<<i+1<<" (i1,i2,i3,i4)= "<<"("<<i1<<","<<i2<<","<<i3<<","<<i4<<")____"<<endl;

        fprintf(f,"%d\t%d\t%d\t%d\t%d\n",i+1,i1,i2,i3,i4);
    }
    fclose(f);
    cout<<"____diagnostic file written____"<<endl;
#endif

    return true;
}
