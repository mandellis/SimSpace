//! ----------------
//! custom includes
//! ----------------
#include "particlesemitter.h"
#include "ng_meshvs_datasourceface.h"

//! ----
//! OCC
//! ----
#include <BRepTools.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_Surface.hxx>
#include <gp_Pnt.hxx>

//! ----
//! C++
//! ----
#include <random>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
particlesEmitter::particlesEmitter(const occHandle(Ng_MeshVS_DataSourceFace) &theFaceMesh,
                                   double theIntensity,
                                   double theInternalTime):
    myFaceMesh(theFaceMesh), myIntensity(theIntensity),myInternalTime(theInternalTime)
{
    //! --------------
    //! discrete face
    //! --------------
    myDiscreteFace = std::make_shared<meshFace>();
    myDiscreteFace->setGeometry(theFaceMesh);
}

//! ----------------------------
//! function: resetInternalTime
//! details:
//! ----------------------------
void particlesEmitter::resetInternalTime()
{
    myInternalTime = 0.0;
}

//! --------------------------
//! function: getInternalTime
//! details:
//! --------------------------
double particlesEmitter::getInternalTime() const
{
    return myInternalTime;
}

//! ----------------------
//! function: setLocation
//! details:
//! ----------------------
void particlesEmitter::setLocation(const occHandle(Ng_MeshVS_DataSourceFace) &theFaceMesh)
{
    myFaceMesh = theFaceMesh;

    //! --------------
    //! discrete face
    //! --------------
    myDiscreteFace = std::make_shared<meshFace>();
    myDiscreteFace->setGeometry(theFaceMesh);
}

//! -----------------------
//! function: setIntensity
//! details:
//! -----------------------
void particlesEmitter::setIntensity(double theIntensity)
{
    myIntensity = theIntensity;
}

//! -----------------------
//! function: getIntensity
//! details:
//! -----------------------
double particlesEmitter::getIntensity() const
{
    return myIntensity;
}

//! ---------------------------
//! function: getBirthPosition
//! details:
//! ---------------------------
bool particlesEmitter::sampleRandomBirthPositionOnSource(double &x, double &y, double &z) const
{
    double p[3];
    myDiscreteFace->randomPointOn(p);
    x = p[0]; y = p[1]; z = p[2];
    return true;
}

//! ----------------------------
//! function: generateParticles
//! details:
//! ----------------------------
void particlesEmitter::generateParticles(std::vector<particle> &newParticles, double time)
{
    double deltaT = time-myInternalTime;
    const double accuracy = deltaT/50.0;
    if(fabs(deltaT)>accuracy) return;     // no particle generation

    //! -------------------------------------------------
    //! time increase for the very first position update
    //! -------------------------------------------------
    double dt = deltaT/10.0;
    int NbParticles = int(myIntensity*deltaT);
    for(int i=0; i<NbParticles; i++)
    {
        double x,y,z;
        this->sampleRandomBirthPositionOnSource(x,y,z);

        //! ------------------------------
        //! sample velocity ... to do ...
        //! ------------------------------
        double vx = 1e5; double vy = 1e5; double vz = 1e5;

        x += vx*dt; y += vy*dt; z += vz*dt;
        particle aParticle(x,y,z,vx,vy,vz,9.11e-31,1.6e-19);
        newParticles.push_back(aParticle);
    }
    //! --------------------------
    //! set the new internal time
    //! --------------------------
    myInternalTime = time;
}
