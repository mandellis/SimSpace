#ifndef PARTICLESEMITTER_H
#define PARTICLESEMITTER_H

//! ----------------
//! custom includes
//! ----------------
#include "occhandle.h"
#include "ng_meshvs_datasourceface.h"
#include "meshface.h"
#include "particle.h"
#include "poissonsolver.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>
#include <Geom_Surface.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <memory>
#include <iostream>
using namespace std;

class particlesEmitter
{
public:

    particlesEmitter(const occHandle(Ng_MeshVS_DataSourceFace) &theFaceMesh,
                     double theIntensity = 0.0,
                     double theInternalTime = 0.0);

    particlesEmitter(const particlesEmitter &rhs)
    {
        myFaceMesh = new Ng_MeshVS_DataSourceFace(rhs.myFaceMesh);
        myIntensity = rhs.myIntensity;
        myInternalTime = rhs.myInternalTime;
        myDiscreteFace = rhs.myDiscreteFace;
    }

    void operator = (const particlesEmitter &rhs)
    {
        myFaceMesh = new Ng_MeshVS_DataSourceFace(rhs.myFaceMesh);
        myIntensity = rhs.myIntensity;
        myInternalTime = rhs.myInternalTime;
        myDiscreteFace = rhs.myDiscreteFace;
    }

public:

    void setLocation(const occHandle(Ng_MeshVS_DataSourceFace) &theFaceMesh);
    void setIntensity(double theIntensity); // particles per unit time
    double getIntensity() const;
    double getInternalTime() const;
    void resetInternalTime();
    void generateParticles(std::vector<particle> &newParticles, double time);

private:

    occHandle(Ng_MeshVS_DataSourceFace) myFaceMesh;
    double myIntensity;
    double myInternalTime;
    std::shared_ptr<meshFace> myDiscreteFace;

private:

    bool sampleRandomBirthPositionOnSource(double &x, double &y, double &z) const;
    void setInternalTime(double time) { myInternalTime = time; }
};

#endif // PARTICLESEMITTER_H
