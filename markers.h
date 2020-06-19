#ifndef MARKERS_H
#define MARKERS_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! OCC
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <AIS_Shape.hxx>
#include <AIS_ColoredShape.hxx>

//! Qt
#include <QVector>

//! custom includes
#include "ais_spheremarker.h"
#include "ais_arrowmarker.h"
#include "ais_doublearrowmarker.h"
#include "ais_customtrihedron.h"

#include "handle_ais_doublearrowmarker_reg.h"
#include "handle_ais_trihedron_reg.h"
#include "handle_ais_customtrihedron_reg.h"
#include "handle_ais_arrowmarker_reg.h"
#include "handle_ais_curvedarrowmarker_reg.h"
#include "handle_ais_spheremarker_reg.h"

class markers
{
private:

public:

    markers(){;}

    //! -------------------------------------------------------
    //! P: base point of the arrow (centre of the base circle)
    //! V: direction, R: radius, H: height
    //! -------------------------------------------------------
    static AIS_ArrowMarker_handle_reg buildArrowMarker(const gp_Pnt &P, const gp_Dir &V, double D);

    //! ----------------------
    //! build two arrows =><=
    //! ----------------------
    static AIS_DoubleArrowMarker_handle_reg buildDoubleArrowMarker(const gp_Pnt &P, const gp_Dir &V, double D);

    //! ---------
    //! C center
    //! ---------
    static AIS_SphereMarker_handle_reg buildSphereMarker(gp_Pnt C, double radius);

    //! -------------------------
    //! build a custom trihedron
    //! -------------------------
    static AIS_CustomTrihedron_handle_reg buildCustomTrihedron(QVector<double> Origin,
                                                               QVector<QVector<double>> directionalData,
                                                               double axisLength);

    //! ----------------------------------
    //! build a trihedron using OCC stuff
    //! ----------------------------------
    static AIS_Trihedron_handle_reg buildTrihedron(QVector<double> Origin, QVector<QVector<double>> directionalData, int axisLength=20);

    //! ---------------------
    //! build a curved arrow
    //! ---------------------
    static AIS_CurvedArrowMarker_handle_reg buildCurvedArrow(const gp_Ax2 &axes = gp_Ax2(), double Rin = 20, bool invert=false);
};

#endif // MARKERS_H
