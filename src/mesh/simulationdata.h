#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H

#include <QMetaType>
#include <QMap>
#include <vector>

//! -----------------------------------------------------------------------
//! key:   global node ID (node ID of the node within the volume mesh
//!        for a 3D analysis, or within the surface mesh for a 2D analysis
//! value: can be scalar, or 2D/3D vectorial
//! -----------------------------------------------------------------------
typedef QMap<int,std::vector<double>> distribution;
Q_DECLARE_METATYPE(distribution)

//! --------------------------------------------------------------------
//! key:    time index - needs a linked list for linking the index to a
//!         (set,substep) pair or a physical time
//! value:  a distribution
//! --------------------------------------------------------------------
typedef QList<distribution> timeHistoryOfDistributions;
Q_DECLARE_METATYPE(timeHistoryOfDistributions)

#endif // SIMULATIONDATA_H
