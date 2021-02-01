#ifndef MESHUVPROJECTION_H
#define MESHUVPROJECTION_H

#include <ext/eigen/Eigen/Dense>
#include <TopoDS_Face.hxx>

class meshUVprojection
{
public:

    meshUVprojection();

    static bool projectUV(const TopoDS_Face &aFace, const Eigen::MatrixXd &V, Eigen::MatrixXd &VProj);
};

#endif // MESHUVPROJECTION_H
