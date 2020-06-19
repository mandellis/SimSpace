#ifndef MESHPOINT_H
#define MESHPOINT_H

//! -------------
//! a mesh point
//! -------------
struct meshPoint
{
    int ID;
    double x,y,z;

    meshPoint(double xP, double, yP, double zP, int pointID = -1)
    {
        ID = pointID;
        x = xP; y = yP; z = zP;
    }

    meshPoint()
    {
        ID = -1;
        x = y = z = 0.0;
    }

    inline bool operator == (const meshPoint &rhs)
    {
        if(x==rhs.x && y==rhs.y && z == rhs.z) return true;
        return false;
    }

    /*
    inline bool operator < (const meshPoint &p0, const meshPoint &p1) const
    {
        if(pow(p0.x,2)+pow(p0.y,2)+pow(p0.z,2)+p0.x+2*p0.y+3*p0.z<
                pow(p1.x,2)+pow(p1.y,2)+pow(p1.z,2)+p1.x+2*p1.y+3*p1.z) return true;
        return false;
    }
    */

    inline bool operator < (const meshPoint &p1) const
    {
        if(pow(x,2)+pow(y,2)+pow(z,2)+x+2*y+3*z<
                pow(p1.x,2)+pow(p1.y,2)+pow(p1.z,2)+p1.x+2*p1.y+3*p1.z) return true;
        return false;
    }
};

#endif // MESHPOINT_H
