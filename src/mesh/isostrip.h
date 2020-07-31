#ifndef ISOSTRIP_H
#define ISOSTRIP_H

//! ------
//! point
//! ------
struct isoStripPoint
{
    double x,y,z;
    double val;
    isoStripPoint(double ax=0, double ay=0, double az=0, double aVal = 0) { x = ax; y = ay; z = az; val = aVal; }
    isoStripPoint(const isoStripPoint &aP) { x = aP.x; y = aP.y; z = aP.z; val = aP.val; }
    isoStripPoint operator =(const isoStripPoint &aP) { x = aP.x; y = aP.y; z = aP.z; val = aP.val; return *this;}
    bool operator == (const isoStripPoint &aP) const { if(x == aP.x && y == aP.y && z == aP.z) return true; return false; }
};

//! ---------
//! isostrip
//! ---------
struct isoStrip
{
    double vmin,vmax;

    isoStrip(double avmin = 0, double avmax = 100):vmin(avmin),vmax(avmax){;}
    isoStrip(const isoStrip &rhs) {vmin = rhs.vmin; vmax = rhs.vmax; }
    bool operator == (isoStrip &rhs) { if(vmin==rhs.vmin && vmax==rhs.vmax) return true; return false; }
    isoStrip operator = (const isoStrip &rhs) {vmin=rhs.vmin; vmax=rhs.vmax; return *this; }
    bool operator < (const isoStrip &rhs) const { if(center()<rhs.center()) return true; return false; }
    double center() const {return 0.5*(vmin+vmax);}
    double width() const {return vmax-vmin;}
    bool contains(double v) const { if (v>=vmin && v<=vmax) return true; return false; }
};

#endif // ISOSTRIP_H
