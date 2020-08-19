#ifndef ISOSTRIP_H
#define ISOSTRIP_H

//! ------
//! point
//! ------
struct isoStripPoint
{
    double x,y,z;
    double val;
    int qualify;
    isoStripPoint(double ax=0, double ay=0, double az=0, double aVal = 0, int aQualify = 0)
    {
        x = ax; y = ay; z = az; val = aVal;
        qualify = aQualify;
    }
    isoStripPoint(const isoStripPoint &aP)
    {
        x = aP.x; y = aP.y; z = aP.z;
        val = aP.val;
        qualify = aP.qualify;
    }
    isoStripPoint operator = (const isoStripPoint &aP)
    {
        x = aP.x; y = aP.y; z = aP.z;
        val = aP.val;
        qualify = aP.qualify;
        return *this;
    }
    bool operator == (const isoStripPoint &aP) const
    {
        if(x == aP.x && y == aP.y && z == aP.z && qualify == aP.qualify) return true;
        return false;
    }

    void setTop() { qualify = 1; }
    void setBottom() { qualify = -1; }
    bool isBottom() const { return (qualify == -1? true: false); }
    bool isTop() const { return (qualify == 1? true: false); }
    int type() const { return qualify; }
};

/*
struct isoStripPoint
{
    double x,y,z;
    double val;
    int qualify;
    isoStripPoint(double ax=0, double ay=0, double az=0, double aVal = 0, int aQualify = 0) { x = ax; y = ay; z = az; val = aVal; qualify = aQualify; }
    isoStripPoint(const isoStripPoint &aP)
    {
        x = aP.x; y = aP.y; z = aP.z;
        val = aP.val;
        qualify = aP.qualify;
    }
    isoStripPoint operator =(const isoStripPoint &aP)
    {
        x = aP.x; y = aP.y; z = aP.z;
        val = aP.val;
        qualify=aP.qualify;
        return *this;
    }
    bool operator == (const isoStripPoint &aP) const
    {
        if(x == aP.x && y == aP.y && z == aP.z && val == aP.val && qualify == aP.qualify) return true;
        //if(x == aP.x && y == aP.y && z == aP.z) return true;
        return false;
    }
    inline void setTop() { qualify = 1; }
    inline void setBottom() { qualify = -1; }
    inline bool isBottom() const { return (qualify == -1? true: false); }
    inline bool isTop() const { return (qualify==1? true: false); }
    inline int type() const { return qualify; }
};
*/

//! ---------
//! isostrip
//! ---------
struct isoStrip
{
    double vmin,vmax;

    isoStrip(double avmin = -100, double avmax = 100):vmin(avmin),vmax(avmax){;}
    isoStrip(const isoStrip &rhs) {vmin = rhs.vmin; vmax = rhs.vmax; }
    bool operator == (isoStrip &rhs) { if(vmin==rhs.vmin && vmax==rhs.vmax) return true; return false; }
    isoStrip operator = (const isoStrip &rhs) {vmin=rhs.vmin; vmax=rhs.vmax; return *this; }
    bool operator < (const isoStrip &rhs) const { if(center()<rhs.center()) return true; return false; }
    double center() const {return 0.5*(vmin+vmax);}
    double width() const {return vmax-vmin;}
    bool contains(double v) const { if (v>=vmin && v<=vmax) return true; return false; }

    bool isPointAtBottom(const isoStripPoint &aPoint) { return (aPoint.val==vmin? true: false); }
    bool isPointAtTop(const isoStripPoint &aPoint) { return (aPoint.val==vmax? true: false); }
};

#endif // ISOSTRIP_H
