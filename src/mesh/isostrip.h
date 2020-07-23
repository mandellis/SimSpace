#ifndef ISOSTRIP_H
#define ISOSTRIP_H

struct isoStrip
{
    double xmin,xmax;

    isoStrip(double axmin = 0, double axmax = 100):xmin(axmin),xmax(axmax){;}
    isoStrip(const isoStrip &rhs) {xmin = rhs.xmin; xmax = rhs.xmax; }
    bool operator == (isoStrip &rhs) { if(xmin==rhs.xmin && xmax==rhs.xmax) return true; return false; }
    isoStrip operator = (const isoStrip &rhs) {xmin=rhs.xmin; xmax=rhs.xmax; return *this; }
    bool operator < (const isoStrip &rhs) const { if(center()<rhs.center()) return true; return false; }
    double center() const {return 0.5*(xmin+xmax);}
    double width() const {return xmax-xmin;}
    bool contains(double x) const { if (x>=xmin && x<=xmax) return true; return false; }
};

#endif // ISOSTRIP_H
