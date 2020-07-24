#ifndef ISOSTRIP_H
#define ISOSTRIP_H

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
    bool contains(double x) const { if (x>=vmin && x<=vmax) return true; return false; }
};

#endif // ISOSTRIP_H
