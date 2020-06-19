#ifndef CGAL_TOOLS_H
#define CGAL_TOOLS_H

#include <ng_meshvs_datasourceface.h>

#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

//! ------------------------------------------------------------------------
//! http://jamesgregson.blogspot.com/2012/05/example-code-for-building.html
//! A modifier creating a triangle with the incremental builder
//! ------------------------------------------------------------------------
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS>
{
public:

    std::vector<double> &coords;
    std::vector<int>    &tris;

    polyhedron_builder( std::vector<double> &_coords, std::vector<int> &_tris ) : coords(_coords), tris(_tris) {}
    void operator()( HDS& hds)
    {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( coords.size()/3, tris.size()/3 );

        // add the polyhedron vertices
        for( int i=0; i<(int)coords.size(); i+=3 )
        {
            B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
        }

        // add the polyhedron triangles
        for( int i=0; i<(int)tris.size(); i+=3 )
        {
            B.begin_facet();
            B.add_vertex_to_facet( tris[i+0] );
            B.add_vertex_to_facet( tris[i+1] );
            B.add_vertex_to_facet( tris[i+2] );
            B.end_facet();
        }

        // finish up the surface
        B.end_surface();
    }
};
#endif // CGAL_TOOLS_H
