// Copyright (c) 1997-2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14/Optimisation_basic/include/CGAL/Optimisation/Access_dimension_3.h $
// $Id: Access_dimension_3.h 0698f79 %aI Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_ACCESS_DIMENSION_3_H
#define CGAL_OPTIMISATION_ACCESS_DIMENSION_3_H

namespace CGAL {

// Class declaration
// =================
template < class R_ >
class Access_dimension_3;

// Class interface
// ===============
template < class R_ >
class Access_dimension_3 {
  public:
    // self
    typedef  R_                         R;
    typedef  Access_dimension_3<R>      Self;

    // types
    typedef  typename R::Point_3        Point;

    // unary function class types
    typedef  int                        result_type;
    typedef  Point                      argument_type;

    // creation
    Access_dimension_3( ) { }

    // operations
    int  operator() ( const Point& p) const { return p.dimension(); }
};

} //namespace CGAL

#endif // CGAL_OPTIMISATION_ACCESS_DIMENSION_3_H

// ===== EOF ==================================================================
