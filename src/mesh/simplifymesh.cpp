//! ----------------
//! custom includes
//! ----------------
#include "simplifymesh.h"

//! ----
//! OCC
//! ----
#include <OSD_Path.hxx>
#include <RWStl.hxx>
#include "extendedrwstl.h"
#include <StlMesh_Mesh.hxx>

#define NUMBER_OF_ITERATIONS 5000
#define MESH_UPDATE_RATE 5

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
SimplifyMesh::SimplifyMesh():fixedPoints(std::vector<mesh::tolerantPoint>())
{
    ;
}

//! -------------------------
//! function: setFixedPoints
//! details:
//! -------------------------
void SimplifyMesh::setFixedPoints(const std::vector<mesh::tolerantPoint> &points)
{
    myTolerance = 0;
    for(std::vector<mesh::tolerantPoint>::const_iterator it = points.cbegin(); it!=points.cend(); it++)
    {
        const mesh::tolerantPoint &aPoint = *it;
        fixedPoints.push_back(aPoint);
        myTolerance +=aPoint.eps;
    }
    myTolerance /= points.size();
}

//! ------------------------------------------------------
//! Main simplification function
//!
//! target_count  : target nr. of triangles
//! agressiveness : sharpness to increase the threshold.
//!                 5..8 are good numbers
//!                 more iterations yield higher quality
//! -----------------------------------------------------
double calc_threshold(int iteration)
{
    double val = 1e-10*pow(double(iteration+1),3);
    return val;
}

void SimplifyMesh::simplify_mesh(int target_count, double agressiveness, bool verbose)
{
    // init
    loopi(0,triangles.size())
    {
        triangles[i].deleted=0;
    }

    //! --------------------
    //! main iteration loop
    //! --------------------
    int deleted_triangles=0;
    std::vector<int> deleted0,deleted1;
    int triangle_count=int(triangles.size());
    for (int iteration = 0; iteration < NUMBER_OF_ITERATIONS; iteration++)
    {
        if(triangle_count-deleted_triangles<=target_count)break;

        // update mesh once in a while
        if(iteration % MESH_UPDATE_RATE==0)
        {
            update_mesh(iteration);
        }

        // clear dirty flag
        loopi(0,triangles.size()) triangles[i].dirty=0;

        //! --------------------------------------------------------------
        //! All triangles with edges below the threshold will be removed
        //! The following numbers works well for most models.
        //! If it does not, try to adjust the 3 parameters
        //! --------------------------------------------------------------
        //double threshold = 0.000000001*pow(double(iteration+3),agressiveness);
        double threshold = calc_threshold(iteration);

        // target number of triangles reached? Then break
        if ((verbose) && (iteration% MESH_UPDATE_RATE ==0))
        {
            printf("iteration %d - triangles %d threshold %g\n",iteration,triangle_count-deleted_triangles, threshold);
        }

        // remove vertices & mark deleted triangles
        loopi(0,triangles.size())
        {
            Triangle &t=triangles[i];
            if(t.err[3]>threshold) continue;
            if(t.deleted) continue;
            if(t.dirty) continue;

            loopj(0,3)if(t.err[j]<threshold)
            {

                int i0=t.v[ j     ]; Vertex &v0 = vertices[i0];
                int i1=t.v[(j+1)%3]; Vertex &v1 = vertices[i1];

                //! -------------
                //! Border check
                //! -------------
                //! https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification/issues/14
                //if(v0.border != v1.border)  continue;
                if(v0.border || v1.border) continue;

                // Compute vertex to collapse to
                vec3f p;
                calculate_error(i0,i1,p);
                deleted0.resize(v0.tcount); // normals temporarily
                deleted1.resize(v1.tcount); // normals temporarily
                // don't remove if flipped
                if(flipped(p,i0,i1,v0,v1,deleted0)) continue;

                if(flipped(p,i1,i0,v1,v0,deleted1)) continue;

                // not flipped, so remove edge
                v0.p=p;
                v0.q=v1.q+v0.q;
                int tstart=refs.size();

                update_triangles(i0,v0,deleted0,deleted_triangles);
                update_triangles(i0,v1,deleted1,deleted_triangles);

                int tcount=refs.size()-tstart;

                if(tcount<=v0.tcount)
                {
                    // save ram
                    if(tcount)memcpy(&refs[v0.tstart],&refs[tstart],tcount*sizeof(Ref));
                }
                else
                    // append
                    v0.tstart=tstart;

                v0.tcount=tcount;
                break;
            }
            // done?
            if(triangle_count-deleted_triangles<=target_count)break;
        }
    }
    // clean up mesh
    compact_mesh();
}

//! ---------------------------------
//! function: simplify_mesh_lossless
//! details:
//! ---------------------------------
void SimplifyMesh::simplify_mesh_lossless(double threshold, bool verbose)
{
    // init
    loopi(0,triangles.size()) triangles[i].deleted=0;

    // main iteration loop
    int deleted_triangles=0;
    std::vector<int> deleted0,deleted1;

    //!int triangle_count=triangles.size();

    for (int iteration = 0; iteration < NUMBER_OF_ITERATIONS; iteration ++)
    {
        // update mesh constantly
        update_mesh(iteration);

        // clear dirty flag
        loopi(0,triangles.size()) triangles[i].dirty=0;

        //! -------------------------------------------------------------
        //! All triangles with edges below the threshold will be removed
        //! The following numbers works well for most models.
        //! If it does not, try to adjust the 3 parameters
        //! -------------------------------------------------------------
        //double threshold = DBL_EPSILON; //1.0E-3 EPS;
        if (verbose) {
            //printf("lossless iteration %d\n", iteration);
            printf("lossless iteration %d - triangles %d threshold %g\n",iteration,triangles.size(), threshold);
        }

        // remove vertices & mark deleted triangles
        loopi(0,triangles.size())
        {
            Triangle &t=triangles[i];
            if(t.err[3]>threshold) continue;
            if(t.deleted) continue;
            if(t.dirty) continue;

            loopj(0,3)if(t.err[j]<threshold)
            {
                int i0=t.v[ j     ]; Vertex &v0 = vertices[i0];
                int i1=t.v[(j+1)%3]; Vertex &v1 = vertices[i1];

                //! -------------
                //! Border check
                //! -------------
                //! https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification/issues/14
                //if(v0.border != v1.border)  continue;
                if(v0.border || v1.border) continue;

                // Compute vertex to collapse to
                vec3f p;
                calculate_error(i0,i1,p);

                deleted0.resize(v0.tcount); // normals temporarily
                deleted1.resize(v1.tcount); // normals temporarily

                // don't remove if flipped
                if( flipped(p,i0,i1,v0,v1,deleted0) ) continue;
                if( flipped(p,i1,i0,v1,v0,deleted1) ) continue;

                // not flipped, so remove edge
                v0.p=p;
                v0.q=v1.q+v0.q;
                int tstart=refs.size();

                update_triangles(i0,v0,deleted0,deleted_triangles);
                update_triangles(i0,v1,deleted1,deleted_triangles);

                int tcount=refs.size()-tstart;

                if(tcount<=v0.tcount)
                {
                    // save ram
                    if(tcount)memcpy(&refs[v0.tstart],&refs[tstart],tcount*sizeof(Ref));
                }
                else
                    // append
                    v0.tstart=tstart;

                v0.tcount=tcount;
                break;
            }
        }
        if(deleted_triangles<=0)break;
        deleted_triangles=0;
    }
    // for each iteration
    // clean up mesh
    compact_mesh();
}

//! --------------------------------------------------------------
//! function: flipped
//! details:  Check if a triangle flips when this edge is removed
//! --------------------------------------------------------------
bool SimplifyMesh::flipped(vec3f p,int i0,int i1,Vertex &v0,Vertex &v1,std::vector<int> &deleted)
{
    loopk(0,v0.tcount)
    {
        Triangle &t=triangles[refs[v0.tstart+k].tid];
        if(t.deleted)continue;

        int s=refs[v0.tstart+k].tvertex;
        int id1=t.v[(s+1)%3];
        int id2=t.v[(s+2)%3];

        if(id1==i1 || id2==i1) // delete ?
        {
            deleted[k]=1;
            continue;
        }
        vec3f d1 = vertices[id1].p-p; d1.normalize();
        vec3f d2 = vertices[id2].p-p; d2.normalize();
        if(fabs(d1.dot(d2))>0.99999) return true;
        vec3f n;
        n.cross(d1,d2);
        n.normalize();
        deleted[k]=0;
        if(n.dot(t.n)<0.2) return true;
    }
    return false;
}

//! ---------------------------------------------------------------------
//! Update triangle connections and edge error after a edge is collapsed
//! ---------------------------------------------------------------------
void SimplifyMesh::update_triangles(int i0,Vertex &v,std::vector<int> &deleted,int &deleted_triangles)
{
    vec3f p;
    loopk(0,v.tcount)
    {
        Ref &r=refs[v.tstart+k];
        Triangle &t=triangles[r.tid];
        if(t.deleted)continue;
        if(deleted[k])
        {
            t.deleted=1;
            deleted_triangles++;
            continue;
        }
        t.v[r.tvertex]=i0;
        t.dirty=1;
        t.err[0]=calculate_error(t.v[0],t.v[1],p);
        t.err[1]=calculate_error(t.v[1],t.v[2],p);
        t.err[2]=calculate_error(t.v[2],t.v[0],p);
        t.err[3]=min(t.err[0],min(t.err[1],t.err[2]));
        refs.push_back(r);
    }
}

//! ---------------------------------------------------------------
//! compact triangles, compute edge error and build reference list
//! ---------------------------------------------------------------
void SimplifyMesh::update_mesh(int iteration)
{
    if(iteration>0) // compact triangles
    {
        int dst=0;
        loopi(0,triangles.size())
                if(!triangles[i].deleted)
        {
            triangles[dst++]=triangles[i];
        }
        triangles.resize(dst);
    }

    //! --------------------------------------------------------
    //! Init Quadrics by Plane & Edge Errors
    //!
    //! required at the beginning ( iteration == 0 )
    //! recomputing during the simplification is not required,
    //! but mostly improves the result for closed meshes
    //! -------------------------------------------------------
    if( iteration == 0 )
    {
        loopi(0,vertices.size())
                vertices[i].q=SymmetricMatrix(0.0);

        loopi(0,triangles.size())
        {
            Triangle &t=triangles[i];
            vec3f n,p[3];
            loopj(0,3) p[j]=vertices[t.v[j]].p;
            n.cross(p[1]-p[0],p[2]-p[0]);
            n.normalize();
            t.n=n;
            loopj(0,3) vertices[t.v[j]].q =
                    vertices[t.v[j]].q+SymmetricMatrix(n.x,n.y,n.z,-n.dot(p[0]));
        }
        loopi(0,triangles.size())
        {
            // Calc Edge Error
            Triangle &t=triangles[i];vec3f p;
            loopj(0,3) t.err[j]=calculate_error(t.v[j],t.v[(j+1)%3],p);
            t.err[3]=min(t.err[0],min(t.err[1],t.err[2]));
        }
    }

    //! -----------------------
    //! Init Reference ID list
    //! -----------------------
    loopi(0,vertices.size())
    {
        vertices[i].tstart=0;
        vertices[i].tcount=0;
    }
    loopi(0,triangles.size())
    {
        Triangle &t=triangles[i];
        loopj(0,3) vertices[t.v[j]].tcount++;
    }
    int tstart=0;
    loopi(0,vertices.size())
    {
        Vertex &v=vertices[i];
        v.tstart=tstart;
        tstart+=v.tcount;
        v.tcount=0;
    }

    //! -----------------
    //! Write References
    //! -----------------
    refs.resize(triangles.size()*3);
    loopi(0,triangles.size())
    {
        Triangle &t=triangles[i];
        loopj(0,3)
        {
            Vertex &v=vertices[t.v[j]];
            refs[v.tstart+v.tcount].tid=i;
            refs[v.tstart+v.tcount].tvertex=j;
            v.tcount++;
        }
    }

    //! ------------------------------------------
    //! Identify boundary : vertices[].border=0,1
    //! ------------------------------------------
    if( iteration == 0 )
    {
        std::vector<int> vcount,vids;

        loopi(0,vertices.size())
        {
            if(vertices[i].border==1) continue; //! GB
            vertices[i].border=0;
        }

        loopi(0,vertices.size())
        {
            Vertex &v=vertices[i];
            vcount.clear();
            vids.clear();
            loopj(0,v.tcount)
            {
                int k=refs[v.tstart+j].tid;
                Triangle &t=triangles[k];
                loopk(0,3)
                {
                    int ofs=0,id=t.v[k];
                    while(ofs<vcount.size())
                    {
                        if(vids[ofs]==id)break;
                        ofs++;
                    }
                    if(ofs==vcount.size())
                    {
                        vcount.push_back(1);
                        vids.push_back(id);
                    }
                    else
                        vcount[ofs]++;
                }
            }
            loopj(0,vcount.size()) if(vcount[j]==1)
                    vertices[vids[j]].border=1;
        }
    }
}

//! ------------------------------------
//! Finally compact mesh before exiting
//! ------------------------------------
void SimplifyMesh::compact_mesh()
{
    int dst=0;
    loopi(0,vertices.size())
    {
        vertices[i].tcount=0;
    }
    loopi(0,triangles.size())
            if(!triangles[i].deleted)
    {
        Triangle &t=triangles[i];
        triangles[dst++]=t;
        loopj(0,3)vertices[t.v[j]].tcount=1;
    }
    triangles.resize(dst);
    dst=0;
    loopi(0,vertices.size())
            if(vertices[i].tcount)
    {
        vertices[i].tstart=dst;
        vertices[dst].p=vertices[i].p;
        dst++;
    }
    loopi(0,triangles.size())
    {
        Triangle &t=triangles[i];
        loopj(0,3)t.v[j]=vertices[t.v[j]].tstart;
    }
    vertices.resize(dst);
}

//! -------------------------------------------
//! function: vertex_error()
//! details:  error between vertex and Quadric
//! --------------------------------------------
double SimplifyMesh::vertex_error(SymmetricMatrix q, double x, double y, double z)
{
    return   q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y
            + 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
}

//! -----------------------------
//! function: calculate_error
//! details:  Error for one edge
//! -----------------------------
double SimplifyMesh::calculate_error(int id_v1, int id_v2, vec3f &p_result)
{
    // compute interpolated vertex

    SymmetricMatrix q = vertices[id_v1].q + vertices[id_v2].q;
    bool border = vertices[id_v1].border & vertices[id_v2].border;
    double error=0;
    double det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);
    if ( det != 0 && !border )
    {

        // q_delta is invertible
        p_result.x = -1/det*(q.det(1, 2, 3, 4, 5, 6, 5, 7 , 8));	// vx = A41/det(q_delta)
        p_result.y =  1/det*(q.det(0, 2, 3, 1, 5, 6, 2, 7 , 8));	// vy = A42/det(q_delta)
        p_result.z = -1/det*(q.det(0, 1, 3, 1, 4, 6, 2, 5,  8));	// vz = A43/det(q_delta)

        error = vertex_error(q, p_result.x, p_result.y, p_result.z);
    }
    else
    {
        // det = 0 -> try to find best result
        vec3f p1=vertices[id_v1].p;
        vec3f p2=vertices[id_v2].p;
        vec3f p3=(p1+p2)/2;
        double error1 = vertex_error(q, p1.x,p1.y,p1.z);
        double error2 = vertex_error(q, p2.x,p2.y,p2.z);
        double error3 = vertex_error(q, p3.x,p3.y,p3.z);
        error = min(error1, min(error2, error3));
        if (error1 == error) p_result=p1;
        if (error2 == error) p_result=p2;
        if (error3 == error) p_result=p3;
    }
    return error;
}

//! -------------------
//! function: load_Stl
//! details:
//! -------------------
void SimplifyMesh::load_Stl(const char* fileName, int &NbNodes, int &NbElements)
{
    cout<<"SimplifyMesh::load_Stl()->____function called. File: "<<fileName<<"____"<<endl;
    vertices.clear();
    triangles.clear();

    Triangle t;
    Vertex v;

    OSD_Path filePath;
    TCollection_AsciiString aName(fileName);
    filePath.SetName(aName);

    const occHandle(StlMesh_Mesh) &StlMesh = ExtendedRWStl::ReadAscii(filePath,triangleTagMap);

    const TColgp_SequenceOfXYZ &seqPoints = StlMesh->Vertices();
    const StlMesh_SequenceOfMeshTriangle &seqTriangles = StlMesh->Triangles();

    //! ----------------------
    //! fill the vertex array
    //! ----------------------
    //FILE *f = fopen("D:/testboundary.txt","w");
    for(TColgp_SequenceOfXYZ::const_iterator it = seqPoints.cbegin(); it != seqPoints.cend(); ++it)
    {
        gp_XYZ PXYZ = *it;
        v.border = 0;           //! the point can be modified
        v.p.x=PXYZ.X();
        v.p.y=PXYZ.Y();
        v.p.z=PXYZ.Z();
        mesh::tolerantPoint curPoint(PXYZ.X(),PXYZ.Y(),PXYZ.Z(),myTolerance);
        if(std::find(fixedPoints.begin(), fixedPoints.end(), curPoint)!= fixedPoints.end())
        {
            cout<<"SimplifyMesh::load_Stl()->____fixing point: ("<<PXYZ.X()<<", "<<PXYZ.Y()<<", "<<PXYZ.Z()<<")____"<<endl;
            v.border = 1;  //! the point cannot be modified
            //fprintf(f,"%lf\t%lf\t%lf\n",PXYZ.X(),PXYZ.Y(),PXYZ.Z());
        }

        //! ------------
        //! for testing
        //! ------------
        //v.border = 0;

        vertices.push_back(v);
    }
    //fclose(f);

    //! ------------------------
    //! fill the triangle array
    //! ------------------------
    int c = 1;
    for(StlMesh_SequenceOfMeshTriangle::const_iterator it = seqTriangles.cbegin(); it != seqTriangles.cend(); ++it)
    {
        const occHandle(StlMesh_MeshTriangle) &aTriangle = *it;
        int V1,V2,V3;
        double n1,n2,n3;
        aTriangle->GetVertexAndOrientation(V1,V2,V3,n1,n2,n3);

        t.v[0] = V1-1;
        t.v[1] = V2-1;
        t.v[2] = V3-1;
        t.n.x=n1;
        t.n.y=n2;
        t.n.z=n3;
        t.tag = triangleTagMap.value(c);

        triangles.push_back(t);
        //cout<<"____loading triangle ("<<V1<<", "<<V2<<", "<<V3<<") surface tag: "<<t.tag<<"____"<<endl;
        c++;
    }
    NbNodes = int(vertices.size());
    NbElements = int(triangles.size());
    cout<<"SimplifyMesh::load_Stl()->____stl file loaded. Points "<<vertices.size()<<" triangles: "<<triangles.size()<<"____"<<endl<<endl;
}

//! --------------------
//! function: write_Stl
//! details:
//! --------------------
void SimplifyMesh::write_Stl(const char *fileName)
{
    cout<<"SimplifyMesh::write_Stl()->____function called. File name: \""<<fileName<<"____"<<endl;

    occHandle(StlMesh_Mesh) StlMesh = new StlMesh_Mesh();
    StlMesh->AddDomain();
    for(std::vector<Vertex>::iterator it = vertices.begin(); it!=vertices.end(); ++it)
    {
        const Vertex &aVertex = *it;
        StlMesh->AddOnlyNewVertex(aVertex.p.x,aVertex.p.y,aVertex.p.z);
    }

    //! -------------------------------------------------------------------
    //! some triangles have been eliminated, so the (triangleRang,tag) map
    //! must be recomputed
    //! -------------------------------------------------------------------
    QMap<int,int> triangleTagMapPurged;
    for(std::vector<Triangle>::iterator it = triangles.begin(); it!=triangles.end(); ++it)
    {
        const Triangle &aTriangle = *it;
        if(!aTriangle.deleted)
        {
            int rank = StlMesh->AddTriangle(aTriangle.v[0]+1,aTriangle.v[1]+1,aTriangle.v[2]+1,aTriangle.n.x,aTriangle.n.y,aTriangle.n.z);
            //! recompute the map (triangle nb, triangle tag);
            //cout<<"____inserting in the new map: triangle nb: "<<rank<<" tag: "<<aTriangle.tag<<"____"<<endl;
            triangleTagMapPurged.insert(rank,aTriangle.tag);
        }
    }
    //cout<<"____"<<StlMesh->Triangles().Length()<<" triangles "<<StlMesh->Vertices().Length()<<" vertices added to StlMesh_Mesh____"<<endl;

    //! -----------------------------
    //! write the extended .stl file
    //! -----------------------------
    OSD_Path filePath;
    TCollection_AsciiString aName(fileName);
    filePath.SetName(aName);

    ExtendedRWStl::WriteAscii(StlMesh,triangleTagMapPurged,filePath);
}
