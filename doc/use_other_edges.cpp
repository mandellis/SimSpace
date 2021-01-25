// ngscxx -c use_other_edges.cpp; ngsld use_other_edges.o -lnglib -lstl -lmesh

#include <stlgeom.hpp>
#include <meshing.hpp>

namespace netgen
{
  void STLMeshing (STLGeometry & geom,
                   Mesh & mesh);
}

using namespace netgen;

int main() {
  // read the edges in from somewhere... 
  Array<pair<Point<3>,Point<3>>> edges;
  edges.Append({{1.80522, 7.52768, 3.3329}, {5.80522, 7.52768, 3.3329}});
  ifstream ist("part1.stl");
  auto geo = STLGeometry::Load(ist);
  auto mesh = make_shared<Mesh>();
  MeshingParameters mp;
  cout << "STLMeshing done" << endl;
  mesh->SetGlobalH(mp.maxh);
  mesh->SetLocalH(geo->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
                  geo->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
                  0.3);

  for(int i=1; i<=geo->GetNE(); i++)
    // geo->GetTopEdge(i).SetStatus(ED_EXCLUDED);
    geo->GetTopEdge(i).SetStatus(ED_UNDEFINED);

  for(auto pair : edges)
    geo->GetTopEdge(geo->GetTopEdgeNum(geo->GetPointNum(pair.first),geo->GetPointNum(pair.second))).SetStatus(ED_CONFIRMED);
  cout << "Build edges" << endl;

  geo->BuildEdges();
  geo->edgesfound = 1;
  geo->MakeAtlas(*mesh);
  geo->CalcFaceNums();
  geo->AddFaceEdges();
  geo->LinkEdges();
  
  mesh->ClearFaceDescriptors();
  for (int i = 1; i <= geo->GetNOFaces(); i++)
    mesh->AddFaceDescriptor (FaceDescriptor (i, 1, 0, 0));
  mp.perfstepsstart = MESHCONST_MESHSURFACE;
  mp.perfstepsend = MESHCONST_OPTSURFACE;
  geo->GenerateMesh(mesh, mp);
  mesh->Save("a.vol.gz");
  return 0;
}
