#include <array>
#include <cmath>
#include <algorithm>

#include <GEL/CGLA/CGLA.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/graph_io.h>
#include <GEL/Geometry/graph_util.h>
#include <GEL/Geometry/SphereDelaunay.h>

using namespace Geometry;
using namespace CGLA;
using namespace HMesh;
using namespace std;
using NodeID = AMGraph::NodeID;
using NodeSet = AMGraph::NodeSet;

// Initialize global arrays

map<NodeID, int> val2deg;
map<pair<NodeID,NodeID>, int> branchdeg;
map<pair<NodeID,NodeID>, HMesh::FaceID> branchface;
map<pair<NodeID,NodeID>, HMesh::FaceID> branch_best_face;
map<pair<NodeID,NodeID>, HMesh::VertexID> branch_best_vertex;
map<pair<NodeID,NodeID>, HMesh::VertexID> one_ring_vertex;
map<pair<NodeID,NodeID>, CGLA::Vec3d> branch2vert;
map<FaceID, VertexID> face_vertex;
map<FaceID, VertexID> one_ring_face_vertex;
map<FaceID, int> val2_faces;

void clear_global_arrays() {
  val2deg.clear();
  branchdeg.clear();
  branchface.clear();
  branch_best_face.clear();
  branch_best_vertex.clear();
  branch2vert.clear();
  face_vertex.clear();
  one_ring_face_vertex.clear();
  return;
}

//Graph util functions

vector<NodeID> next_neighbours(const Geometry::AMGraph3D& g, NodeID prev, NodeID curr) {

    vector<NodeID> neighbour_list;
    auto N = g.neighbors(curr);
    for (auto next: N) {
        if(next != prev)
            neighbour_list.push_back(next);
    }
    return neighbour_list;

}

NodeID next_jn(const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {
  auto N = next_neighbours(g, n, nn);
  NodeID curr_node = nn;
  NodeID prev_node = n;

  while(true) {
      auto curr_nbs = next_neighbours(g, prev_node, curr_node);
      if(curr_nbs.size() != 1) {
          return curr_node;
      }
      prev_node = curr_node;
      curr_node = curr_nbs[0];
  }

}

double graph_length(const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {
  vector<NodeID> branch_path;
  NodeID curr_node = nn;
  NodeID prev_node = n;

  double graph_len = 0.0;
  int leaf_flag = 0;

  while(true) {
      auto curr_nbs = next_neighbours(g, prev_node, curr_node);
      if(curr_nbs.size() > 1)
          break;
      else if(curr_nbs.size() == 0) {
          branch_path.push_back(curr_node); graph_len += abs(length(g.pos[curr_node] - g.pos[prev_node]));leaf_flag = 1; break;
      }
      else {
          branch_path.push_back(curr_node);
          graph_len += abs(length(g.pos[curr_node] - g.pos[prev_node]));
          prev_node = curr_node;
          curr_node = curr_nbs[0];
      }
  }

  return graph_len;

}

//Mesh util functions

void id_preserving_cc(HMesh::Manifold& m_in) {

    vector<FaceID> base_faces;
    int Invalid = -1;
    HalfEdgeAttributeVector<int> htouched(m_in.allocated_halfedges(), Invalid);
    int Valid = 1;
    map<FaceID,VertexID> face2centerv;
    vector<HalfEdgeID> base_edges;

    vector<HalfEdgeID> new_edges;


    for(auto f: m_in.faces())
        base_faces.push_back(f);

    for(auto h: m_in.halfedges())
        base_edges.push_back(h);

    for(auto f: base_faces)
        if(m_in.in_use(f)) {
            VertexID center_v = m_in.split_face_by_vertex(f);
            for(Walker w = m_in.walker(center_v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                new_edges.push_back(w.halfedge());
                face2centerv.insert(std::make_pair(w.face(), center_v));
            }
        }

    FaceAttributeVector<int> ftouched(m_in.allocated_faces(), Invalid);

    for (auto h: base_edges) {
        FaceID f1 = m_in.walker(h).face();
        FaceID f2 = m_in.walker(h).opp().face();
        if(ftouched[f1]==Invalid && ftouched[f2]==Invalid) {
            VertexID opp_v = m_in.split_edge(h);
            ftouched[f1] = Valid;
            ftouched[f2] = Valid;
            m_in.split_face_by_edge(f1, face2centerv.find(f1)->second, opp_v);
            m_in.split_face_by_edge(f2, face2centerv.find(f2)->second, opp_v);
         }

    }
    for (auto h_dissolve: new_edges)
        if(m_in.in_use(h_dissolve))
            m_in.merge_faces(m_in.walker(h_dissolve).face(), h_dissolve);
    return;

}

void quad_mesh_leaves(HMesh::Manifold& m) {

    vector<FaceID> base_faces;
    vector<HalfEdgeID> new_edges;

    for(auto f: m.faces())
        if(no_edges(m, f) != 4  || (val2_faces.count(f) && one_ring_face_vertex.count(f))) {
             HalfEdgeID ref_h;

             VertexID ref_v  = one_ring_face_vertex[f];

             if(ref_v == InvalidVertexID)
                 continue;

            VertexID center_v = m.split_face_by_vertex(f);
            int counter = 0;
            int dissolve_flag = 0;
            for(Walker w = m.walker(center_v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                if(m.walker(w.halfedge()).vertex() == ref_v || m.walker(w.halfedge()).opp().vertex() == ref_v) {
                    dissolve_flag = counter%2;
                    ref_h = w.halfedge();
                }
                counter++;
            }
            counter = 0;
            for(Walker w = m.walker(center_v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                if(counter%2 != dissolve_flag)
                    new_edges.push_back(w.halfedge());
                counter++;
            }

        }


    for (auto h_dissolve: new_edges)
        if(m.in_use(h_dissolve))
            m.merge_faces(m.walker(h_dissolve).face(), h_dissolve);


    return;


}

VertexID split_LIE(Manifold& mani, HalfEdgeID h) {
    Walker w = mani.walker(h);
    VertexID v = w.next().vertex();
    VertexID vo = w.opp().next().vertex();
    FaceID f = w.face();
    FaceID fo = w.opp().face();

    if(v == vo) {
        return InvalidVertexID;
    }
    VertexID vid = mani.split_edge(h);

    mani.split_face_by_edge(f, vid, v);
    mani.split_face_by_edge(fo, vid, vo);
    return vid;
}

bool check_planar(const HMesh::Manifold &m, HalfEdgeID h) {

    FaceID face_1 = m.walker(h).face();
    FaceID face_2 = m.walker(h).opp().face();

    Vec3d normal_1 = normal(m, face_1);
    Vec3d normal_2 = normal(m, face_2);
    Vec3d center_1 = centre(m, face_1);
    Vec3d center_2 = centre(m, face_2);

    float dot_val_1 = dot(normal_1, center_2 - center_1);
    float dot_val_2 = dot(normal_2, center_1 - center_2);

    if(dot_val_1 > 0 && dot_val_2 > 0)
        return true;

    float dot_val = dot(normal_1 , normal_2);

    if(dot_val < 0)
        return false;

    dot_val = abs(dot_val);

    if(dot_val < 0.75)
        return false;
    else
        return true;
}

bool check_convex(const HMesh::Manifold &m, HalfEdgeID h) {

    auto v1 = m.walker(h).vertex();
    auto v2 = m.walker(h).opp().vertex();

    auto v3 = m.walker(h).next().vertex();
    auto v4 = m.walker(h).opp().prev().prev().vertex();

    Vec3d edge_vec_1 = m.pos(v3) - m.pos(v1);
    Vec3d edge_vec_2 = m.pos(v1) - m.pos(v4);

    Vec3d edge_vec_3 = m.pos(v3) - m.pos(v2);
    Vec3d edge_vec_4 = m.pos(v2) - m.pos(v4);

    Vec3d diag_vec = m.pos(v2) - m.pos(v1);

    float angle_1 = acos(dot(normalize(diag_vec),normalize(-edge_vec_3))) + acos(dot(normalize(diag_vec),normalize(edge_vec_4)));

    float angle_2 = acos(dot(normalize(diag_vec),normalize(edge_vec_1))) + acos(dot(normalize(diag_vec),normalize(-edge_vec_2)));

    angle_1 *= 180/3.14;
    angle_2 *= 180/3.14;

    if (angle_1 > 180 || angle_2 > 180)
        return false;
    else
        return true;
}

void stellate_face_set_retopo(HMesh::Manifold &m, const HMesh::FaceSet& fs, const HMesh::HalfEdgeSet& hs) {

  HalfEdgeSet interior_edges;

  HalfEdgeSet visited;

  VertexSet aux_vertices;

  for (auto h : hs) {
    if(visited.find(h) != visited.end() || visited.find(m.walker(h).opp().halfedge()) != visited.end())
      continue;

    VertexID interior_vertex = m.walker(h).next().vertex();
    VertexID interior_vertex_opp = m.walker(h).prev().prev().vertex();

    if(aux_vertices.find(interior_vertex) != aux_vertices.end() || aux_vertices.find(interior_vertex_opp) != aux_vertices.end()) {
      m.flip_edge(h);
      visited.insert(h);
      visited.insert(m.walker(h).opp().halfedge());
      continue;
    }


    visited.insert(h);
    visited.insert(m.walker(h).opp().halfedge());
    VertexID curr_v = split_LIE(m, h);
    aux_vertices.insert(curr_v);
    for (auto h_flip : hs) {
      if(visited.find(h_flip) == visited.end() && (m.walker(h_flip).next().vertex() == curr_v || m.walker(h_flip).prev().prev().vertex() == curr_v) && check_planar(m,h_flip)) {
        m.flip_edge(h_flip);
        visited.insert(h_flip);
        visited.insert(m.walker(h_flip).opp().halfedge());
      }
    }

  }
  for (auto h_flip : m.halfedges()) {
    VertexID interior_vertex = m.walker(h_flip).next().vertex();
    if(visited.find(h_flip) == visited.end() && aux_vertices.find(interior_vertex) != aux_vertices.end() && check_planar(m, h_flip)) {
      m.flip_edge(h_flip);
      visited.insert(h_flip);
      visited.insert(m.walker(h_flip).opp().halfedge());
    }
  }

}

vector<FaceSet> retopologize_planar_regions(HMesh::Manifold &m) {

// find triangles with low dihedral angle between them

    vector<FaceSet> planar_regions;

    vector<HalfEdgeSet> planar_edges;

    // identify planar regions

    FaceSet global_visited;

    for (auto f : m.faces()) {

        if(global_visited.find(f) != global_visited.end())
            continue;

        FaceSet planar_set;
        HalfEdgeSet planar_edge_set;

        FaceSet visited;

        queue<FaceID> Q;

        Q.push(f);

        while(!Q.empty()) {

            FaceID curr_f = Q.front();

            HalfEdgeSet edge_set;


            for(Walker w = m.walker(curr_f); !w.full_circle(); w = w.circulate_face_ccw()) {
                edge_set.insert(w.halfedge());
            }

            for (auto h : edge_set)
                if(check_planar(m,h)) {
                  planar_edge_set.insert(h);
                  planar_edge_set.insert(m.walker(h).opp().halfedge());
                  FaceID f1 = m.walker(h).face();
                  FaceID f2 = m.walker(h).opp().face();

                  planar_set.insert(f1);
                  planar_set.insert(f2);

                  global_visited.insert(f1);
                  global_visited.insert(f2);

                  if(visited.find(f1) == visited.end())
                    Q.push(f1);
                  if(visited.find(f2) == visited.end())
                    Q.push(f2);
                }

            visited.insert(curr_f);
            Q.pop();
       }

       planar_regions.push_back(planar_set);
       planar_edges.push_back(planar_edge_set);
    }
    int id = 0;
    for (auto fset : planar_regions) {
      if(fset.size() > 2)
        stellate_face_set_retopo(m, fset, planar_edges[id]);
      id++;
    }

    return planar_regions;

}

//Graph - Mesh relationship Functions

VertexID branch2vertex (HMesh::Manifold &m_out, const Geometry::AMGraph3D& g,
                        NodeID n, NodeID nn, const Util::AttribVec<NodeID, FaceSet>& node2fs) {

    Vec3d vert_pos = branch2vert.find(std::make_pair(n,nn))->second;

    for (auto v: m_out.vertices())
        if(sqr_length(m_out.pos(v) - vert_pos) == 0)
            return v;

    return InvalidVertexID;

}

void init_branch_degree(HMesh::Manifold &m, const Geometry::AMGraph3D& g,
                        const Util::AttribVec<NodeID, FaceSet>& node2fs) {


    for (auto n:g.node_ids()) {
        auto N = g.neighbors(n);

        //for all branch nodes

        if(N.size() > 2) {


        // for each outgoing arc

            for (auto nn: N) {

                int src_branch_degree = valency(m, branch2vertex(m, g, n,nn, node2fs));
                vector<NodeID> branch_path;
                NodeID curr_node = nn;
                NodeID prev_node = n;

                int leaf_flag = 0;

              //traverse val 2 nodes to next branch node


                while(true) {
                    auto curr_nbs = next_neighbours(g, prev_node, curr_node);
                    if(curr_nbs.size() > 1)
                        break;
                    else if(curr_nbs.size() == 0) {
                        branch_path.push_back(curr_node); leaf_flag = 1; break;
                    }
                    else {
                        branch_path.push_back(curr_node);
                        prev_node = curr_node;
                        curr_node = curr_nbs[0];
                    }
                }

                //pick lower degree

                int dest_branch_degree;
                if(leaf_flag == 1) {
                    dest_branch_degree = src_branch_degree;
                }
                else
                    dest_branch_degree = valency(m, branch2vertex(m, g, curr_node, prev_node, node2fs));


                int path_degree = 0;
                int jn_degree = 0;

                if(dest_branch_degree < src_branch_degree) {
                    path_degree = (dest_branch_degree)*2;
                    jn_degree = dest_branch_degree - 1;
                }

                else if (dest_branch_degree == src_branch_degree) {
                    path_degree = dest_branch_degree*2;
                    jn_degree = dest_branch_degree;
                }

                else {
                    jn_degree = src_branch_degree - 1;
                    path_degree = (src_branch_degree)*2;
                }

                auto key = std::make_pair(n,nn);
                branchdeg.insert(std::make_pair(key,jn_degree));
                for (auto val2node : branch_path)
                    val2deg.insert(std::make_pair(val2node,path_degree));

            }
        }
    }

    // for junction-less graphs

    bool has_junction = false;

    for (auto n: g.node_ids()) {
      if(g.valence(n) > 2)
        has_junction = true;
    }

    if(!has_junction)
      for (auto n : g.node_ids())
        if(g.valence(n) <= 2)
          if(val2deg.find(n) == val2deg.end())
            val2deg.insert(std::make_pair(n,4));

}

FaceID branch2face (HMesh::Manifold &m_out,
                    const Geometry::AMGraph3D& g, NodeID n, NodeID nn,
                    Util::AttribVec<NodeID, FaceSet>& node2fs) {

    VertexID v = branch2vertex(m_out, g, n, nn, node2fs);
    vector<FaceID> face_set;

    double d_max = FLT_MAX;
    FaceID f_max = InvalidFaceID;
    Vec3d pn = g.pos[n];
    Vec3d pnn = g.pos[nn];
    Vec3d v_n_nn = pnn - pn;


    for(Walker w = m_out.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
        face_set.push_back(w.face());
    }

    for(auto f: face_set) {
        double d = dot(v_n_nn, normal(m_out, f));

        Vec3d face_normal = normal(m_out, f);
        Vec3d face_center = centre(m_out,f);
        float face_plane_d = dot(face_normal, face_center);

        float intersection_x = (face_plane_d - dot(face_normal, pn)) / dot(face_normal, v_n_nn);

        Vec3d intersection_pt = pn + intersection_x*v_n_nn;

        d = sqr_length(face_center - intersection_pt);

        if(d < d_max) {
            f_max = f;
            d_max = d;
        }
    }
    if(g.neighbors(n).size()>2)
        node2fs[n].erase(f_max);
    return f_max;


}

double face_dist(HMesh::Manifold &m_out, const Geometry::AMGraph3D& g, NodeID n, NodeID nn, FaceID f) {
    Vec3d pn = g.pos[n];
    Vec3d pnn = g.pos[nn];
    Vec3d v_n_nn = pnn - pn;


    if(!m_out.in_use(f))
        return 1;

    Vec3d face_normal = normal(m_out, f);
    Vec3d face_center = centre(m_out,f);
    float face_plane_d = dot(face_normal, face_center);

    float intersection_x = 0;
    if(dot(face_normal, v_n_nn) != 0)
       intersection_x = (face_plane_d - dot(face_normal, pn)) / dot(face_normal, v_n_nn);

    Vec3d intersection_pt = pn + intersection_x*v_n_nn;

    float d = sqr_length(face_center - intersection_pt);

    return d;

}

void init_branch_face_pairs(HMesh::Manifold &m, const Geometry::AMGraph3D& g,
                            Util::AttribVec<NodeID, FaceSet>& node2fs) {

    for (auto n:g.node_ids()) {

        auto N = g.neighbors(n);

        //for all branch nodes

        if(N.size() > 2) {

            // for each outgoing arc

            for (auto nn: N) {

                auto key = std::make_pair(n,nn);

                FaceID f = branch2face(m, g, n, nn, node2fs);
                branch_best_face.insert(std::make_pair(key,f));

                VertexID v = branch2vertex(m, g, n, nn, node2fs);
                branch_best_vertex.insert(std::make_pair(key,v));

                one_ring_vertex.insert(std::make_pair(key,InvalidVertexID));

            }
        }

    }
}

double graph_unique_score(HMesh::Manifold &m, const Geometry::AMGraph3D& g, NodeID n, NodeID nn, VertexID curr_vert,
                          const Util::AttribVec<NodeID, FaceSet>& node2fs) {
  auto curr_nbs = next_neighbours(g, nn, n);
  vector<double> graph_lens;
  for (auto nb : curr_nbs) {
    double curr_len = graph_length(g, n, nb);
    graph_lens.push_back(curr_len);
  }
  double max_diff = FLT_MAX;
  for (auto nb : curr_nbs) {
    if(branch2vertex(m, g, n, nb, node2fs) == curr_vert) {
        double curr_len = graph_length(g, n, nb);
        for (auto graph_len : graph_lens) {
          double curr_diff = abs(curr_len - graph_len);
          if(curr_diff < max_diff && curr_diff != 0)
              max_diff = curr_diff;
        }
    return max_diff;
  }
}
return -1;
}

std::map<VertexID,int> val_diff_map(HMesh::Manifold &m, const Geometry::AMGraph3D& g) {

    Util::AttribVec<NodeID, FaceSet> node2fs;
    map<VertexID,int> vertex2val;
    for (auto n:g.node_ids()) {
        auto N = g.neighbors(n);

            //for all branch nodes

        if(N.size() > 2) {


            // for each outgoing arc

            for (auto nn: N) {

                VertexID src_vertex = branch2vertex(m, g, n,nn, node2fs);
                int src_branch_degree = valency(m, src_vertex);
                vector<NodeID> branch_path;
                NodeID curr_node = nn;
                NodeID prev_node = n;

                int leaf_flag = 0;

                //traverse val 2 nodes route to next branch node


                while(true) {
                    auto curr_nbs = next_neighbours(g, prev_node, curr_node);
                    if(curr_nbs.size() > 1)
                        break;
                    else if(curr_nbs.size() == 0) {
                        branch_path.push_back(curr_node); leaf_flag = 1; break;
                    }
                    else {
                        branch_path.push_back(curr_node);
                        prev_node = curr_node;
                        curr_node = curr_nbs[0];
                    }
                }

                int dest_branch_degree;
                if(leaf_flag == 1) {
                    dest_branch_degree = 0;
                    vertex2val.insert(std::make_pair(src_vertex, 0));
                }
                else {
                    dest_branch_degree = valency(m, branch2vertex(m, g, curr_node, prev_node, node2fs));
                    vertex2val.insert(std::make_pair(src_vertex, dest_branch_degree - src_branch_degree));
                }
            }
        }
    }
    return vertex2val;
}

//Functions for constructing / editing mesh elements from skeletal nodes

vector<Vec3d> get_face_points(int n) {

    vector<Vec3d> face_vertices;
    double h = 0.5;
    double angle = 0;

    for(int i =0; i<n; i++) {

        Vec3d face_vertex = Vec3d(0, h*cos(angle), h*sin(angle));
        face_vertices.push_back(face_vertex);
        angle+=2*22.0/(7.0*n);
    }
    return face_vertices;
}

vector<FaceID> create_face_pair(Manifold& m, const Vec3d& pos, const Mat3x3d& _R, int axis, int num_sides) {
    if(num_sides == 0) {
        vector<FaceID>fvec;
        return fvec;
    }

    vector<Vec3d> face_points = get_face_points(num_sides);
    Mat3x3d R = _R;
    double det = determinant(R);
    if(abs(det))
        if(det<0) {
            Mat3x3d M = identity_Mat3x3d();
            M[2][2] = -1;
            R = R * M;
        }

    vector<FaceID> fvec;
        vector<Vec3d> front_pts;
        for(int i = 0; i < num_sides; i++)
        {
            Vec3d _p = face_points[i];
            Vec3d p(0);
            p[(0+axis)%3] += _p[0];
            p[(1+axis)%3] += _p[1];
            p[(2+axis)%3] += _p[2];
            front_pts.push_back(R*p+pos);
        }
        fvec.push_back(m.add_face(front_pts));

        vector<Vec3d> back_pts;
        for (int i = 0; i < num_sides; i++) {
            int curr_index = 1 - i;
            if(curr_index < 0 )
                curr_index = num_sides + curr_index;
            Vec3d _p = face_points[curr_index];
            Vec3d p(0);
            p[(0+axis)%3] += _p[0];
            p[(1+axis)%3] += _p[1];
            p[(2+axis)%3] += _p[2];
            back_pts.push_back(R*p+pos);
        }


        fvec.push_back(m.add_face(back_pts));

    return fvec;
}

void val2nodes_to_boxes(const Geometry::AMGraph3D& g, HMesh::Manifold& mani,
                        Util::AttribVec<NodeID, FaceSet>& n2fs,
                        const vector<double>& r) {
    Vec3d c(0);
    for(auto n: g.node_ids())
        c += g.pos[n];
    c /= g.no_nodes();
    double min_dist=DBL_MAX;
    NodeID middle_node = *begin(g.node_ids());
    for(auto n: g.node_ids())
        if(g.valence(n)>2)
        {
            double d = sqr_length(g.pos[n]-c);
            if(d < min_dist) {
                min_dist = d;
                middle_node = n;
            }
        }
    if(g.valence(middle_node) == 0)
      for (auto n : g.node_ids())
        if(g.valence(n) > 0)
        {
          double d = sqr_length(g.pos[n] - c);
          if(d < min_dist) {
            min_dist = d;
            middle_node = n;
          }
        }
    Util::AttribVec<NodeID, int> touched(g.no_nodes(),0);
    Util::AttribVec<NodeID, Mat3x3d> warp_frame(g.no_nodes(),identity_Mat3x3d());

    queue<NodeID> Q;
    Q.push(middle_node);

    if(g.valence(middle_node) > 2)
      touched[middle_node] = 1;

    while(!Q.empty()) {
        NodeID n = Q.front();
        Q.pop();
        for(auto m : g.neighbors(n))
            if (!touched[m]) {
                Q.push(m);
                touched[m] = 1;
                Vec3d v = g.pos[m]-g.pos[n];
                if(g.neighbors(m).size() <= 2) {
                  auto node_list = next_neighbours(g, n, m);
                  Vec3d nb_v(0);
                  for(auto m_nb : node_list)
                    nb_v = g.pos[m_nb];
                  if(nb_v != Vec3d(0)) {
                    v = 0.5*(v + (nb_v - g.pos[m]));
                  }
                }
                Mat3x3d M = warp_frame[n];
                Vec3d warp_v = M * v;

                double max_sgn = sign(warp_v[0]);
                double max_val = abs(warp_v[0]);
                int max_idx = 0;
                for(int i=1;i<3;++i) {
                    if(abs(warp_v[i])>max_val) {
                        max_sgn = sign(warp_v[i]);
                        max_val = abs(warp_v[i]);
                        max_idx = i;
                    }
                }
                auto v_target = max_sgn * normalize(v);
                Quatd q;
                q.make_rot(M[max_idx], v_target);
                M = transpose(q.get_Mat3x3d() * transpose(M));
                warp_frame[m] = M;

                if(g.neighbors(m).size()<=2) {
                    Vec3d s(r[m]);
                    Mat3x3d S = scaling_Mat3x3d(s);
                    auto face_list = create_face_pair(mani, g.pos[m], transpose(M)*S, max_idx, val2deg.find(m)->second);
                    stitch_mesh(mani, 1e-10);
                    for(auto f: face_list) {
                            n2fs[m].insert(f);
                            val2_faces.insert(std::make_pair(f, 1));
                    }
                }


            }
    }

}

int add_ghosts(const vector<Vec3i>& tris, vector<Vec3d>& pts) {
    
    /* This function creates extra points to add to the BNP vertices for a branch node.
     These extra points are called ghost points because they do not correspond to an outgoing
     edge.
     
     To do this, we categorize nodes according to three categories
     
     Type A (flat): All outgoing edges lie more or less in a plane for this type of nodes, and
     the BNP becomes very flat.
     
     Type B (semi flat): All outgoing edges again lie more or less in a plane, but on one side of
     the BNP one or more additional outgoing edges emanate.
     
     Type C (general): This type comprises everything else. For type C nodes, it is hard to say
     something meaningful about the configuration of outgoing edges.
     
     The idea behind this function is to add two ghost points to Type A nodes, a single for Type B
     and 0 for Type C. This heuristic is based on the observation that for Type A this will lead to
     valence 4 vertices generally since most of the vertices will be connected to two neighbouring
     points and the two ghosts. For Type B, the same is true if there is a single outgoing edge
     perpendicular to the edges that form a flat region.
     
     Procedure:
     Initially, we add a ghost point for every triangle in the
     initial BNP where the smallest pairwise dot product between the vertices
     is less than -0.1. The ghost point is the normal of the triangle.
     Intuitively, this adds what corresponds to a new outgoing edge in a direction
     of the sphere that is otherwise not well covered.
     */
    vector<Vec3d> ghost_pts;
    for(auto t: tris) {
        const auto& p0 = pts[t[0]];
        const auto& p1 = pts[t[1]];
        const auto& p2 = pts[t[2]];
        Vec3d v1 = pts[t[1]]-pts[t[0]];
        Vec3d v2 = pts[t[2]]-pts[t[0]];
        double l = min(dot(p0, p1), min( dot(p1,p2), dot(p2,p0)));
        if (l<-0.1) {
            ghost_pts.push_back(normalize(cross(v1, v2)));
        }
    }
    
    /* Next, we cluster the ghost points. This is because in flatish
     configurations we could have several quite similar ghost points.
     The threshold -0.1 is very loose to avoid adding too many ghosts */
    vector<int> cluster_id(ghost_pts.size(), -1);
    int max_id = 0;
    for(int i=0;i<ghost_pts.size(); ++i) {
        if (cluster_id[i] == -1) {
            cluster_id[i] = max_id++;
        }
        for(int j=i+1; j<ghost_pts.size(); ++j) {
            if (cluster_id[j] == -1) {
                if (dot(ghost_pts[i], ghost_pts[j]) > -0.1)
                    cluster_id[j] = cluster_id[i];
            }
        }
    }
    
    vector<Vec3d> ghost_pts_new(max_id, Vec3d(0));
    for(int i=0;i<ghost_pts.size(); ++i) {
        ghost_pts_new[cluster_id[i]] += ghost_pts[i];
    }
      
    /* Finally, we cull ghost points too close to an existing non-ghost point.
     The threshold of 0.4 allows ghost points to be a little closer to other
     points than each other. */
    ghost_pts.resize(0);
    for(auto& p: ghost_pts_new) {
        p.normalize();
        vector<double> dots;
        for(const auto& p_orig: pts)
            dots.push_back(dot(p,p_orig));
        if(*max_element(begin(dots), end(dots))<0.4)
            ghost_pts.push_back(p);
    }
    
    /* If there are more than three ghost points, it is a Type C BNP, and
     we do nothing. */
    if(ghost_pts.size()>2)
        return 0;
        
    for (auto g: ghost_pts)
        pts.push_back(g);
    
    return ghost_pts.size();
}

void construct_bnps(HMesh::Manifold &m_out, const Geometry::AMGraph3D& g, Util::AttribVec<NodeID, FaceSet>& node2fs, vector<double> r_arr) {

  map<int, pair<NodeID,NodeID>> spts2branch;
  map <int, VertexID> spts2vertexid;

  double base_r = 0.5 * g.average_edge_length();

  double r = base_r;

  for (auto n: g.node_ids()) {

          auto N = g.neighbors(n);
          if(N.size()>2) {
              Manifold m;
              int node_vertex_count =0;
              Vec3d pn = g.pos[n];

              r = r_arr[n];

              auto project_to_sphere = [&]() {
                  for(int iter=0;iter<1;++iter) {
                      auto new_pos = m.positions_attribute_vector();
                      for(auto v: m.vertices())
                          new_pos[v] = normalize(1*normal(m,v) + m.pos(v));
                      m.positions_attribute_vector() = new_pos;
                  }
                  for(auto v: m.vertices())
                      m.pos(v) = normalize(m.pos(v))*r + pn;
              };
              vector<Vec3d> spts;
              int spts_vertex_count = 0;

              spts2branch.clear();
              spts2vertexid.clear();

              for (auto nn: N) {
                  Vec3d pnn = g.pos[nn];
                  spts.push_back(normalize(pnn-pn));

                  auto spts_value = std::make_pair(n,nn);
                  auto spts_key = spts_vertex_count;
                  spts2branch.insert(std::make_pair(spts_key,spts_value));
                  spts_vertex_count++;

              }
              bool ghost_added = false;

              std::vector<CGLA::Vec3i> stris = SphereDelaunay(spts);

              if (!ghost_added)
                  if (add_ghosts(stris, spts)>0) {
                      stris = SphereDelaunay(spts);
                  }

              for(auto tri: stris) {
                  vector<Vec3d> triangle_pts;
                  for(int i=0;i<3; ++i) {
                      triangle_pts.push_back(spts[tri[i]]);
                      node_vertex_count++;

                  }
                  m.add_face(triangle_pts);
              }
              stitch_mesh(m, 1e-10);

              m.cleanup();

              for(auto v: m.vertices())
                  for(int i = 0; i < spts.size(); i++)
                      if(sqr_length(m.pos(v) - spts[i]) < 0.0001)
                          spts2vertexid.insert(std::make_pair(i, v));


              if(N.size() > 3 && !ghost_added)
                vector<FaceSet> planar_regions = retopologize_planar_regions(m);

              taubin_smooth(m, 1);

              project_to_sphere();

              for(int iter=0;iter<10;++iter) {
                  auto new_pos = m.positions_attribute_vector();
                  for(auto v: m.vertices())
                      new_pos[v] = 0.01*r*normal(m,v) + m.pos(v);
                  m.positions_attribute_vector() = new_pos;
              }

              for(int i = 0; i < spts.size(); i++) {
                  auto key = spts2branch.find(i)->second;
                  auto value = m.pos(spts2vertexid.find(i)->second);
                  branch2vert.insert(std::make_pair(key,value));
              }

              m.cleanup();

              size_t no_faces_before = m_out.no_faces();



              m_out.merge(m);
              for(auto f: m_out.faces())
                  if(f.index >= no_faces_before)
                      node2fs[n].insert(f);
          }
    }
    m_out.cleanup();
    stitch_mesh(m_out, 1e-10);
}

void refine_BNPs(HMesh::Manifold &m, const Geometry::AMGraph3D& g,
                 const Util::AttribVec<NodeID, FaceSet>& node2fs) {

    map<VertexID, int> vertex2valdiff = val_diff_map(m,g);
    bool work_done = false;
    int count = 0;
    VertexSet new_vertices;
    VertexSet end_vertices;
    do {

        count++;

        if(count > 1)
            return;

        work_done = false;
        for (auto n:g.node_ids()) {
            auto N = g.neighbors(n);

            //for all branch nodes

            if(N.size() > 2) {

                for (auto nn: N) {

                    VertexID src_vertex = branch2vertex(m, g, n, nn, node2fs);

                    int refine_number = vertex2valdiff.find(src_vertex)->second;

                    HalfEdgeID h_split = InvalidHalfEdgeID;

                    VertexID end_vertex_split = InvalidVertexID;

                    if(refine_number > 1 || refine_number < 1)
                        continue;

                    for (int i = 0; i < refine_number; i++) {
                        double best_split_score = FLT_MAX;
                        h_split = InvalidHalfEdgeID;

                        circulate_vertex_ccw(m, branch2vertex(m,g,n,nn,node2fs), std::function<void(HalfEdgeID)>([&](HalfEdgeID h) {
                            double curr_split_score;
                            HalfEdgeID curr_link_edge;
                            if(m.walker(h).vertex() == branch2vertex(m,g,n,nn,node2fs))
                                curr_link_edge = m.walker(h).prev().halfedge();
                            else if(m.walker(h).opp().vertex() == branch2vertex(m,g,n,nn,node2fs))
                                curr_link_edge = m.walker(h).next().halfedge();
                            VertexID new_vertex = m.walker(curr_link_edge).opp().next().vertex();
                            if(vertex2valdiff.find(new_vertex) != vertex2valdiff.end()) {
                                if(vertex2valdiff.find(new_vertex)->second == 0) {
                                    double unique_score = graph_unique_score(m, g, n, nn, new_vertex, node2fs);
                                    curr_split_score = refine_number + 1.0 / (1 + abs(unique_score));
                                }
                                else
                                    curr_split_score = abs(refine_number - vertex2valdiff.find(new_vertex)->second);
                                if(new_vertices.find(m.walker(curr_link_edge).vertex()) != new_vertices.end() || new_vertices.find(m.walker(curr_link_edge).opp().vertex()) != new_vertices.end())
                                    curr_split_score += 1;
                                if(end_vertices.find(new_vertex) != end_vertices.end()) {
                                    curr_split_score += 1;
                                }
                                if(curr_split_score < best_split_score) {
                                    best_split_score = curr_split_score;
                                    h_split = curr_link_edge;
                                    end_vertex_split = new_vertex;
                                }
                            }
                        }));
                        if(h_split != InvalidHalfEdgeID) {
                            VertexID split_vertex = split_LIE(m, h_split);
                            new_vertices.insert(split_vertex);
                            end_vertices.insert(end_vertex_split);
                            if(split_vertex != InvalidVertexID) {
                                work_done = true;
                                vertex2valdiff = val_diff_map(m,g);
                            }
                        }
                    }
                }
            }
        }

    } while (work_done);
    return;
}

void merge_branch_faces(HMesh::Manifold &m, const Geometry::AMGraph3D& g,
                        const Util::AttribVec<NodeID, FaceSet>& node2fs) {

    VertexID v;
    FaceID f, face_1, face_2;
    int branch_degree;
    HalfEdgeID boundary_edge_1;
    HalfEdgeID boundary_edge_2;


    for (auto n:g.node_ids()) {
        auto N = g.neighbors(n);

        //for all branch nodes

        if(N.size() > 2) {

            // for each outgoing arc

            for (auto nn: N) {

                auto key = std::make_pair(n,nn);

                branch_degree = branchdeg.find(key)->second;

                f = branch_best_face.find(key)->second;

                v = branch_best_vertex.find(key)->second;

                if(valency(m,v) == branch_degree) {

                    HalfEdgeID ref_he;
                    VertexID ref_v;

                    for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw()) {
                        ref_he = w.halfedge();
                        if(m.walker(ref_he).vertex() == v)
                            ref_v = m.walker(ref_he).opp().vertex();
                        else if (m.walker(ref_he).opp().vertex() == v)
                            ref_v = m.walker(ref_he).vertex();
                    }

                    FaceID f = m.merge_one_ring(v);

                    if(m.in_use(ref_v)) {
                        one_ring_vertex[key] = ref_v;
                        one_ring_face_vertex[f] = ref_v;
                    }
                    else {
                        one_ring_vertex[key] = InvalidVertexID;
                        one_ring_face_vertex[f] = InvalidVertexID;
                    }

                    branchface.insert(std::make_pair(key,f));
                    branch_best_vertex[key] = InvalidVertexID;
                    continue;

                }

                branchface.insert(std::make_pair(key,f));

                for(int i = 0; i < branch_degree - 1; i++) {

                    Walker w_f = m.walker(f);
                    Walker w_h = m.walker(w_f.halfedge());


                    HalfEdgeID w_start = w_f.halfedge();

                    do {

                        if(w_h.vertex() == v) {
                            boundary_edge_1 = w_h.halfedge();
                        }
                        if(w_h.opp().vertex() == v) {
                            boundary_edge_2 = w_h.halfedge();
                        }
                        w_h=w_h.next();

                    } while(w_h.halfedge() != w_start);

                    Walker b_e_1 = m.walker(boundary_edge_1);
                    Walker b_e_2 = m.walker(boundary_edge_2);
                    face_1 = b_e_1.opp().face();
                    face_2 = b_e_2.opp().face();


                    if( ! m.in_use(face_1)) {
                        m.merge_faces(f, boundary_edge_2);
                        continue;
                    }
                    if( ! m.in_use(face_2)) {
                        m.merge_faces(f, boundary_edge_1);
                        continue;
                    }


                    if(face_dist(m, g, n, nn, face_1) < face_dist(m, g, n, nn, face_2))
                        m.merge_faces(f, boundary_edge_1);
                    else
                        m.merge_faces(f, boundary_edge_2);

                    }
              }
        }

    }

    return;

}

//Bridging Functions

FaceID rotate_bridge_face_set_once(HMesh::Manifold& m, FaceID f0) {
    VertexID central_vertex_0 = face_vertex[f0];
    HalfEdgeID bd_edge;

    vector<VertexID> vloop0;
    circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
        vloop0.push_back(v);
    }) );


    circulate_face_ccw(m, f0, std::function<void(HalfEdgeID)>([&](HalfEdgeID h){
        if(m.walker(h).vertex() == central_vertex_0)
            bd_edge = h;
    }));



    size_t L= vloop0.size();
    VertexID split_vertex;

    for(int i = 0; i < L; i++) {
        if(vloop0[i] == central_vertex_0)
            split_vertex = vloop0[(i+3)%L];
    }

    float dot_sum = 0;
    float len = 0;

    FaceID new_face = f0;
    if(L == 4) {
        if(m.walker(bd_edge).face() == f0)
            new_face = m.walker(bd_edge).opp().face();
        else
            new_face = m.walker(bd_edge).face();
    }
    else {
       new_face = m.split_face_by_edge(f0, central_vertex_0, split_vertex);
       if(m.in_use(bd_edge))
           m.merge_faces(m.walker(bd_edge).face(),bd_edge);
      }

    f0 = new_face;
    face_vertex[new_face] = central_vertex_0;
    return f0;
}

double compute_torsion_score(HMesh::Manifold& m, FaceID f0, FaceID f1) {

  VertexID central_vertex_0 = face_vertex[f0];
  VertexID central_vertex_1 = face_vertex[f1];
  vector<VertexID> vloop0;
  circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
      vloop0.push_back(v);
  }) );
  vector<VertexID> vloop1;
  circulate_face_ccw(m, f1, std::function<void(VertexID)>([&](VertexID v){
      vloop1.push_back(v);
  }) );
  float len = 0;
  float dot_sum = FLT_MAX;
  size_t L= vloop0.size();

  if(central_vertex_0 != InvalidVertexID && central_vertex_1 != InvalidVertexID) {

    for(int j_off = 0; j_off < L; ++j_off) {
       bool center_match = false;
       for(int i=0;i<L;++i)
          if(vloop0[i] == central_vertex_0 && vloop1[(L + j_off - i)%L] == central_vertex_1)
              center_match = true;
       if(center_match) {
         len = 0;
          for(int i=0;i<L;i++) {
              Vec3d bridge_edge_i, bridge_edge_j;
              len += abs(length(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L])));
              bridge_edge_i = normalize(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L]));
              for(int j = 0; j < L; j++) {
                bridge_edge_j = normalize(m.pos(vloop0[j]) - m.pos(vloop1[(L+j_off - j)%L]));
                double curr_dot_sum = abs(dot(bridge_edge_i, bridge_edge_j));
                if(curr_dot_sum < dot_sum)
                  dot_sum = curr_dot_sum;
              }
          }
       }
     }
     return dot_sum;
   }
   else {
     float min_len = FLT_MAX;
     int j_off_min_len = -1;
     for(int j_off = 0; j_off < L; j_off = j_off + 1) {
         len = 0;
         for(int i=0;i<L;++i) {
             Vec3d bridge_edge_i, bridge_edge_j;
             bridge_edge_i = normalize(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L]));
             for(int j = 0; j < L; j++) {
               bridge_edge_j = normalize(m.pos(vloop0[j]) - m.pos(vloop1[(L+j_off - j)%L]));
               double curr_dot_sum = abs(dot(bridge_edge_i, bridge_edge_j));
               if(curr_dot_sum < dot_sum)
                 dot_sum = curr_dot_sum;
             }
             len += sqr_length(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L]));
         }
         if(len < min_len)   {
             j_off_min_len = j_off;
             min_len = len;
         }
     }
     return dot_sum;
   }
   return len;
}

FaceID rotate_bridge_face_set(HMesh::Manifold& m, FaceID f0, FaceID f1,
                              const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {
    VertexID central_vertex_0 = face_vertex[f0];
    VertexID central_vertex_1 = face_vertex[f1];
    int central_valency = valency(m, central_vertex_0);
    Vec3d cv_edge = normalize(m.pos(central_vertex_0) - m.pos(central_vertex_1));
    float max_dot_sum = -FLT_MAX;
    double min_len = FLT_MAX;

    for(int iter = 0; iter <= 2*central_valency; iter++) {

        vector<VertexID> vloop0;
        circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
            vloop0.push_back(v);
        }) );

        vector<VertexID> vloop1;
        circulate_face_ccw(m, f1, std::function<void(VertexID)>([&](VertexID v){
            vloop1.push_back(v);
        }) );

        HalfEdgeID bd_edge;
        circulate_face_ccw(m, f0, std::function<void(HalfEdgeID)>([&](HalfEdgeID h){
            if(m.walker(h).vertex() == central_vertex_0)
                bd_edge = h;
        }) );

        size_t L= vloop0.size();
        VertexID split_vertex;

        for(int i = 0; i < L; i++) {
            if(vloop0[i] == central_vertex_0)
                split_vertex = vloop0[(i+3)%L];
        }

        float dot_sum = 0;
        float len = 0;
        for(int j_off = 0; j_off < L; ++j_off) {
            bool center_match = false;
            for(int i=0;i<L;++i) {
                if(vloop0[i] == central_vertex_0 && vloop1[(L + j_off - i)%L] == central_vertex_1)
                    center_match = true;
            }
            if(center_match) {
                len = -FLT_MAX;
                dot_sum = FLT_MAX;
                Vec3d bridge_edge_i, bridge_edge_j;
                for(int i=0;i<L;i++) {
                    if(abs(length(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L]))) > len)
                      len = abs(length(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L])));
                    bridge_edge_i = normalize(m.pos(vloop0[i]) - m.pos(vloop1[(L+j_off - i)%L]));
                    double curr_dot_sum = abs(dot(bridge_edge_i, cv_edge));// + abs(dot(bridge_edge_i, v_n_nn));
                    if(curr_dot_sum < dot_sum)
                        dot_sum = curr_dot_sum;
                    for(int j=0;j<L;j++) {
                        bridge_edge_j = normalize(m.pos(vloop0[j]) - m.pos(vloop1[(L+j_off - j)%L]));
                        curr_dot_sum = abs(dot(bridge_edge_i, bridge_edge_j));
                        if(curr_dot_sum < dot_sum)
                          dot_sum = curr_dot_sum;
                    }
                }

            }
        }
        dot_sum = dot_sum;
        if(iter < central_valency){
            if(dot_sum > max_dot_sum)
                max_dot_sum = dot_sum;
            if(len < min_len)
                min_len = len;
        }
        else if (iter >= central_valency){
            //if(dot_sum == max_dot_sum)
            //    return f0;
            if(len == min_len)
                return f0;

       }

       FaceID new_face = f0;
       if(L == 4) {
           if(m.walker(bd_edge).face() == f0)
               new_face = m.walker(bd_edge).opp().face();
           else
               new_face = m.walker(bd_edge).face();
       }
       else {
       new_face = m.split_face_by_edge(f0, central_vertex_0, split_vertex);

       if(m.in_use(bd_edge))
           m.merge_faces(m.walker(bd_edge).face(),bd_edge);
      }

        f0 = new_face;

    }
    return f0;
}

vector<pair<VertexID, VertexID>> face_match_careful(HMesh::Manifold& m, FaceID &f0, FaceID &f1,
                                                    const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {

    vector<pair<VertexID, VertexID> > connections;
    if(!m.in_use(f0) || !m.in_use(f1))
        return connections;

    vector<VertexID> loop0;
    circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
        loop0.push_back(v);
    }));

    vector<VertexID> loop1;
    circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
        loop1.push_back(v);
    }));

    size_t L0= loop0.size();
    size_t L1= loop1.size();

    if (L0 != L1)
        return connections;

    size_t L = L0;

    if(face_vertex[f0] == InvalidVertexID || face_vertex[f1] == InvalidVertexID) {

        FaceID new_face = InvalidFaceID;

        if(face_vertex[f0] != InvalidVertexID) {
          if(g.valence(n) > 2) {
            int central_valency = valency(m, face_vertex[f0]);
            int rot_iters = L - 2 + (central_valency - 1)*2;
            double min_torsion_score = FLT_MAX;
            new_face = f0;
            for (int rot = 0; rot < 2*rot_iters; rot++) {
              new_face = rotate_bridge_face_set_once(m, new_face);
              double curr_torsion_score = compute_torsion_score(m, new_face, f1);
              if(curr_torsion_score < min_torsion_score)
                min_torsion_score = curr_torsion_score;
            }
            for (int rot = 0; rot < 2*rot_iters; rot++) {
              new_face = rotate_bridge_face_set_once(m, new_face);
              double curr_torsion_score = compute_torsion_score(m, new_face, f1);
              if(curr_torsion_score == min_torsion_score)
                break;
            }
            if(new_face != InvalidFaceID) {
              f0 = new_face;
              loop0.clear();
              circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
                  loop0.push_back(v);
              }));
            }
          }
        }

        if(face_vertex[f1] != InvalidVertexID) {
          if(g.valence(nn) > 2) {
            int central_valency = valency(m, face_vertex[f1]);
            int rot_iters = L - 2 + (central_valency - 1)*2;
            double min_torsion_score = FLT_MAX;
            new_face = f1;
            compute_torsion_score(m, f0, new_face);
            for (int rot = 0; rot < 2*rot_iters; rot++) {
              new_face = rotate_bridge_face_set_once(m, new_face);
              double curr_torsion_score = compute_torsion_score(m, f0, new_face);
              if(curr_torsion_score < min_torsion_score)
                min_torsion_score = curr_torsion_score;
            }
            for (int rot = 0; rot < 2*rot_iters; rot++) {
              new_face = rotate_bridge_face_set_once(m, new_face);
              double curr_torsion_score = compute_torsion_score(m, f0, new_face);
              if(curr_torsion_score == min_torsion_score)
                break;
            }
            if(new_face != InvalidFaceID) {
              f1 = new_face;
              loop1.clear();
              circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
                  loop1.push_back(v);
              }));
            }
          }
        }
        float min_len = FLT_MAX;
        int j_off_min_len = -1;
        for(int j_off = 0; j_off < L; j_off = j_off + 1) {
            float len = 0;
            for(int i=0;i<L;++i)
                len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
            if(len < min_len)   {
                j_off_min_len = j_off;
                min_len = len;
            }
        }

        for (int i = 0; i < L; i++) {

            VertexID v0 = loop0[i];
            VertexID v1 = loop1[(L + j_off_min_len - i)%L];

            if(face_vertex[f1] == v1) {
                Walker w = m.walker(f0);
                face_vertex[f0] = v0;
                face_vertex[m.walker(w.halfedge()).opp().face()] = v0;
            }

            else if (face_vertex[f0] == v0) {

                Walker w = m.walker(f1);
                face_vertex[f1] = v1;
                face_vertex[m.walker(w.halfedge()).opp().face()] = v1;
            }

        }

    for(int i=0;i<L;++i)
      connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));
    }

    else {
        VertexID central_vertex_0 = face_vertex[f0];
        VertexID central_vertex_1 = face_vertex[f1];
        int central_valency = valency(m, central_vertex_0);
        int rot_iters = L - 2 + (central_valency - 1)*2;
        double min_torsion_score = -FLT_MAX;
        FaceID new_face = f0;

        for (int rot = 0; rot < 2*rot_iters; rot++) {
          new_face = rotate_bridge_face_set_once(m, new_face);
          double curr_torsion_score = compute_torsion_score(m, new_face, f1);
          if(curr_torsion_score > min_torsion_score)
            min_torsion_score = curr_torsion_score;
        }

        for (int rot = 0; rot < 2*rot_iters; rot++) {
          new_face = rotate_bridge_face_set_once(m, new_face);
          double curr_torsion_score = compute_torsion_score(m, new_face, f1);
          if(curr_torsion_score == min_torsion_score)
            break;
        }

        if(new_face != InvalidFaceID)
            f0 = new_face;
        else {
            return connections;
        }

        if(f0 == InvalidFaceID) {
            return connections;
        }


        loop0.clear();
        circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
            loop0.push_back(v);
        }) );

        loop1.clear();
        circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
            loop1.push_back(v);
        }) );

        int j_off_min_len = -1;
        for(int j_off = 0; j_off < L; j_off = j_off + 1) {
            bool center_match = false;
            for(int i=0;i<L;++i) {
                if(loop0[i] == central_vertex_0 && loop1[(L + j_off - i)%L] == central_vertex_1)
                    center_match = true;

            }
            if(center_match) {
                j_off_min_len = j_off;
            }
        }

    for(int i=0;i<L;++i)
        connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));

    }


    return connections;
}


vector<pair<VertexID, VertexID>> face_match_one_ring(HMesh::Manifold& m, FaceID &f0, FaceID &f1,
                                                     const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {

    vector<pair<VertexID, VertexID> > connections;
    if(!m.in_use(f0) || !m.in_use(f1))
        return connections;

    VertexID face_vertex_0 = one_ring_face_vertex[f0];
    VertexID face_vertex_1 = one_ring_face_vertex[f1];

    bool fv_flag = false;

    Vec3d v_n_nn = normalize(g.pos[n] - g.pos[nn]);


    if(face_vertex_0 != InvalidVertexID && face_vertex_1 != InvalidVertexID)
      fv_flag = true;

    int loop0_index = 0, loop1_index = 0;

    vector<VertexID> loop0;

    int count = 0;

    circulate_face_ccw(m, f0, std::function<void(VertexID)>([&](VertexID v){
        loop0.push_back(v);
        if(v == face_vertex_0)
          loop0_index = count;
        count++;
    }) );

    vector<VertexID> loop1;
    count = 0;

    circulate_face_ccw(m, f1, std::function<void(VertexID)>( [&](VertexID v) {
        loop1.push_back(v);
        if(v == face_vertex_1)
          loop1_index = count;
        count++;
    }) );

    size_t L0= loop0.size();
    size_t L1= loop1.size();

    if (L0 != L1)
        return connections;

    size_t L = L0;

    int j_off_min_len = -1;

    if(one_ring_face_vertex[f0] == InvalidVertexID || one_ring_face_vertex[f1] == InvalidVertexID) {

      float min_len = FLT_MAX;

      for(int j_off = 0; j_off < L; j_off = j_off + 1) {
        Vec3d bridge_edge_i, bridge_edge_j;
        float len = 0;

        for(int i=0;i<L;++i) {
          len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
        }
           if(len < min_len)   {
            j_off_min_len = j_off;
            min_len = len;
           }
      }

      int found_flag = 0;

      //if(one_ring_face_vertex[f0] != InvalidVertexID || one_ring_face_vertex[f1] != InvalidVertexID)
      //  return connections;

      for (int i = 0; i < L; i++) {

        VertexID v0 = loop0[i];
        VertexID v1 = loop1[(L + j_off_min_len - i)%L];

        if(face_vertex_1 == v1) {
            Walker w = m.walker(f0);
            one_ring_face_vertex[f0] = v0;
            one_ring_face_vertex[m.walker(w.halfedge()).opp().face()] = v0;
            found_flag = 1;
        }

        else if (face_vertex_0 == v0) {
            Walker w = m.walker(f1);
            one_ring_face_vertex[f1] = v1;
            one_ring_face_vertex[m.walker(w.halfedge()).opp().face()] = v1;
            found_flag = 1;
        }
      }

      if(found_flag == 0 && (one_ring_face_vertex[f0] != InvalidVertexID || one_ring_face_vertex[f1] != InvalidVertexID)) {
        connections.clear();
        return connections;
      }
    }
    else {

      for(int j_off = 0; j_off < L; j_off = j_off + 1) {

        bool center_match = false;

        for(int i=0;i<L;++i) {
            if(loop0[i] == one_ring_face_vertex[f0] && loop1[(L + j_off - i)%L] == one_ring_face_vertex[f1])
                center_match = true;

        }
        if(center_match) {
            j_off_min_len = j_off;
        }

      }
      float min_len = FLT_MAX;
      for(int j_off = j_off_min_len; j_off < 2*L; j_off = j_off + 2) {

        //Vec3d bridge_edge_i, bridge_edge_j;
        //float dot_sum = FLT_MAX;
        float len = 0;

        for(int i=0;i<L;++i) {
          len += sqr_length(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
          /**bridge_edge_i = normalize(m.pos(loop0[i]) - m.pos(loop1[(L+j_off - i)%L]));
          double curr_dot_sum = abs(dot(bridge_edge_i, v_n_nn));// + abs(dot(bridge_edge_i, v_n_nn));
          if(curr_dot_sum < dot_sum)
              dot_sum = curr_dot_sum;
          for(int j=0;j<L;j++) {
              bridge_edge_j = normalize(m.pos(loop0[j]) - m.pos(loop1[(L+j_off - j)%L]));
              curr_dot_sum = abs(dot(bridge_edge_i, bridge_edge_j));
              if(curr_dot_sum < dot_sum)
                dot_sum = curr_dot_sum;
        }**/
      }

      if(len < min_len)   {
       j_off_min_len = j_off;
       min_len = len;
      }


/**      if(dot_sum > max_dot_sum) {
          j_off_min_len = j_off;
          max_dot_sum = dot_sum;
      }**/
     }
   }

   for(int i=0;i<L;++i)
    connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));


   return connections;
}

FaceID find_bridge_face(HMesh::Manifold &m_out,
                        const Geometry::AMGraph3D& g, NodeID start_node, NodeID next_node, Util::AttribVec<NodeID, FaceSet>& node2fs) {

  FaceID f0;

  auto best_face = [&](NodeID n, NodeID nn) {
      Vec3d pn = g.pos[n];
      Vec3d pnn = g.pos[nn];
      Vec3d v_n_nn = pnn - pn;
      double d_max = -1000;
      FaceID f_max = InvalidFaceID;
      v_n_nn = normalize(v_n_nn);
      for(auto f: node2fs[n]) {
          double d = dot(v_n_nn, normal(m_out, f));
          if(d> d_max) {
              f_max = f;
              d_max = d;
          }

      }
      node2fs[n].erase(f_max);
      return f_max;
  };

  if(g.neighbors(start_node).size()>2)  {

    auto key = std::make_pair(start_node,next_node);
    f0 = branchface.find(key)->second;
    VertexID v0 = branch_best_vertex.find(key)->second;
    if(v0 != InvalidVertexID)
        face_vertex[f0] = v0;

  }
  else
    f0 = best_face(start_node,next_node);

  return f0;

}

vector<pair<VertexID, VertexID>> find_bridge_connections(HMesh::Manifold &m_out, FaceID &f0, FaceID &f1,
                                                         const Geometry::AMGraph3D& g, NodeID n, NodeID nn) {
  vector<pair<VertexID, VertexID>> connections;

  if(f0 == InvalidFaceID || f1 == InvalidFaceID)
    return connections;

  if(face_vertex[f0] == InvalidVertexID && face_vertex[f1] == InvalidVertexID)
    connections = face_match_one_ring(m_out, f0, f1, g, n, nn);

  else if (f0 != InvalidFaceID && f1 != InvalidFaceID)
    connections = face_match_careful(m_out, f0, f1, g, n, nn);

  return connections;

}

//Setup Global arrays

void init_graph_arrays(HMesh::Manifold &m_out, const Geometry::AMGraph3D& g, Util::AttribVec<NodeID, FaceSet>& node2fs) {
  init_branch_degree(m_out, g, node2fs);
  init_branch_face_pairs(m_out, g, node2fs);
  merge_branch_faces(m_out, g, node2fs);
}

//Main functions

HMesh::Manifold graph_to_FEQ(const Geometry::AMGraph3D& g, const vector<double>& _node_radii) {

    double r = 0.5 * g.average_edge_length();
    Manifold m_out;
    Util::AttribVec<NodeID, FaceSet> node2fs;

    clear_global_arrays();

    vector node_radii = _node_radii;
    node_radii.resize(g.no_nodes());
    for(auto n : g.node_ids())
      if(node_radii[n] == 0.0)
          node_radii[n] = r;

    construct_bnps(m_out, g, node2fs, node_radii);
    //refine_BNPs(m_out, g, node2fs);
    id_preserving_cc(m_out);
    init_graph_arrays(m_out, g, node2fs);

    FaceAttributeVector<int> ftouched(m_out.allocated_faces(),-1);
    val2nodes_to_boxes(g, m_out, node2fs, node_radii);

    for(auto f_id: m_out.faces()) {
        face_vertex[f_id] = InvalidVertexID;
        if(one_ring_face_vertex.find(f_id) == one_ring_face_vertex.end())
            one_ring_face_vertex[f_id] = InvalidVertexID;
    }

    bool has_junction = false;

    for (auto n: g.node_ids()) {
      if(g.valence(n) > 2)
        has_junction = true;
    }

    for (auto n: g.node_ids()) {
        FaceID f0 = InvalidFaceID;
        VertexID v0, v1;

        auto N = g.neighbors(n);

        if (N.size()<=2 && has_junction)
            continue;

        for(auto nn: N) {
            auto key = std::make_pair(n,nn);
            f0 = branchface.find(key)->second;

            if(branchdeg.find(key)->second < 1 && has_junction)
                continue;

            NodeID start_node = n;
            NodeID next_node = nn;

            vector<NodeID> nbd_list = next_neighbours(g, start_node, next_node);

            do {

                FaceID f0 = find_bridge_face(m_out, g, start_node, next_node, node2fs);
                FaceID f1 = find_bridge_face(m_out, g, next_node, start_node, node2fs);

                nbd_list = next_neighbours(g, start_node, next_node);


                if(g.valence(next_node) > g.valence(start_node)) {
                  auto connections = find_bridge_connections(m_out, f1, f0, g, next_node, start_node);
                  if(connections.size()!=0) {
                    m_out.bridge_faces(f1,f0,connections);
                    ftouched[f0] = 1;
                    ftouched[f1] = 1;
                  }
                }
                else {
                  auto connections = find_bridge_connections(m_out, f0, f1, g, start_node, next_node);
                  if(connections.size()!=0) {
                    m_out.bridge_faces(f0,f1,connections);
                    ftouched[f0] = 1;
                    ftouched[f1] = 1;
                  }
                }

                start_node = next_node;
                if(nbd_list.size()==1)
                  next_node = nbd_list[0];

            }  while(nbd_list.size()==1);
        }
    }

    quad_mesh_leaves(m_out);
    return m_out;
}
