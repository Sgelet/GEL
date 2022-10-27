//
// Created by usr on 18/07/2022.
//

#include <PyGEL/graph_functions.h>
#include "GEL/Geometry/graph_io.h"
#include "GEL/Geometry/Graph.h"

#include "GEL/HMesh/Manifold.h"
#include "chrono"
#include "GEL/Geometry/graph_skeletonize.h"
#include "GEL/Geometry/DynCon.h"

#ifndef MULTISCALE
#define MULTISCALE 1
#endif

using Graph = Geometry::AMGraph3D;
using SamplingType = Geometry::SamplingType;

void skeletonize(Graph &g, Graph* skel_ptr, SamplingType sampling){

    auto t0 = std::chrono::high_resolution_clock::now();

    Geometry::NodeSetVec seps;
    if constexpr (MULTISCALE) seps = Geometry::multiscale_local_separators(g,sampling,64,0.13);
    else seps = Geometry::local_separators(g,SamplingType::Basic);

    auto t1 = std::chrono::high_resolution_clock::now();

    auto [skel, _]  = skeleton_from_node_set_vec(g, seps);
    *skel_ptr = skel;

    std::cout << "Time to generate: "<<(t1-t0).count() * 1e-9 << std::endl;
    std::cout << "#####################" << std::endl;
}

void graph_from_mesh(HMesh::Manifold* m, Geometry::AMGraph3D* g) {
    HMesh::VertexAttributeVector<Geometry::AMGraph::NodeID> v2n;

    for(auto v : m->vertices())
        v2n[v] = g->add_node(m->pos(v));
    for(auto h: m->halfedges()) {
        HMesh::Walker w = m->walker(h);
        if(h<w.opp().halfedge())
            g->connect_nodes(v2n[w.opp().vertex()], v2n[w.vertex()]);
    }
}

// Count Leaf nodes and cycles
std::pair<uint,uint> count_topology(Geometry::AMGraph3D& g){
    // Nr. of chordless cycles/genus of graph
    g.cleanup();
    uint genus = g.no_edges()-g.no_nodes()+1;
    uint leafs = 0;
    for(auto n: g.node_ids()){
        if(g.edges(n).size()==1) leafs++;
    }
    return {genus,leafs};
}

// Do Haussdorff stuff

int main(int argc, char* argv[]){
    bool success;
    auto g = Graph();

    std::string path = "../package/fertility.off";
    std::string out_path;
    std::string view_mode = "NONE";
    std::string alt_model = "../package/wood_statue.off";   // A different model to show in case the input is a graph. Really only used with wsv.

    if(argc >= 2) path = argv[1];
    if(argc >= 3) out_path = argv[2];

    bool graph_is_mesh = false;

    HMesh::Manifold m;
    // If source is mesh, load and convert to graph
    if(path.compare(path.size()-6,6,".graph")){
        m = HMesh::Manifold();
        success = HMesh::load(path, m);
        if (!success) {
            std::cout << "ERROR : Could not load mesh!";
            return 0;
        }
        graph_from_mesh(&m,&g);
        graph_is_mesh = true;
    } else {
        g = Geometry::graph_load(path);
    }

    success = !g.empty();
    if (!success) std::cout << "ERROR : Graph is empty!" << std::endl;

    std::cout << "Generating skeleton: " <<path<<" saved at "<<out_path<< std::endl;
    std::cout << "Vertices: "<<g.no_nodes()<<"\nEdges: "<<g.no_edges()<<std::endl;
    std::cout << "#####################" << std::endl;

    // Make skeleton
    auto skel = Geometry::AMGraph3D(); // Target graph
    skeletonize(g,&skel,SamplingType::Advanced);

    if(!out_path.empty()) Geometry::graph_save(out_path,skel);

    return 0;
}