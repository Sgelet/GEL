#ifndef dyncon_hpp
#define dyncon_hpp

#include <iostream>
#include <list>
#include <stack>
#include <unordered_map>
#include <set>
#include <map>
#include <queue>

namespace Geometry {
    enum BBT { Splay, Treap };
    template<typename T, BBT TT>
    class DynCon {
    private:
        struct UndirectedEdge {
            // like a pair but (first, second) is the same as (second, first).
            T first, second;
            // Location in adjacency list.
            typename std::list<T>::iterator adjListFirst, adjListSecond;

            UndirectedEdge(const T a, const T b) {
                first = std::min(a, b);
                second = std::max(a, b);
            }

            bool operator==(const UndirectedEdge &e) const {
                return (first == e.first && second == e.second);
            }
        };

        struct UndirectedEdgeHash {
            size_t operator()(const UndirectedEdge &e) const {
                uint32_t v = e.first, w = e.second;
                return std::hash<uint64_t>()((uint64_t) v << 32 | w);
            }
        };

        struct Node {
            T u, v; // endpoints and size of subtree
            bool adjT, adjNT, marked; // subtree has tree/nontree to consider, node itself should be considered
            size_t size,l, r, p; // left,right,parent
            uint priority;

            Node(T x, T y) {
                u = x;
                v = y;
                size = 1;
                adjNT = false;
                marked = x < y;
                adjT = marked;
                l = r = p = -1;
                priority = 0;
            }
        };

        class Sequence{
        private:
            std::vector<Node> v;

            inline void make_rc(size_t c, size_t p){
                if(c!=-1) v[c].p = p;
                v[p].r = c;
            }
            inline void make_lc(size_t c, size_t p){
                if(c!=-1) v[c].p = p;
                v[p].l = c;
            }
            inline void disc_rc(size_t t){
                if(!has_r(t)) return;
                v[rc(t)].p = -1;
                v[t].r = -1;
            }
            inline void disc_lc(size_t t){
                if(!has_l(t)) return;
                v[lc(t)].p = -1;
                v[t].l = -1;
            }

            void update_data(size_t t) {
                if (t!=-1) {
                    v[t].size = 1 + count(lc(t)) + count(rc(t));
                    v[t].adjNT = hasAdjNT(lc(t)) || hasAdjNT(rc(t)) || (v[t].u == v[t].v && v[t].marked);
                    v[t].adjT = hasAdjT(lc(t)) || hasAdjT(rc(t)) || (v[t].u < v[t].v && v[t].marked);
                }
            }

            void propagate_data(size_t t) {
                bool pastNT;
                bool pastT;
                size_t pastSize;
                while (t!=-1) {
                    pastNT = v[t].adjNT;
                    pastT = v[t].adjT;
                    pastSize = v[t].size;
                    update_data(t);
                    if ((pastNT != v[t].adjNT) || (pastT != v[t].adjT) || (pastSize != v[t].size)) t = p(t);
                    else break;
                }
            }

            void rotateR(size_t t) {
                size_t l = lc(t);
                v[t].l = rc(l);
                if (has_r(l)) v[rc(l)].p = t;
                v[l].p = p(t);
                if (has_p(t)) {
                    if (lc(p(t)) == t) v[p(t)].l = l;
                    else v[p(t)].r = l;
                }
                v[l].r = t;
                v[t].p = l;
                update_data(t);
                update_data(l);
            }

            void rotateL(size_t t) {
                size_t r = rc(t);
                v[t].r = lc(r);
                if (has_l(r)) v[lc(r)].p = t;
                v[r].p = p(t);
                if (has_p(t)) {
                    if (lc(p(t)) == t) v[p(t)].l = r;
                    else v[p(t)].r = r;
                }
                v[r].l = t;
                v[t].p = r;
                update_data(t);
                update_data(r);
            }

            void splay(size_t t) {
                while (has_p(t)) {
                    if (!has_p(p(t))) {
                        if (lc(p(t)) == t) rotateR(p(t));
                        else rotateL(p(t));
                    } else if (lc(p(t)) == t && lc(p(p(t))) == p(t)) {
                        rotateR(p(p(t)));
                        rotateR(p(t));
                    } else if (rc(p(t)) == t && rc(p(p(t))) == p(t)) {
                        rotateL(p(p(t)));
                        rotateL(p(t));
                    } else if (rc(p(t)) == t && lc(p(p(t))) == p(t)) {
                        rotateL(p(t));
                        rotateR(p(t));
                    } else {
                        rotateR(p(t));
                        rotateL(p(t));
                    }
                }
            }


            void s_join(size_t& t, size_t l, size_t r) {
                if(r==-1) {t = l; return;}
                if(l==-1) {t = r; return;}

                size_t max = l;
                while (has_r(max)) max = rc(max);
                splay(max);
                v[max].r = r;
                v[r].p = max;
                update_data(max);
                t = max;
            }

            void s_split(size_t t, size_t &l, size_t &r) {
                splay(t);
                l = lc(t);
                v[t].l = -1;
                r = t;
                if (l!=-1) {
                    v[l].p = -1;
                    update_data(l);
                }
                update_data(r);
            }

            void t_join(size_t& t, size_t l, size_t r) {
                if(r==-1){t = l; return;}
                if(l==-1){t = r; return;}


                //std::cout<<"Joining: "<<std::endl;
                //verify_children(l);
                //std::cout<<std::endl<<"and: "<<std::endl;
                //verify_children(r);
                //std::cout<<std::endl;


                size_t root;
                bool go_r = v[l].priority > v[r].priority;
                bool went_r = go_r;
                if(go_r){
                    root = l;
                    v[root].size += v[r].size;
                    v[root].adjNT |= v[r].adjNT;
                    v[root].adjT |= v[r].adjT;
                    l = rc(l);
                } else {
                    root = r;
                    v[root].size += v[l].size;
                    v[root].adjNT |= v[l].adjNT;
                    v[root].adjT |= v[l].adjT;
                    r = lc(r);
                }
                t = root;
                while(l != -1 && r != -1){
                    go_r = v[l].priority > v[r].priority;
                    if(go_r){
                        if(went_r) make_rc(l,root);
                        else make_lc(l,root);
                        root = l;
                        v[root].size += v[r].size;
                        v[root].adjNT |= v[r].adjNT;
                        v[root].adjT |= v[r].adjT;
                        l = rc(l);
                    } else {
                        if (went_r) make_rc(r, root);
                        else make_lc(r, root);
                        root = r;
                        v[root].size += v[l].size;
                        v[root].adjNT |= v[l].adjNT;
                        v[root].adjT |= v[l].adjT;
                        r = lc(r);
                    }
                    went_r = go_r;
                }
                if(l!=-1){
                    if(went_r) make_rc(l,root);
                    else make_lc(l,root);
                } else {
                    if (went_r) make_rc(r, root);
                    else make_lc(r, root);
                }

                //std::cout<<"Join produced: "<<std::endl;
                //verify_children(t);
                //std::cout<<std::endl;

            }

            void t_split(size_t t, size_t &l, size_t &r) {

                //std::cout<<"Splitting: "<<std::endl;
                //verify_children(find_representative(t));
                //std::cout<<std::endl;

                l = lc(t);
                r = t;
                disc_lc(t);
                update_data(t);
                bool went_r = true;
                bool going_r;
                while(has_p(t)){ // Traverse up
                    going_r = lc(p(t))==t;
                    t = p(t);
                    if(going_r && went_r) r=t; // Continue R
                    else if(!(going_r||went_r)) l=t; // Continue L
                    else if(going_r) { // Switch from L to R
                        disc_lc(t);
                        make_lc(r,t);
                        r = t;
                    } else { // Switch from R to L
                        disc_rc(t);
                        make_rc(l,t);
                        l = t;
                    }
                    update_data(t);
                    went_r = going_r;
                }
                update_data(t);
                //std::cout << "Split has L: "<<std::endl;
                //verify_children(l);
                //std::cout<<std::endl<<"Split has R: "<<std::endl;
                //verify_children(r);
                //std::cout<<std::endl;

            }

        public:
            explicit Sequence(size_t capacity){
                v.reserve(capacity);
            }

            inline size_t p(size_t t){return v[t].p;}
            inline size_t lc(size_t t){return v[t].l;}
            inline size_t rc(size_t t){return v[t].r;}

            inline bool has_p(size_t t){return v[t].p != -1;}
            inline bool has_l(size_t t){return v[t].l != -1;}
            inline bool has_r(size_t t){return v[t].r != -1;}

            size_t add(T p, T q){
                v.emplace_back(Node(p,q));
                if constexpr(TT == Treap) v.back().priority = std::hash<uint32_t>()((uint32_t) p << 16 | q);
                return v.size()-1;
            }

            size_t find_representative(size_t t) {
                if(t==-1) return -1;
                while (has_p(t)) {
                    //std::cout << "("<<t->u<<","<<t->v << ")-";
                    t = p(t);
                }
                return t;
            }

            int count(size_t t) {
                return t!=-1 ? v[t].size : 0;
            }

            bool hasAdjNT(size_t t) {
                return t != -1 && v[t].adjNT;
            }

            bool hasAdjT(size_t t) {
                return t != -1 && v[t].adjT;
            }

            void join(size_t& t, size_t l, size_t r) {
                if constexpr(TT == Splay) s_join(t,l,r);
                else if constexpr(TT == Treap) t_join(t,l,r);
            }

            void split(size_t t, size_t &l, size_t &r) {
                if constexpr(TT == Splay) s_split(t,l,r);
                else if constexpr(TT == Treap) t_split(t,l,r);
            }

            size_t remove_first(size_t t) {
                size_t res = rc(t);
                if(res != -1) disc_rc(t);
                if(has_p(t)){
                    make_lc(res,p(t));
                    propagate_data(p(t));
                    res = p(t);
                }
                v[t].l = v[t].r = v[t].p = -1;
                //verify_children(find_representative(res));
                return find_representative(res);
            }

            void mark(size_t t, bool val){
                v[t].marked = val;
                propagate_data(t);
            }

            size_t reroot(size_t t) {
                size_t A, B;
                split(t, A, B);
                join(B,B, A);
                return B;
            }

            void verify_children(size_t t) {
                if (t==-1) return;
                if (has_p(t) && !(lc(p(t)) == t || rc(p(t)) == t))
                    std::cout << "ERROR P " << "(" << v[t].u << "," << v[t].v << ")" << "," << "(" << v[p(t)].u << "," << v[p(t)].v
                              << ")" << std::endl;
                if (has_l(t) && (p(lc(t)) != t))
                    std::cout << "ERROR L " << "(" << v[t].u << "," << v[t].v << ")" << "," << "(" << v[lc(t)].u << "," << v[lc(t)].v
                              << ")" << std::endl;
                if (has_r(t) && (p(rc(t)) != t))
                    std::cout << "ERROR R " << "(" << v[t].u << "," << v[t].v << ")" << "," << "(" << v[rc(t)].u << "," << v[rc(t)].v
                              << ")" << std::endl;
                if (v[t].size != 1 + count(lc(t)) + count(rc(t))){
                    std::cout << "ERROR SIZE" << std::endl;
                }
                verify_children(lc(t));
                //std::cout<<"("<<v[t].u<<","<<v[t].v<<")";
                verify_children(rc(t));
            }

            Node access(size_t t){
                return v[t];
            }
        };


        class EulerTourForest {
        private:
            std::map<std::pair<T, T>, size_t> v_map;
            Sequence& stw;

        public:
            explicit EulerTourForest(DynCon::Sequence& _stw) : stw(_stw){}

            size_t add(T u, T v){
                size_t t = stw.add(u,v);
                v_map[std::pair<T, T>(u, v)] = t;
                return t;
            }

            size_t find_tree(T u, T v) {
                auto search = v_map.find(std::pair<T, T>(u, v));
                if (search == v_map.end()) return -1;
                else return search->second;
            }

            size_t find_tree(T x) {
                return find_tree(x, x);
            }

            size_t link(T u, T v) {
                size_t t_u = find_tree(u);
                size_t t_v = find_tree(v);

                if (t_u == -1) t_u = add(u,u);
                if (t_v == -1) t_v = add(v,v);

                //stw.splay(t_u);
                //stw.splay(t_v);
                //if (stw.p(t_u) == t_v || (stw.has_p(t_u) && stw.p(stw.p(u)) == t_v)) return t_v;

                if(stw.find_representative(t_u) == stw.find_representative(t_v)) return t_v;

                t_u = stw.reroot(t_u);
                t_v = stw.reroot(t_v);

                auto forward = add(u,v);
                auto retreat = add(v, u);

                stw.join(t_u, t_u, forward);
                stw.join(t_v, t_v, retreat);
                stw.join(t_u, t_u, t_v);

                return t_u;
            }

            // Cuts the edge between u and v. A and B will be the roots of the new trees after the cut (not the first element in ETT sequence).
            void cut(T u, T v, size_t &A, size_t &B) {
                auto s_advance = v_map.find(std::pair<T, T>(u, v));
                auto s_retreat = v_map.find(std::pair<T, T>(v, u));
                if (s_advance == v_map.end() || s_retreat == v_map.end()) {
                    A = B = -1;
                    return;
                }
                size_t J, K, L;
                size_t t_advance = s_advance->second;
                size_t t_retreat = s_retreat->second;

                stw.split(t_advance, J, L);
                //L = stw.rc(L);
                L = stw.remove_first(t_advance);
                v_map.erase(s_advance);

                if (stw.find_representative(t_retreat) == J) {
                    stw.split(t_retreat, J, K);
                    K = stw.remove_first(t_retreat);
                } else {
                    stw.split(t_retreat, K, L);
                    L = stw.remove_first(t_retreat);
                }

                //stw.remove_first(t_retreat);
                v_map.erase(s_retreat);
                /*
                std::cout << "J: "<<std::endl;
                stw.verify_children(J);
                std::cout<<std::endl;
                std::cout << "K: "<<std::endl;
                stw.verify_children(K);
                std::cout<<std::endl;
                std::cout << "L: "<<std::endl;
                stw.verify_children(L);
                std::cout<<std::endl;
                */
                A = K;
                stw.join(B, J, L);
            }

            bool is_connected(T u, T v) {
                size_t t_u = find_tree(u);
                size_t t_v = find_tree(v);

                if (!(t_u!=-1 && t_v!=-1)) return false;

                //stw.splay(t_u);
                //stw.splay(t_v);
                //return (stw.p(t_u) == t_v || (stw.has_p(t_u) && stw.p(stw.p(t_u)) == t_v));
                return stw.find_representative(t_u) == stw.find_representative(t_v);
            }

            bool edge_exists(T u, T v) {
                return find_tree(u, v) != -1;
            }

            void mark(T x, T y, bool val) {
                auto t = find_tree(x, y);
                if (t == -1) {t = add(x,y);}

                stw.mark(t,val);
            }

            void mark(T x, bool val) { mark(x, x, val); }
        };

        EulerTourForest forest;

        // An incident list for each vertex. Uses adjacency lists since we then only have to store one int.
        std::map<T, std::list<T>*> adjacencyLists;

        std::unordered_map<UndirectedEdge, size_t, UndirectedEdgeHash> edgeSet;
        std::vector<UndirectedEdge> ev;

        Sequence stw = Sequence(64);

        size_t get_rep_tree(T v){
            size_t rep = forest.find_tree(v);
            if (rep == -1) return -1;
            return stw.find_representative(rep);
        }

        // Adds a non-tree edge with the given specifications
        void addNonTree(T v, T w, size_t e) {
            if (v > w) std::swap(v, w);

            // Mark that the corresponding nodes in forest has adjacent non-tree edges
            forest.mark(v, true);
            forest.mark(w, true);

            // Add e to adjacency lists.
            if (adjacencyLists.count(v) == 0) {
                adjacencyLists[v] = new std::list<T>();
            }
            if (adjacencyLists.count(w) == 0) {
                adjacencyLists[w] = new std::list<T>();
            }
            auto listV = adjacencyLists.find(v)->second;
            auto listW = adjacencyLists.find(w)->second;

            listV->push_front(w);
            listW->push_front(v);

            // Store iterators to positions in adjacency lists for O(1) removal
            ev[e].adjListFirst = listV->begin();
            ev[e].adjListSecond = listW->begin();
        }

        void disconnect_nontree(size_t e) {
            int v = *ev[e].adjListSecond;
            int w = *ev[e].adjListFirst;
            auto listV = adjacencyLists.find(v)->second;
            auto listW = adjacencyLists.find(w)->second;
            listV->erase(ev[e].adjListFirst);
            listW->erase(ev[e].adjListSecond);
            if (listV->empty()) forest.mark(v, false);
            if (listW->empty()) forest.mark(w, false);
        }

        bool disconnect(T v, T w){
            if (v > w) std::swap(v, w);

            // Early return if edge doesn't exist
            auto search = edgeSet.find(UndirectedEdge(v, w));
            if (search == edgeSet.end()) return false;

            auto e = search->second;

            // Remove from edgeset
            edgeSet.erase(ev[e]);

            // If edge is non-tree, removal cannot break connectivity
            if (!forest.edge_exists(v, w)) {
                // But we must update adjacency lists
                disconnect_nontree(e);
                return false;
            }

            size_t V, W;

            forest.cut(v, w, V, W);

            return true;
        }

        bool reconnect(T v, T w){
            if(is_connected(v,w)) return true;

            size_t replacement = -1;

            size_t V = get_rep_tree(v);
            size_t W = get_rep_tree(w);

            // Ensure that |V| < |W|. We only need to use V.
            if (stw.access(V).size > stw.access(W).size) std::swap(V, W);

            // If there are no adjacent non-tree edges in V no need for searching
            if (!stw.hasAdjNT(V)) {
                return false;
            }

            // Iterative level order traversal of V and non-tree edges
            std::queue<size_t> queue;
            size_t current = V;
            queue.push(current);


            while (replacement==-1 && !queue.empty()) { // Process
                current = queue.front();
                queue.pop();

                if (stw.has_l(current) && stw.hasAdjNT(stw.lc(current))) queue.push(stw.lc(current));
                if (stw.has_r(current) && stw.hasAdjNT(stw.rc(current))) queue.push(stw.rc(current));

                if (stw.access(current).marked &&
                    stw.access(current).u == stw.access(current).v) { // Node represents vertex with adjacent non-tree edges
                    auto list = adjacencyLists.find(stw.access(current).u)->second;
                    auto iter = list->begin();
                    while (iter != list->end()) { // Destructive iteration of adjacent non-tree edges
                        T candidate = *iter;
                        iter++;
                        auto c_e = edgeSet.find(UndirectedEdge(stw.access(current).u, candidate))->second;
                        if (stw.find_representative(forest.find_tree(candidate)) != V) { // Non-tree edge goes to W
                            disconnect_nontree(c_e);
                            replacement = c_e;
                            forest.link(ev[c_e].first, ev[c_e].second);
                            break;
                        }
                    }
                }
            }
            if (replacement!=-1) {
                return true;
            }
            return false;
        }

    public:
        DynCon() : forest(EulerTourForest(stw)) {}

        ~DynCon() {
            for (auto l: adjacencyLists) {
                delete l.second;
            }
        }

        // Insert the edge going from v to w
        int insert(T v, T w) {
            // Swap v and w so that an edge (v, w) is the same as (w, v).
            if (v > w) std::swap(v, w);

            // Edge already exists - early return
            if (edgeSet.find(UndirectedEdge(v, w)) != edgeSet.end()) {
                return 0;
            }

            auto e = ev.size();
            ev.emplace_back(UndirectedEdge(v,w));
            edgeSet[ev[e]] = e;

            // if v & w disconnected in F_0 -> insert as tree edge in F_0.
            if (!is_connected(v, w)) {
                forest.link(v, w);
                return 2;
            } else { // Add as non-tree otherwise
                addNonTree(v, w, e);
                return 1;
            }
        }

        // Returns false if no replacement edge found
        bool remove(T v, T w) {
            if(!disconnect(v,w)) return true;
            return reconnect(v,w);
        }

        // Batch removes every edge adjacent to given vertex
        // Returns false if neighbourhood was not reconnected
        bool remove(T v, const std::vector<T>& adj){
            std::vector<T> adj_tree;
            for(auto w: adj){
                if(disconnect(v,w)) adj_tree.push_back(w);
            }
            uint reconnected = 0;
            for(uint i=0; i<adj_tree.size(); ++i){
                for(uint j=i+1; j<adj_tree.size(); ++j){
                    if(reconnect(adj_tree[i],adj_tree[j])) reconnected++;
                }
            }
            return reconnected == adj_tree.size();
        }

        // Returns true if there exists a path going between v and w.
        bool is_connected(T v, T w) {
            // Check connected in F_0.
            return forest.is_connected(v, w);
        }

        // Returns representative of component of v
        T get_representative(T v) {
            auto rep = get_rep_tree(v);
            if(rep == -1){
                return v;
            }
            return stw.access(rep).u;
        }

        T get_size(T v) {
            // Returns size of component of v
            // Structure actually stores |V|+2*|E| nodes
            // Means (size+2)/3 vertices are in structure
            size_t rep = forest.find_tree(v);
            if (rep == -1) return 1;
            return (stw.access(stw.find_representative(rep)).size + 2) / 3;
        }

        void print_tree(T v) {
            in_order(find_representative(forest.find_tree(v)));
        }
    };
}

#endif //dyncon_hpp
