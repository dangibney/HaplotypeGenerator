#include <iostream>
#include <stdlib.h>
#include <zlib.h>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <omp.h>
#include <climits>
#include <random>
#include "gfa-priv.h"
#include "gfa.h"
#include "kseq.h"
#include "igraph.h"

// CPP
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>


#define PROB_RECOMBINATION 5*1.0/5000000.0

using namespace std;


typedef struct edge_attributes_struct{
    igraph_integer_t edge_idx;
    igraph_real_t label;
}edge_attributes;


void print_variation_graph(const igraph_t* graph){

    cout << "Printing Variation Graph" << endl;
    igraph_vs_t vs;
    igraph_vit_t vit;

    igraph_vs_all(&vs);
    igraph_vit_create(graph, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        cout << IGRAPH_VIT_GET(vit) << ": " << endl;

        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges, 0);
        igraph_incident(graph, &edges, IGRAPH_VIT_GET(vit), IGRAPH_OUT);
        for(long int i = 0; i < igraph_vector_int_size(&edges); i++){
            igraph_integer_t start, end;
            igraph_edge(graph, igraph_vector_int_get(&edges, i), &start, &end);
            igraph_integer_t eid;
            igraph_get_eid(graph, &eid, start, end, IGRAPH_DIRECTED, true);
            cout << "\t(" << start << ", " << end << ", " 
                << "label: " << igraph_cattribute_EAN(graph, "label", igraph_vector_int_get(&edges, i))
                << ", eid: " << eid << ")" << endl;   

        }
        igraph_vector_int_destroy(&edges);
        IGRAPH_VIT_NEXT(vit);
        cout << endl;
    }
    cout<< endl;
}


void create_variation_graph(std::vector<int>* start,
                            std::vector<int>* end,
                            std::vector<std::string> *labels, 
                            igraph_t* graph){
   
    int num_vertices = -1;

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int start_vertex, end_vertex;
    char edge_label;

    vector<edge_attributes> attributes_vector;

    //cout << start->size() << endl;
    for(int i = 0; i < start->size(); i++)
    {
        start_vertex = start->at(i);
        end_vertex = end->at(i);

        if(start_vertex == end_vertex){
            cout << "SELF-LOOP DETECTED" << endl;
            return;
        }
        edge_label = labels->at(i)[0];
        //cout << i << ": " << start_vertex << " " << end_vertex << " " << edge_label << endl;
        if(start_vertex+1 > num_vertices){
            num_vertices = start_vertex + 1;
        }
        if(end_vertex+1 > num_vertices){
            num_vertices= end_vertex + 1;
        }

        igraph_real_t label;
        if(edge_label == 'A'){
            label = 0;
        }else if(edge_label == 'T'){
            label = 1;
        }else if(edge_label == 'C'){
            label = 2;
        }else if(edge_label == 'G'){
            label = 3;
        }else{
            label = 4;
        }

        //cout << "label:" << label << endl;
        attributes_vector.push_back({i, label});

        //add edge
        igraph_vector_int_push_back(&edges, start_vertex);
        igraph_vector_int_push_back(&edges, end_vertex);
    }

    //cout << "number of vertices: " << num_vertices << endl;
    bool directed = true;
    igraph_empty(graph, num_vertices, directed);
    igraph_add_edges(graph, &edges, NULL);
    igraph_vector_int_destroy(&edges);

    //cout << igraph_ecount(graph) << endl;
    for(edge_attributes a: attributes_vector){
        //cout << a.edge_idx << " " << a.label << endl;
        igraph_cattribute_EAN_set(graph, "label", a.edge_idx, a.label);
    }

}


// coverts adjacency list to edge list for edge labeled graph
// creates new haplotypes walks based on this
void convert_graph(std::vector<std::vector<int>> *adj_list,
                    std::vector<std::string> *node_seq, 
                    std::vector<std::vector<int>> *paths,
                    std::vector<int> *start,
                    std::vector<int> *end,
                    std::vector<std::string> *label,
                    std::vector<std::vector<int>> *haplotypes,
                    std::vector<pair<int,int>> *inverse_map
                    ){

    
    std::vector<std::vector<int>> vertex_expansion;

    int new_vertex_id = adj_list->size();

    for(int i = 0; i < adj_list->size(); i++){

        std::vector<int> expansion;

        //cout << "vertex: " << i << ", seq: " << node_seq->at(i) << endl;
        // vertex sequence length is 1, multiple outputs
        if(node_seq->at(i).length() == 1 && adj_list->at(i).size() > 0){
            for(int j = 0; j < adj_list->at(i).size(); j++){
                start->push_back(i);
                end->push_back(adj_list->at(i)[j]);
                label->push_back(node_seq->at(i));
            }
            expansion.push_back(i);
        // vertex sequence length is 1, no outputs
        }else if(node_seq->at(i).length() == 1 && adj_list->at(i).size() == 0){
            // add new sink vertex
            start->push_back(i);
            end->push_back(new_vertex_id);
            label->push_back(node_seq->at(i));
            expansion.push_back(i);
            expansion.push_back(new_vertex_id);
            new_vertex_id++;
        
        }else{
            //cout << "Longer label detected" << endl;
            // vertex label is longer than 1
            start->push_back(i);
            end->push_back(new_vertex_id);
            label->push_back(node_seq->at(i).substr(0,1));

            expansion.push_back(i);
            expansion.push_back(new_vertex_id);
            for(int k = 1; k < node_seq->at(i).length()-1; k++){
                start->push_back(new_vertex_id);
                
                new_vertex_id++;
                end->push_back(new_vertex_id);
                label->push_back(node_seq->at(i).substr(k,1));
                expansion.push_back(new_vertex_id);
                
            }

            //outgoing edges
            if(adj_list->at(i).size() > 0){
                for(int j = 0; j < adj_list->at(i).size(); j++){
                    start->push_back(new_vertex_id);
                    end->push_back(adj_list->at(i)[j]);
                    label->push_back(node_seq->at(i).substr(node_seq->at(i).length()-1, 1));
                }
                new_vertex_id++;
            }else{
                start->push_back(new_vertex_id);
                new_vertex_id++;
                end->push_back(new_vertex_id);
                expansion.push_back(new_vertex_id);
                label->push_back(node_seq->at(i).substr(node_seq->at(i).length()-1, 1));
            }

        }
        vertex_expansion.push_back(expansion);
    }

    int num_vertices = 0;
    for(int i = 0; i < vertex_expansion.size(); i++){
        num_vertices += vertex_expansion[i].size();
    }

    inverse_map->resize(num_vertices, {0,0});
    for(int i = 0; i < vertex_expansion.size(); i++){
        for(int j = 0; j < vertex_expansion[i].size(); j++){
            inverse_map->at(vertex_expansion[i][j]) =  {i, j};
        }
    }


    // update haplotypes with expansions
    for(int i = 0; i < paths->size(); i++){
        std::vector<int> haplotype;
        for(int j = 0; j < paths->at(i).size(); j++){

            int vertex = paths->at(i)[j];
            for(int k = 0; k < vertex_expansion[vertex].size(); k++){
                haplotype.push_back(vertex_expansion[vertex][k]); 
            }

        }
        haplotypes->push_back(haplotype);
    }
    
}


string generated_hap(igraph_t* graph, 
                    std::vector<std::vector<int>> paths, 
                    vector<string> labels,
                    int target_len
                    ){
    
    cout << "generating haplotypes per vertex list" << endl;
    vector<vector<int>> haps_per_vertex;
    for(int i = 0; i < igraph_vcount(graph); i++){
        haps_per_vertex.push_back({});
    }

    int number_of_haplotypes = paths.size();
    for(int h = 0; h < number_of_haplotypes; h++){
        vector<int> hap = paths.at(h);
        for(int j = 0; j < hap.size(); j++){
            haps_per_vertex[hap.at(j)].push_back(h);
        }
    }


    //cout << "Number of haplotypes: " << number_of_haplotypes  << endl;
    // randomly choose a haplotype start
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist_int(0,number_of_haplotypes-1); // distribution in range [1, 6]

    int hap = dist_int(rng);

    cout << "Starting on haplotype: " << hap + 1<< endl;
    int vertex = paths.at(hap).at(0);
    auto it_current = find(paths.at(hap).begin(), paths.at(hap).end(), vertex);
    //cout << "Vertex " << vertex << endl;

    std::random_device dev2;
    std::mt19937 rng2(dev2());
    std::uniform_real_distribution<> dist_real(0.0,1.0);

    igraph_integer_t out_degree;
    igraph_degree_1(graph, &out_degree, vertex, IGRAPH_OUT, true);

    int number_of_recombinations = 0;
    string out;
    int len = 0;
    int number_of_times_recombination_possible = 0;
    while(out_degree > 0 && len < target_len){

        len++;

        //check if recombination is possible, i.e,
        igraph_vector_int_t neis;
        igraph_vector_int_init(&neis, 0);
        vector<pair<int,int>> possible_recombinations;
        int preferred_neighbor = -1;
        igraph_neighbors(graph, &neis, vertex, IGRAPH_OUT);

        for(int i = 0; i < out_degree; i++){
            int neighbor = VECTOR(neis)[i];
            for(int h : haps_per_vertex.at(neighbor)){
                if(h != hap){
                    possible_recombinations.push_back(make_pair(neighbor, h));
                }else{
                    auto next = it_current;
                    if (next != paths.at(hap).end() && std::next(next) != paths.at(hap).end()) {
                        ++next; // Increment the iterator first
                        if (*next == neighbor) {
                            preferred_neighbor = neighbor;
                        } else {
                            possible_recombinations.push_back(make_pair(neighbor, h));
                        }
                    } else {
                        possible_recombinations.push_back(make_pair(neighbor, h));
                    }
                }
            }  
        }

        if(possible_recombinations.size() > 0){

            number_of_times_recombination_possible++;
            float random_float = dist_real(rng2);
            if(random_float <= PROB_RECOMBINATION){
                uniform_int_distribution<std::mt19937::result_type> dist_int(0,possible_recombinations.size()-1);
                int selection = dist_int(rng);
                int next_vertex = possible_recombinations.at(selection).first;
                hap = possible_recombinations.at(selection).second;
                
                igraph_integer_t eid;
                igraph_get_eid(graph, &eid, vertex, next_vertex, IGRAPH_DIRECTED, true);
                out += labels.at(eid);

                vertex = next_vertex;

                it_current = find(paths.at(hap).begin(), paths.at(hap).end(), vertex);

                cout << "change to: " << hap << endl;
                
                number_of_recombinations++;

            }else if (random_float > PROB_RECOMBINATION && preferred_neighbor != -1){
                int next_vertex = preferred_neighbor;

                igraph_integer_t eid;
                igraph_get_eid(graph, &eid, vertex, next_vertex, IGRAPH_DIRECTED, true);
                out += labels.at(eid);

                vertex = next_vertex;
                ++it_current;

            }else{
                cout << "Ending path early, no recombination and no way to continue on current path" << endl;
                break;
            }
        }
        else if (possible_recombinations.size() <= 0 && preferred_neighbor != -1){
            int next_vertex = preferred_neighbor;

            igraph_integer_t eid;
            igraph_get_eid(graph, &eid, vertex, next_vertex, IGRAPH_DIRECTED, true);
            out += labels.at(eid);

            vertex = next_vertex;
            ++it_current;

        }else{
            cout << "Ending path early, no possible recombination, and no preferred neighbor" << endl;
            break;
        }

        igraph_degree_1(graph, &out_degree, vertex, IGRAPH_OUT, true);
    }

    cout << "Finished" << endl;

    cout << "Length of haplotype generated: " << out.size() << endl;
    cout << "Number of times recombination possible: " << number_of_times_recombination_possible << endl;
    cout << "Number of recombinations: " << number_of_recombinations << endl;

    return out;
}



int main(int argc, char *argv[]) {

    // if argv is less than 3, print usage
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <gfa_file>  <hap out_file>  <target_length (optional)>\n", argv[0]);
        return 1;
    }

    // read the gfa file, reads file and hap file
    cout << "\nReading gfa file..." << endl;
    std::string gfa_file = argv[1];
    std::string hap_file = argv[2];
    

    gfa_t *g = gfa_read(gfa_file.c_str());
    if (g == 0) {
        fprintf(stderr, "[E::%s] failed to load the GFA file\n", __func__);
        return 1;
    }

    // Read the graph
    int32_t n_vtx;
    std::vector<int> node_len;
    std::vector<std::string> node_seq;
    std::vector<std::string> node_name;
    std::vector<std::vector<int>> adj_list;
    std::vector<std::vector<int>> haps;
    std::vector<std::vector<int>> paths;
    std::vector<std::vector<int>> rev_paths;
    int num_walks;
    int lin_ref = 0;
    bool is_gfa_v12 = false;
    uint v;
    n_vtx = gfa_n_vtx(g);
    // Resize node_len
    node_len.resize(n_vtx/2, 0);
    node_seq.resize(n_vtx/2, "");
    node_name.resize(n_vtx/2, "");
    std::vector<std::vector<uint32_t>> adj_;
    adj_.resize(n_vtx); 
    /* Node_len */
    for (int v = 0; v < n_vtx/2; v++)
    {
        gfa_seg_t segment = (g)->seg[v];
        int len =  segment.len;
        node_len[v] = len;
        node_seq[v] = segment.seq;
        node_name[v] = segment.name;
    }
    // look for all the edges , if the sum of all the 
    // edges are zero then, that's a linear reference
    u_int num_edges = 0;
    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        int v_ = av->v_lv >> 32;
        // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::endl; 
        for (int i = 0; i < n_edges; i++)
        {
            num_edges++;
        }
    }


    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        if (num_edges == 0) // Linear Reference
        {
            lin_ref = 1; // Mark as a linear reference
        }else
        {
            int v_ = av->v_lv >> 32;
            // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::endl; 
            for (int i = 0; i < n_edges; i++)
            {
                uint w = av[i].w;
                adj_[v_].push_back(w);
            }
        }
    }

    // Copy the adjacency list only for forward strand
    n_vtx /= 2;
    adj_list.resize(n_vtx);
    for (int i = 0; i < 2 * n_vtx; i++)
    {
        if (i % 2 == 0) // Only forward strand
        {
            for (auto v : adj_[i])
            {
                adj_list[i/2].push_back(v/2);
            }
        }
    }

    adj_.clear();

    // Read walks
    if (g->n_walk > 0) is_gfa_v12 = true;
    num_walks = g->n_walk; // count number of walks
    haps.resize(n_vtx);
    paths.resize(num_walks);
    rev_paths.resize(num_walks);

    // Fill the paths
    for (size_t w = 0; w < g->n_walk; w++)
    //for (size_t w = 1; w < g->n_walk; w++)
    {
       
        rev_paths[w].resize(n_vtx, -1);
        int32_t idx = 0;
        for (size_t n = 0; n < g->walk[w].n_v; n++)
        {
            int v = g->walk[w].v[n];
            if (v%2 != 0) {
                fprintf(stderr, "Error: Walk %d has reverse strand vertices %d\n", w, v);
                exit(1);
            }
            //cout << w << ", " << v << endl;
			v /= 2; // for forward strand
            haps[v].push_back(w); // Map forward walk to haplotype
            paths[w].push_back(v); // Map forward walk to path
            rev_paths[w][v] = idx++;
        }
    }

    /*
    Graph -> adj_list (adj_list[v] =  set of vertices connected to v)
    node sequence -> node_seq
    paths -> paths [paths[idx] = vtx]
    haps -> haps [haps[vtx] = set of walks]
    num_walks -> count of haplotypes
    */


    cout << "Number of haplotypes: " << paths.size() << endl;
    
    std::vector<int> start;
    std::vector<int> end;
    std::vector<std::string> label;
    std::vector<std::vector<int>> haplotypes;
    std::vector<pair<int,int>> inverse_map;

    
    cout << "Converting to edge labeled graph..." << endl;
    convert_graph(&adj_list, &node_seq, &paths, &start, &end, &label, &haplotypes, &inverse_map);

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_t variation_graph;

    create_variation_graph(&start, &end, &label, &variation_graph);

    cout << "Variation graph constructed." << endl;

    int target_len;
    if(argc< 4){
        target_len = igraph_vcount(&variation_graph);
    }else{
        target_len = atoi(argv[3]);
    }

    string out = generated_hap(&variation_graph, haplotypes, label, target_len);

    ofstream outfile(hap_file);
    outfile << ">generated_hap" << endl;
    outfile << out;
    outfile.close();

    return 0;
}
