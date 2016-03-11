#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <set>
#include <stdlib.h>
#include <assert.h>

using namespace std;

class vertex{
  public:
    int id;
    int pid;

    int oedge_id;
    int ori_id;
};

class edge{
  public:
    int begin;
    int end;
    bool is_ori;
    bool is_cut;
    int partition1; // The first partition it belongs
    int partition2; // The second partition it belongs

    int part_id;
};

class graph{
  public:
    int vertex_num;
    int edge_num;
    vertex * vertexes;
    edge * edges;
};

class partition_score{
  public:
    int id;
    int score;
};

class partition{
  public:
    int id;
    int vertex_num;
    int edge_num;
    vector<vertex> vlist;

    // transformed partition
    int ori_vertex_num;
    vector<int> ovlist;
    int ori_edge_num;
    vector<int> oelist;
    vector<partition_score> waitinglist;
};

class original_vertex{
  public:
    int begin;
    int end;
    vector<partition_score> plist;

    int part_id;
};

void updatePartScore(vector<partition_score> * list, graph * g, vertex * v, bool is_first)
{
    vector<partition_score>::iterator it;

    for(it = list->begin(); it != list->end(); it++)
    {
        if(it->id == v->pid)
        {
            assert(v->oedge_id >= 0);
            if(g->edges[v->oedge_id].is_cut)
            {
                it->score += 1;
                if(is_first)
                {
                    if(v == &g->vertexes[g->edges[v->oedge_id].begin-1])
                        updatePartScore(list, g, &g->vertexes[g->edges[v->oedge_id].end-1], false);
                    else if(v == &g->vertexes[g->edges[v->oedge_id].end-1])
                        updatePartScore(list, g, &g->vertexes[g->edges[v->oedge_id].begin-1], false);
                    else
                        assert(false);
                }
            }
            else
                it->score += 2;
            //cout << i << " " << j << endl;
            return;
        }
    }

    // create a new one if not found
    partition_score * ps = new partition_score;
    ps->id = v->pid;
    assert(v->oedge_id >= 0);
    if(g->edges[v->oedge_id].is_cut)
    {
        ps->score = 1;
        if(is_first)
        {
            if(v == &g->vertexes[g->edges[v->oedge_id].begin-1])
                updatePartScore(list, g, &g->vertexes[g->edges[v->oedge_id].end-1], false);
            else if(v == &g->vertexes[g->edges[v->oedge_id].end-1])
                updatePartScore(list, g, &g->vertexes[g->edges[v->oedge_id].begin-1], false);
            else
                assert(false);
        }
    }
    else
        ps->score = 2;
    list->insert(list->end(), *ps);
}

void sortPartScore(vector<partition_score> * list)
{
    int max_score;
    vector<partition_score>::iterator it0, it1, max_idx;

    for(it0 = list->begin(); it0 != list->end(); it0++)
    {
        max_score = -100;
        for(it1 = it0; it1 != list->end(); it1++)
            if(it1->score > max_score)
            {
                max_score = it1->score;
                max_idx = it1;
            }
        if(max_idx == it0)
            continue;

        // swap the elements
        swap(*it0, *max_idx);
    }
}
