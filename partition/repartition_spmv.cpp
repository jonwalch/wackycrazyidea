#include "repartition.h"

//#define SORT_BY_Y

void sortEdgebyBegin(vector<int> * list, graph * g)
{
    vector<int>::iterator it, it1, min_idx;

    for(it = list->begin(); it != list->end(); it++)
    {
        int min_begin = g->vertex_num;
        int min_end = g->vertex_num;
        for(it1 = it; it1 != list->end(); it1++)
        {
            if(g->edges[*it1].begin < min_begin)
            {
                min_begin = g->edges[*it1].begin;
                min_end = g->edges[*it1].end;
                min_idx = it1;
            }
            else if(g->edges[*it1].begin == min_begin && g->edges[*it1].end < min_end)
            {
                min_end = g->edges[*it1].end;
                min_idx = it1;
            }
        }
        if(it == min_idx)
            continue;

        swap(*it, *min_idx);
    }
}

void sortEdgebyEnd(vector<int> * list, graph * g)
{
    vector<int>::iterator it, it1, min_idx;

    for(it = list->begin(); it != list->end(); it++)
    {
        int min_end = g->vertex_num;
        int min_begin = g->vertex_num;
        for(it1 = it; it1 != list->end(); it1++)
        {
            if(g->edges[*it1].end < min_end)
            {
                min_end = g->edges[*it1].end;
                min_begin = g->edges[*it1].begin;
                min_idx = it1;
            }
            else if(g->edges[*it1].end == min_end && g->edges[*it1].begin < min_begin)
            {
                min_begin = g->edges[*it1].begin;
                min_idx = it1;
            }
        }
        if(it == min_idx)
            continue;

        swap(*it, *min_idx);
    }
}

int main(int argc, char *argv[])
{
    /**************************************************/
    // initiation
    /**************************************************/
    if(argc != 5)
    {
        cerr << "There should be four arguments." << endl;
        return 0;
    }

    char filename1[200];
    char filename2[200];
    char filename3[200];
    strcpy(filename1, argv[1]);
    strcpy(filename2, argv[2]);
    strcpy(filename3, argv[3]);
    ifstream gfile(filename1);
    if(gfile == NULL)
    {
        cerr << "Cannot open graph file." << endl;
        return 0;
    }
    ifstream pfile(filename2);
    if(pfile == NULL)
    {
        cerr << "Cannot open partition file." << endl;
        return 0;
    }
    ifstream vmapfile(filename3);
    if(vmapfile == NULL)
    {
        cerr << "Cannot open vertex map file." << endl;
        return 0;
    }
    int edge_per_part = atoi(argv[4]);

    graph g;
    gfile >> g.vertex_num >> g.edge_num;

    // read vertex map
    int ori_vertex_num;
    vmapfile >> ori_vertex_num;
    int * vmap = new int[ori_vertex_num + 1];
    for(int i = 0; i < ori_vertex_num+1; i++)
        vmapfile >> vmap[i];
    int tmp;
    vmapfile >> tmp;
    assert(!vmapfile.good());
    vmapfile.close();

    // read vertexes
    //cout << g.vertex_num << endl;
    g.vertexes = new vertex[g.vertex_num];
    int maxpart = 0;
    for(int i = 0; i < g.vertex_num; i++)
    {
        pfile >> g.vertexes[i].pid;
        g.vertexes[i].id = i + 1;
        g.vertexes[i].oedge_id = -1;
        if(g.vertexes[i].pid > maxpart)
            maxpart = g.vertexes[i].pid;
    }
    //maxpart++; // lld
    pfile.close();
    // create reverse vertex map
    for(int i = 0; i < ori_vertex_num; i++)
    {
        for(int j = vmap[i]; j < vmap[i+1]; j++)
        {
            assert(j-1 >= 0 && j-1 < g.vertex_num);
            g.vertexes[j-1].ori_id = i;
        }
    }

    partition * parts = new partition[maxpart+1];
    for(int i = 0; i < maxpart+1; i++)
    {
        parts[i].id = i;
        parts[i].vertex_num = 0;
        parts[i].edge_num = 0;
        parts[i].vlist.clear();

        parts[i].ori_vertex_num = 0;
        parts[i].ovlist.clear();
        parts[i].ori_edge_num = 0;
        parts[i].oelist.clear();
        parts[i].waitinglist.clear();
    }
    for(int i = 0; i < g.vertex_num; i++)
    {
        parts[g.vertexes[i].pid].vlist.insert(parts[g.vertexes[i].pid].vlist.end(), g.vertexes[i]);
        parts[g.vertexes[i].pid].vertex_num++;
    }

    // calculate data in each partition
    int total_data_loads = 0;
    set<int> ori_total_v_set;
    set<int> * ori_v_set = new set<int>[maxpart+1];
    for(int i = 0; i < maxpart+1; i++)
    {
        ori_v_set[i].clear();
        for(vector<vertex>::iterator it = parts[i].vlist.begin(); it != parts[i].vlist.end(); it++) {
            ori_v_set[i].insert(it->ori_id);
            ori_total_v_set.insert(it->ori_id);
        }
        //cout << "part " << i << " has " << ori_v_set[i].size() << "." << endl;
        total_data_loads += ori_v_set[i].size();
    }
    
    // read edges
    string tmpstr;
    getline(gfile, tmpstr);

    g.edges = new edge[g.edge_num];
    int startpoint = 1;
    int edge_idx = 0;
    int endpoint;
    int num_oricut = 0;
    int num_addcut = 0;
    for(int i = 0; i < g.vertex_num; i++)
    {
        bool first_edge = true;

        getline(gfile, tmpstr);
        stringstream tmpstream;
        tmpstream.str("");
        tmpstream << tmpstr;
        while(tmpstream >> endpoint)
        {
            if(startpoint > endpoint)
            {
                first_edge = false;
                continue;
            }

            g.edges[edge_idx].begin = startpoint;
            g.edges[edge_idx].end = endpoint;
            if(first_edge)
            {
                first_edge = false;
                g.edges[edge_idx].is_ori = true;
                g.vertexes[startpoint-1].oedge_id = edge_idx;
                g.vertexes[endpoint-1].oedge_id = edge_idx;
            }
            else
                g.edges[edge_idx].is_ori = false;
            if(g.vertexes[startpoint-1].pid == g.vertexes[endpoint-1].pid)
            {
                g.edges[edge_idx].is_cut = false;
                parts[g.vertexes[startpoint-1].pid].edge_num++;
            }
            else
            {
                g.edges[edge_idx].is_cut = true;
                if(g.edges[edge_idx].is_ori)
                {
                    num_oricut++;
                    if((ori_v_set[g.vertexes[startpoint-1].pid].find(g.vertexes[endpoint-1].ori_id) == ori_v_set[g.vertexes[startpoint-1].pid].end()) && (ori_v_set[g.vertexes[endpoint-1].pid].find(g.vertexes[startpoint-1].ori_id) == ori_v_set[g.vertexes[endpoint-1].pid].end()))
                        total_data_loads++;
                }
                else {
                    num_addcut++;
                }
            }
            g.edges[edge_idx].partition1 = g.vertexes[startpoint-1].pid;
            g.edges[edge_idx].partition2 = g.vertexes[endpoint-1].pid;
            g.edges[edge_idx].part_id = -1;
            //if(i < 6)
            //    cout << startpoint << " " << endpoint << " " << g.edges[edge_idx].is_ori << endl;
            edge_idx++;
        }
        startpoint++;
    }
    assert(edge_idx == g.edge_num);
    gfile >> endpoint;
    assert(!gfile.good());
    gfile.close();
    cout << "INFO: additional data load is " << total_data_loads - ori_total_v_set.size() << " (" << total_data_loads << " - " << ori_total_v_set.size() << ")" << endl;
    cout << "INFO: original cut is " << num_oricut << ", addtional cut is " << num_addcut << "." << endl;
    //return 0;

    // compute balance
    int max_vertex_num = 0;
    for(int i = 0; i < maxpart+1; i++)
    {
        if(parts[i].vertex_num > max_vertex_num)
            max_vertex_num = parts[i].vertex_num;
        //cout << "Part " << parts[i].id << ": vertex " << parts[i].vertex_num << ", edge " << parts[i].edge_num << endl;
    }
    cout << "Balance factor: " << (double)(max_vertex_num * (maxpart+1)) / (double)g.vertex_num << endl;

    /**************************************************/


    /**************************************************/
    // fill edges into partitions
    /**************************************************/
    vector<int> otherlist;
    otherlist.clear();
    int nocutedge_num = 0;
    int cutedge_num = 0;
    int * part_left = new int[maxpart+1];
    for(int i = 0; i < maxpart+1; i++)
        part_left[i] = 0;
    for(int i = 0; i < g.vertex_num; i++)
    {
        int ori_eid = g.vertexes[i].oedge_id;
        assert(ori_eid >= 0);
        assert(g.edges[ori_eid].is_ori);
        if(g.edges[ori_eid].begin != i+1)
            continue;

        // try to add no cutting edges into partition
        if(!g.edges[ori_eid].is_cut)
        {
            nocutedge_num++;
            int epid = g.edges[ori_eid].partition1;
            if(parts[epid].ori_edge_num < edge_per_part)
            {
                parts[epid].oelist.push_back(ori_eid);
                parts[epid].ori_edge_num++;
            }
            else
            {
                part_left[epid]++;
                otherlist.push_back(ori_eid);
            }
        }
        // add cutting edges into waiting list
        else
        {
            cutedge_num++;
            int epid1 = g.edges[ori_eid].partition1;
            int epid2 = g.edges[ori_eid].partition2;
            assert(epid1 != epid2);
            partition_score ps;
            ps.id = ori_eid;
            ps.score = 1; //TODO
            parts[epid1].waitinglist.push_back(ps);
            parts[epid2].waitinglist.push_back(ps);
        }
    }
    //for(int i = 0; i < maxpart+1; i++)
    //{
    //    if(parts[i].ori_edge_num >= edge_per_part)
    //        cerr << "INFO: partition " << i << " is filled fully in the first process. "
    //             << part_left[i] << " edges are left." << endl;
    //    else
    //        assert(part_left[i] == 0);
    //}
    cerr << "INFO: no cut edge number " << nocutedge_num << "; cut edge number " << cutedge_num << endl;

    // fill cut edges
    int filled_cutedge = 0;
    for(int i = 0; i < maxpart+1; i++)
    {
        if(parts[i].ori_edge_num >= edge_per_part)
            continue;
        sortPartScore(&parts[i].waitinglist);
        //if(parts[i].waitinglist.size() < edge_per_part - parts[i].ori_edge_num)
        //    cerr << "WARN: partition " << i << " doesn't have enough candidates in the waiting list."<< endl;

        vector<partition_score>::iterator it;
        for(it = parts[i].waitinglist.begin(); it != parts[i].waitinglist.end(); it++)
        {
            assert(it->id < g.edge_num);
            if(g.edges[it->id].part_id >= 0)
                continue;
            g.edges[it->id].part_id = i;
            parts[i].oelist.push_back(it->id);
            parts[i].ori_edge_num++;
            filled_cutedge++;
            if(parts[i].ori_edge_num >= edge_per_part)
                break;
        }
    }
    int refill_num = 0;
    // fill the rest edges into otherlist
    for(int i = 0; i < maxpart+1; i++)
    {
        vector<partition_score>::iterator it;
        for(it = parts[i].waitinglist.begin(); it != parts[i].waitinglist.end(); it++)
        {
            if(g.edges[it->id].part_id < 0)
            {
                otherlist.push_back(it->id);
                refill_num++;
            }
        }
    }
    cerr << "INFO: " << filled_cutedge << " cut edges are filled. " << refill_num / 2 << " edges are refilled. (" << otherlist.size() << ")" << endl;

    // fill other edges
    //for(vector<int>::iterator it = otherlist.begin(); it != otherlist.end(); it++)
    //    cout << g.edges[*it].begin << "-" << g.edges[*it].end << " " << g.edges[*it].part_id << "; ";
    //cout << endl;
    for(vector<int>::iterator it = otherlist.begin(); it != otherlist.end(); it++)
    {
        if(g.edges[*it].part_id >= 0)
            continue;
        assert(g.edges[*it].part_id < 0);
        int i;
        for(i = 0; i < maxpart+1; i++)
        {
            if(parts[i].ori_edge_num >= edge_per_part)
                continue;
            g.edges[*it].part_id = i;
            parts[i].oelist.push_back(*it);
            parts[i].ori_edge_num++;
            break;
        }
        if(i == maxpart+1)
        {
            cerr << "ERROR: " << g.edges[*it].begin << "-" << g.edges[*it].end 
                 << " edges cannot be placed into any partition." << endl;
            return 0;
        }
    }

    for(int i = 0; i < maxpart+1; i++)
    {
        if(parts[i].ori_edge_num < edge_per_part && i < maxpart) {
            int rest_edge_num = edge_per_part - parts[i].ori_edge_num;
            for(int j = 0; j < rest_edge_num; j++) {
                assert(!parts[maxpart].oelist.empty() && parts[maxpart].ori_edge_num > 0);
                int mover = parts[maxpart].oelist.back();
                parts[maxpart].oelist.pop_back();
                parts[maxpart].ori_edge_num--;
                g.edges[mover].part_id = i;
                parts[i].oelist.push_back(mover);
                parts[i].ori_edge_num++;
            }
        }

        if(parts[i].ori_edge_num < edge_per_part)
            cerr << "INFO: partition " << i << " shorts for " << edge_per_part - parts[i].ori_edge_num << " edges at the end." << endl;
    }

    // output the results
    ofstream fpfile("final.spmv.part");
    fpfile << g.vertex_num / 2 << endl;
#ifdef SORT_BY_Y
    cerr << "INFO: sort by column";
    for(int i = 0; i < maxpart+1; i++)
    {
        sortEdgebyEnd(&parts[i].oelist, &g);
        cerr << ".";
    }
    cerr << endl;
#else
    cerr << "INFO: sort by row";
    for(int i = 0; i < maxpart+1; i++)
    {
        sortEdgebyBegin(&parts[i].oelist, &g);
        cerr << ".";
    }
    cerr << endl;
#endif
    for(int i = 0; i < maxpart+1; i++)
    {
        for(vector<int>::iterator it = parts[i].oelist.begin(); it != parts[i].oelist.end(); it++)
            fpfile << g.vertexes[g.edges[*it].begin - 1].ori_id + 1 << " "
                   << g.vertexes[g.edges[*it].end - 1].ori_id + 1 << endl;
        if(parts[i].ori_edge_num < edge_per_part)
        {
            if(i != maxpart)
                cerr << "ERROR: illegal results." << endl;
        }
        //fpfile << "ignore me" << endl;
    }

    /**************************************************/

    return 0;
}
