#include "repartition_cusha.h"

#define WEIGHT_GRAPH

int main(int argc, char *argv[])
{
    /**************************************************/
    // initiation
    /**************************************************/
    if(argc != 5)
    {
        cerr << "Usage: repartition_cusha graph_file partition_file map_file edge_per_part" << endl;
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
    int tmp;

    graph g;
    gfile >> g.vertex_num >> g.edge_num;
#ifdef WEIGHT_GRAPH
    gfile >> tmp; // read graph feature (001)
#endif

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
    pfile.close();

    // read vertex map
    int ori_vertex_num;
    vmapfile >> ori_vertex_num;
    int * vmap = new int[ori_vertex_num + 1];
    int * vendmap = new int[ori_vertex_num + 1];
    for(int i = 0; i < ori_vertex_num+1; i++) {
        vmapfile >> vmap[i];
        vmapfile >> vendmap[i];
    }
    vmapfile >> tmp;
    assert(!vmapfile.good());
    vmapfile.close();

    // output the original vertex partition results
    ofstream ovpfile("ori.vertex.part");
    ovpfile << ori_vertex_num << endl;
    for(int i = 0; i < ori_vertex_num; i++) {
        int pid = i / ((ori_vertex_num + maxpart) / (maxpart + 1));
        assert(pid <= maxpart);
        ovpfile << i << " " << pid << endl;
    }
    cout << "original vertex output completes." << endl;

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
        parts[i].ori_edge_num = 0;
        parts[i].ovlist.clear();
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
#ifdef WEIGHT_GRAPH
            tmpstream >> tmp; // read weight
#endif

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

    int max_vertex_num = 0;
    for(int i = 0; i < maxpart+1; i++)
    {
        if(parts[i].vertex_num > max_vertex_num)
            max_vertex_num = parts[i].vertex_num;
        //cout << "Part " << parts[i].id << ": vertex " << parts[i].vertex_num << ", edge " << parts[i].edge_num << endl;
    }
    cout << "Balance factor: " << (double)(max_vertex_num * (maxpart+1)) / (double)g.vertex_num << endl;
    //for(vector<vertex>::iterator it = parts[331].vlist.begin(); it != parts[331].vlist.end(); it++)
    //    cout << it->id << " ";

    /**************************************************/


    /**************************************************/
    // compute score for original vertex to facilitate partition
    /**************************************************/
    original_vertex * ov = new original_vertex[ori_vertex_num];
    int total_ori_edge_num = 0;
    for(int i = 0; i < ori_vertex_num; i++)
    {
        ov[i].begin = vmap[i] - 1;
        ov[i].end = vendmap[i] - 1;
        assert(ov[i].end >= ov[i].begin);
        ov[i].plist.clear();
        ov[i].part_id = -1;
        total_ori_edge_num += ov[i].end - ov[i].begin;

        // full score for this vertex, * 2 since cutting edge scores 1
        int fullscore = (ov[i].end - ov[i].begin) * 2;
        vector<partition_score>::iterator it;

        // look up and insert the partition of this generated vertex in plist
        for(int j = ov[i].begin; j < ov[i].end; j++)
        {
            if(g.vertexes[j].oedge_id < 0)
                cerr << "ERROR: " << i << " " << j << " " << ov[i].end << endl;
            updatePartScore(&ov[i].plist, &g, &g.vertexes[j], true);
        }

        // sort plist
        if(ov[i].plist.size() > 1)
            sortPartScore(&ov[i].plist);

        int total_score = 0;
        for(it = ov[i].plist.begin(); it != ov[i].plist.end(); it++)
        {
            total_score += it->score;
            it->score -= fullscore;
        }
        //assert(total_score == fullscore);

        // add this vertex to the waiting list of partition
        for(it = ov[i].plist.begin(); it != ov[i].plist.end(); it++)
        {
            partition_score * ps = new partition_score;
            ps->id = i;
            ps->score = it->score;
            parts[it->id].waitinglist.insert(parts[it->id].waitinglist.end(), *ps);
        }

        //if(i == 0)
        //if(ov[i].plist.size() > 1)
        //{
        //    for(it = ov[i].plist.begin(); it != ov[i].plist.end(); it++)
        //        cout << it->id << " - " << it->score << "; ";
        //    cout << endl;
        //}
    }
    edge_per_part = (total_ori_edge_num + maxpart) / (maxpart + 1);
    cout << "Total original edge number: " << total_ori_edge_num << ", total partition number: " << maxpart+1 << ", edge per partition: " << edge_per_part << endl;

    /**************************************************/


    /**************************************************/
    // select original vertex for each partition
    /**************************************************/

    // sort the waiting list of each partition
    vector<int> remainingpartlist;
    remainingpartlist.clear();
    for(int i = 0; i < maxpart+1; i++)
    {
        sortPartScore(&parts[i].waitinglist);
        //if(parts[i].waitinglist.size() < vertex_per_part)
        //    cerr << "WARN: partition " << i << " doesn't have enough candidates in the waiting list."<< endl;

        //if(i == 331)
        //{
        //    cout << parts[i].waitinglist.size() << endl;
        //    for(vector<partition_score>::iterator it = parts[i].waitinglist.begin(); it != parts[i].waitinglist.end(); it++)
        //        cout << it->id << " - " << it->score << "; ";
        //    cout << endl;
        //}

        vector<partition_score>::iterator it;
        for(it = parts[i].waitinglist.begin(); it != parts[i].waitinglist.end(); it++)
        {
            assert(it->id < ori_vertex_num);
            if(ov[it->id].part_id >= 0)
                continue;
            if(parts[i].ori_edge_num >= edge_per_part && ov[it->id].end - ov[it->id].begin > 0)
                continue;
            ov[it->id].part_id = i;
            parts[i].ovlist.push_back(it->id);
            parts[i].ori_vertex_num++;
            parts[i].ori_edge_num += ov[it->id].end - ov[it->id].begin;
            //if(parts[i].ori_edge_num >= edge_per_part)
            //    break;
        }

        if(parts[i].ori_edge_num < edge_per_part)
        {
            cout << "Partition " << i << " shorts for " << edge_per_part - parts[i].ori_edge_num << " edges in the 1st round."<< endl;
            remainingpartlist.push_back(i);
        }
        else
            cout << "Partition " << i << " has " << parts[i].ori_edge_num << " edges in the 1st round."<< endl;
    }

    // find vertices which do not have home yet
    int rr_pid = 0;
    for(int i = 0; i < ori_vertex_num; i++)
    {
        if(ov[i].part_id >= 0)
            continue;

        // select partition from the remainingpartlist
        //cout << "Vertex " << i << " doesn't find a home in the 1st round."<< endl;
        assert(remainingpartlist.size() > 0);
        if(ov[i].end - ov[i].begin == 0)
        {
            // round robin
            parts[rr_pid].ovlist.push_back(i);
            parts[rr_pid].ori_vertex_num++;
            ov[i].part_id = rr_pid;
            rr_pid = (rr_pid + 1) % (maxpart + 1);
        }
        else
        {
            vector<int>::iterator it;
            it = remainingpartlist.begin();
            //assert(parts[*it].ori_vertex_num < vertex_per_part);
            assert(parts[*it].ori_edge_num < edge_per_part);
            parts[*it].ovlist.push_back(i);
            parts[*it].ori_vertex_num++;
            parts[*it].ori_edge_num += ov[i].end - ov[i].begin;
            ov[i].part_id = *it;
            if(parts[*it].ori_edge_num >= edge_per_part)
                remainingpartlist.erase(it);
        }
    }

    for(vector<int>::iterator it = remainingpartlist.begin(); it != remainingpartlist.end(); it++)
        cerr << "WARN: partition " << *it << " shorts for " << edge_per_part - parts[*it].ori_edge_num << " edges at the end." << endl;

    // output the results
    ofstream fpfile("final.part");
    //fpfile << (maxpart+1) * vertex_per_part << endl;
    fpfile << ori_vertex_num << endl;
    for(int i = 0; i < maxpart+1; i++)
    {
        for(vector<int>::iterator it = parts[i].ovlist.begin(); it != parts[i].ovlist.end(); it++)
            fpfile << *it << endl;
    }

    // output the vertex partition results
    ofstream fvpfile("final.vertex.part");
    fvpfile << ori_vertex_num << endl;
    for(int i = 0; i < ori_vertex_num; i++)
        fvpfile << i << " " << ov[i].part_id << endl;

    /**************************************************/

    return 0;
}
