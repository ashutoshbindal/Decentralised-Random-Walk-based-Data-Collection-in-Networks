#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<vector>
#include<queue>
#include<set>
#include<sys/time.h>
#include<time.h>
#include<fstream>
#include <sstream>
#include<string>

#define PRECISION 5
#define DELTA pow(10,-PRECISION)
#define INT_MAX 9999999
#define RMAX 30000
#define SIZE 10000
#define NODES 100
using namespace std;
/*Generate the adjacency matrix for given k-regular graph using matlab createrandgraph function, then input that file in creatematrix function, change the no. of nodes in definition NODES and while running the code give a argc count of 3 and given rate as input.*/
vector<int>list;//read adjacency matrix values from matlab exported file
vector<vector<int> >adj_matrix(NODES,vector<int>(NODES,0));//adjacency matrix of our d-regular graph
vector<int>dist(NODES,0);//store the distance from sink
class Node
{
public:
	int id;
	vector <int> buffer;
	vector <int> temp_buf;
	vector <int> rdbuffer;
	vector <int> temp_rdbuf;
	vector <int> round;
	vector <int> adj;
	int hop_dist;
	Node()
	{
		hop_dist = 0;
	}
};
void creatematrix()
{
	ifstream file("deg5_200.txt");// file path to read random no's

    	while(!file.eof())//reading file into 1d array
	     {
		int d;
	        file >> d;
		list.push_back(d);
	     }
	for(int i=0;i<SIZE;i++)//mapping 1d array into 2d adjacency matrix
	{
		int rows=i/NODES;
		int cols=i%NODES;
		adj_matrix[rows][cols]=list[i];
	}
	/*for(int i=0;i<NODES;i++)//display matrix
	{
		for(int j=0;j<NODES;j++)
		{
			cout<<adj_matrix[i][j];
		}
		cout<<endl;
	}*/
}
void addnodes(Node *vertex,int id)//if adj_matrix has 1 add node id in adjacency list
{
	for(int i=0;i<NODES;i++)
	{
		if(adj_matrix[id][i]==1)
		vertex[id].adj.push_back(i);
	}
}

/*Dijkstra's algorithm to find shortest distance*/
// A utility function to find the vertex with minimum distance value, from
// the set of vertices not yet included in shortest path tree
int minDistance(vector<int> dist, bool sptSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;

   for (int v = 0; v < NODES; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;

   return min_index;
}
// A utility function to print the constructed distance array
int printSolution(vector<int> dist, int n)
{
   printf("Vertex   Distance from Source\n");
   for (int i = 0; i < NODES; i++)
      printf("%d \t\t %d\n", i, dist[i]);
}
// Funtion that implements Dijkstra's single source shortest path algorithm
// for a graph represented using adjacency matrix representation
void dijkstra(Node *vertex,vector<vector<int> > graph,vector<int> dist,int src)// The output array.  dist[i] will hold the shortest distance from src to i
{
     bool sptSet[NODES]; // sptSet[i] will true if vertex i is included in shortest
                     // path tree or shortest distance from src to i is finalized

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < NODES; i++)
        dist[i] = INT_MAX, sptSet[i] = false;

     // Distance of source vertex from itself is always 0
     dist[src] = 0;

     // Find shortest path for all vertices
     for (int count = 0; count < NODES-1; count++)
     {
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in first iteration.
       int u = minDistance(dist, sptSet);

       // Mark the picked vertex as processed
       sptSet[u] = true;

       // Update dist value of the adjacent vertices of the picked vertex.
       for (int v = 0; v < NODES; v++)

         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
                                       && dist[u]+graph[u][v] < dist[v])
            dist[v] = dist[u] + graph[u][v];
     }

        //printSolution(dist, NODES);// print the constructed distance array
	for(int i=0;i<NODES;i++)//update hop_distance from sink
	{
		vertex[i].hop_dist=dist[i];
	}
}

void creategraph(Node *vertex)
{
	for(int i=0;i<NODES;i++)//add node ids and neighbors
	{
		vertex[i].id=i;
		addnodes(vertex,i);
	}
	dijkstra(vertex,adj_matrix,dist,0);//find hop_dist
}

/*
This is the main function where simulation is performed
*/
void simulate(Node *vertex,int sink_id,double rate)
{
    //ofstream f2;
    //f2.open ("output1.csv");
    int value=NODES;//500 NODES
    vector<double> start_t;//time when first packet of given round was genearted
    vector<double> end_t;//time when last packet of given round was sunk
    vector<double> total_t;//latency for each round
    vector<double> pkt_arrived;//packets of each round arrived
    vector<double> queue;//to keep track of max queue of each round sunk
    vector<double> queue_diff;//queue size difference between sunk rounds

    for (int i=0;i<INT_MAX;i++)//initialise all to 0
	{
	start_t.push_back(0);
	end_t.push_back(0);
	total_t.push_back(0);
	pkt_arrived.push_back(0);
	queue.push_back(0);
	queue_diff.push_back(0);
        }
    int max2=0;
    int pkt_allowed=value-1;
    double gmax,gid,gx,gy=0;//global max an dits id and coordinates
    int max_round=0;//maximum rounds generated
    int max_round_sunk=0;//maximum rounds sunk
    double check=0;//variable to check packets genrated
    double packets=0;//variable to keep count of packets collected at sink
    vector<double> data;
    ifstream file("output_reg.txt");// file path to read random no's

    	while(!file.eof())
	     {
		double d;
	        file >> d;
		data.push_back(d);
	     }

	int count=0;										// count for number of timeslot
	int find_temp = 0;
	// cout<<"SINK IS\t\t\t"<<sink_id<<endl;					//print sink id
	set <int> s;											// sink set. All the packets that sink recovers is maintained in this
	s.insert(vertex[sink_id].id);		// inserting sink id in list of packets recovered

	vector<vector<double> > round_vertex(RMAX,vector<double>(value,0));//to find vertex id's of rounds sunk
	struct timeval now;
    gettimeofday(&now, NULL);

    //cout<<"Transfering Node"<<"\t\t\t"<<"Symbol Transfered"<<"\t\t\t"<<"Transfered To"<<endl;
	while(count!=100000)	// while no. of time steps is 100000
	{
		for(int i=0;i<value-1;i++)//for packet generation at rate beta
		{ 	double number;
			double v=count*(value-1)+i;
			number=data[v];//cout<<"v:"<<number<<" ";
			if(number<rate)  //(no.<\beta)
			{
				int x=vertex[i+1].round.size();//if(x==0)check++;//cout<<"x:"<<x<<" " <<"s_x:"<<start_t[x]<<" ";
				if(start_t[x]==0)
				{
					start_t[x]=count;//start_t.insert(start_t.begin() + x,count);
					//flag_t.insert(flag_t.begin() + x,1);
				}
				vertex[i+1].buffer.push_back(x+1);	check++;// adding data packet in every nodes buffer
			 	vertex[i+1].rdbuffer.push_back(vertex[i+1].id);
				vertex[i+1].round.push_back(x+1);
				vertex[i+1].round.size()>max_round?max_round=vertex[i+1].round.size():max_round=max_round;
			}
        }


		for(int i=0; i<value; i++)//for packet movement in the network
		{
			//if not sink or buffer is non empty then node will transfer packet
			if((vertex[i].id!=sink_id)&&(!vertex[i].buffer.empty()))
			{
				int adj_size=vertex[i].adj.size(); //degree of node i
				int pos=rand()%(vertex[i].buffer.size()); // selecting random packet from buffer
				int node_select; // variable to select random node from neighbors

				do
				{
					//adj list has itself as neighbor so selecting random node which is other than itself
					node_select=rand()%adj_size;
				}while(vertex[i].id==vertex[i].adj[node_select]);

				if(vertex[i].adj[node_select]==vertex[sink_id].id)
					{//if(vertex[i].rdbuffer[pos]==158)cout<<"yes"<<" ";//cout<<vertex[i].id<<"-"<<vertex[i].buffer[pos]<<" ";
						s.insert(vertex[i].buffer[pos]);		// inserting packet in sink set. i.e. packet is recovered by sink
						packets++;
						int x= vertex[i].buffer.at(pos);//pos has actual round id
						pkt_arrived[x-1]+=1;//cout<<vertex[i].buffer.at(pos)<<" "<<x-1<<" " <<" "<<pkt_arrived[0]<<" "<<pkt_arrived[1]<<endl;
						//round_vertex[x-1][(vertex[i].rdbuffer[pos])]=1;//set up flag if packet ffrom given vertex is received
						if (pkt_arrived[x-1]==pkt_allowed)//all packets of given round colleceted
							{
								end_t[x-1]=count;//cout<<end_t[x-1]<<" ";///end_t.insert(end_t.begin() + (x-1),count);
								max_round_sunk++;max2=1;//if(max_round_sunk==38)cout<<"latency:"<<count<<endl;

							}
						vertex[i].buffer.erase(vertex[i].buffer.begin()+pos);
						vertex[i].rdbuffer.erase(vertex[i].rdbuffer.begin()+pos);
					}
				else
					{
						int degree_ngb = vertex[vertex[i].adj[node_select]].adj.size();   // degree of neighbor node
						double deg_ratio = (double)(adj_size)/(degree_ngb);     // probability of sending
						double accept_prob;
						deg_ratio < 1 ? accept_prob = deg_ratio:accept_prob = 1;//acceptance probability is minimum of deg_ratio and 1
						// if sending node degree is greater than neighbor node degree then packet will always be accepted by neighbor node
						if(accept_prob == 1)
						{
							vertex[vertex[i].adj[node_select]].temp_buf.push_back(vertex[i].buffer.at(pos));
							vertex[i].buffer.erase(vertex[i].buffer.begin()+pos); //erase pkt from sending nodes buffer
							vertex[vertex[i].adj[node_select]].temp_rdbuf.push_back(vertex[i].rdbuffer.at(pos));
							vertex[i].rdbuffer.erase(vertex[i].rdbuffer.begin()+pos); //erase pkt from sending nodes buffer
						}else
						{
							double number=((double)rand()/(double)RAND_MAX);
							if(number < accept_prob)
							{
							// cout << number/100 << "  " << adj_size << " " << degree_ngb << endl;
							vertex[vertex[i].adj[node_select]].temp_buf.push_back(vertex[i].buffer.at(pos));  // pushing packet in temp buffer
							vertex[i].buffer.erase(vertex[i].buffer.begin()+pos); //erase pkt from sending nodes buffer
							vertex[vertex[i].adj[node_select]].temp_rdbuf.push_back(vertex[i].rdbuffer.at(pos));
							vertex[i].rdbuffer.erase(vertex[i].rdbuffer.begin()+pos); //erase pkt from sending nodes buffer
							}
						}
					}
				}
		}
		/*
		since we added packets in temp buffer. In this loop we put all the data from temp buffer to main buffer before the start of the
		next time slot.
		vert_buffer[i] maintains the buffer size of each node in the network
		*/
		double max=0;
		for(int i=0;i<value;i++)
		{
			//since earlier pkt was added in temp_buf now transfer pkt from temp_buf to main buffer and empty temp_buf
			vertex[i].buffer.insert(vertex[i].buffer.end(),vertex[i].temp_buf.begin(),vertex[i].temp_buf.end());
			vertex[i].temp_buf.clear();
			vertex[i].rdbuffer.insert(vertex[i].rdbuffer.end(),vertex[i].temp_rdbuf.begin(),vertex[i].temp_rdbuf.end());
			vertex[i].temp_rdbuf.clear();
			//vert_buffer[i] = vertex[i].buffer.size();
			int(vertex[i].buffer.size())>max?max=int(vertex[i].buffer.size()):max=max;
		}
		//if(max2==1)
		//queue[max_round_sunk-1]=max;
		queue[count]=max;
		max2=0;
		//cout<<"buffer:"<<vertex[193].buffer.size()<<endl;
		//f2<<count<<" "<<max<<endl;
		count++;
		//cout << count << " "<<max<<endl;
		//cout<<endl;
	}

	/*double mean_latency,mean_queue,std_latency,std_queue=0;
	double temp1,temp2,temp3,temp4=0;
	for(int i=0;i<=max_round;i++)//print latency
		{
			if(end_t[i]!=0)
			{
				total_t[i]=end_t[i]-start_t[i]+1;
				temp1+=total_t[i];//cout<<i+1<<" "<<start_t[i]<<" "<<end_t[i]<<" "<<total_t[i]<<endl;
			}//cout<<packets<<" ";//<<check<<" ";
		}mean_latency=temp1/max_round_sunk;
	for(int i=0;i<=max_round;i++)
		{
			if(end_t[i]!=0)
			{
				temp2+=(total_t[i]-mean_latency)*(total_t[i]-mean_latency);
			}//cout<<packets<<" ";//<<check<<" ";
		}std_latency=sqrt(temp2/max_round_sunk);

	for(int i=0;i<max_round_sunk-1;i++)
	{
		queue_diff[i]=queue[i+1]-queue[i];
		temp3+=queue_diff[i];//cout<<" "<<queue_diff[i]<<" ";
	}mean_queue=temp3/(max_round_sunk-1);
	for(int i=0;i<max_round_sunk-1;i++)
	{
		temp4+=(queue_diff[i]-mean_queue)*(queue_diff[i]-mean_queue);
	}std_queue=sqrt(temp4/(max_round_sunk-1));
	for(int i=5001;i<count-1;i++)
	{
		queue_diff[i]=queue[i+1]-queue[i];
		temp3+=queue_diff[i];//cout<<" "<<queue_diff[i]<<" ";
	}mean_queue=temp3/(count-1-5001);
	for(int i=5001;i<count-1;i++)
	{
		temp4+=(queue_diff[i]-mean_queue)*(queue_diff[i]-mean_queue);
	}std_queue=sqrt(temp4/(count-1-5001));*/
	//cout<<"round_generated:"<<max_round<<" " <<"round_sunk:"<<max_round_sunk<<endl;
	//cout << "time taken is " << getTimediff(now) << " ms" <<endl;			//time
	//cout<<count<<" "<<"packets="<<" " <<packets<<" "<<"mean_latency="<<" "<<mean_latency<<" "<<"std_latency="<<" "<<std_latency<<" "<<"mean_queue="<<" "<<mean_queue<<" "<<"std_queue="<<" "<<std_queue<<" "<<"done"<<endl;
	//cout<<"bottleneck id:"<<gid<<" "<<"x:"<<" "<<gx<<"y:"<<gy<<" "<<endl;
	//cout<<"pkt_allowed:"<<pkt_allowed<<" "<<"count:"<<count<<" "<<"rate"<<rate<<" "<<"packets:"<<" " <<packets<<" "<<"mean_queue="<<" "<<mean_queue<<" "<<"std_queue="<<" "<<std_queue<<" "<<"done"<<endl;
//for(int i=0;i<100;i++)cout<<pkt_arrived[i]<<" ";
//print algo rate: packets_sunk: round_generated: round_sunk:
	cout<<"our"<<" "<<rate<<" "<<packets<<" "<<max_round<<" "<<max_round_sunk<<endl;
/*for(int i=0;i<5;i++)
	{
	for(int j=0;j<value;j++)
		{
			if(round_vertex[i][j]!=1)
			cout<<endl<<"round:"<<i<<" "<<"vertex:"<<j<<"x:"<<vertex[j].x<<" "<<"y:"<<vertex[j].y<<"size:"<<vertex[j].adj.size();
		}
		cout<<endl;
	}*/

}
//to run g++ -o oureg fileneame.cpp then ./oureg 3 rate
int main(int argc, char *argv[])
{
	double rate=0.001;
	creatematrix();
	Node *nod;
	nod=new Node[NODES];
	creategraph(nod);
	//cout<<nod[1].adj.size()<<" "<<nod[0].adj[4]<<" "<<nod[0].hop_dist<<" "<<nod[1].hop_dist<<" ";
	simulate(nod,0,rate);//node vector and sink_id
return 0;
}
