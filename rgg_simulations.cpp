#include<bits/stdc++.h>
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
#define NODES 100
/*First run directed diffusion and then our_algo*/
using namespace std;

double prob=1.0;
vector<bool> sources(NODES,0);
ofstream output;

void getrand(int m,int n) //you want m distinct random numbers from first n numbers(1 to n)
{
    for(int i=0;i<NODES;i++)
        sources[i]=0;
    map<int,bool> seen;
    for(int i=0;i<m;i++)
    {
        bool flag=false;
        while(!flag)
        {
            int r=rand()%n;
            if(seen.find(r+1)==seen.end())
            {
                seen[r+1]=true;
                sources[r+1]=true;
                flag=true;
            }
        }
    }
}

unsigned long long getTimediff(struct timeval old_time)
{
	struct timeval now;
    gettimeofday(&now, NULL);

    unsigned long long ret = 1000 * (now.tv_sec - old_time.tv_sec) + (now.tv_usec - old_time.tv_usec) / 1000;

    return ret;
}



bool isgreater(float a, float b)
{
	return *(int *)&a>*(int *)&b;
}

bool islesser(float a, float b)
{
	return *(int *)&a<*(int *)&b;
}

bool isequal(float a,float b)
{
	return fabs(a-b)<=DELTA;
}
/*
Node having x and y axis. With buffer and temporary buffer. Adjacency list as a vector of node pointer.
At the start of each time slot each node selects a random packet from its buffer to send to one of its neighbor. So at start of timeslot
only buffer is full and content of temp_buf is transferred to buffer.
we have taken two buffer, since we assume all nodes send data simultaneously therefore whatever is in buffer at start of timeslot is
maintained as the number of packets in buffer. And in the given time slot node accumulate data in temp_buf. At the end of time slot
when all nodes have sent data. Content of temp_buf is Transfered to buffer.
hop_dist contains distance from sink
*/
class Node
{
public:
	float x,y;
	int id;
	vector <int> buffer;
	vector <int> temp_buf;
	vector <int> temp_time;
	vector <int> rdbuffer;
	vector <int> temp_rdbuf;
	vector <int> round;
	vector <int> adj1;
	vector <int> time_gen;
	vector <Node *> adj;
	int hop_dist;
	Node()
	{
		x = 0;
		y = 0;
		hop_dist = 0;
	}
	Node(float x1, float y1):x(x1),y(y1)
	{

	}
	void display()
	{
		cout<<x<<","<<y<<endl;
	}
	void addNode(Node *ngb)
	{
		adj.push_back(ngb);
		adj1.push_back((*ngb).id);
	}
};

bool canReceive(Node *vertex,int id,int sender, vector<bool> alohaP)//checks whether a node can receive a packet or not
{
    bool can=true;
    if(!alohaP[id])
    {
        for(int i=0;i<vertex[id].adj1.size();i++)
        {
            int x=vertex[id].adj1[i];
            //If any of the neighbors of receiving node i is going to transmit then it results in interference
            if(x!=0 && alohaP[x] && x!=sender)//If the node is not sink and it is going to transmit then it will result in interference
            {
                can=false;
                break;
            }
        }
    }
    else //If the node is going to transmit then it cannot receive
    {
        can=false;
    }
    //cout<<can<<" ";
    return can;
}

/*
quickselect algorithm for median finding in kd tree. For building kd tree
*/
Node quickSelect(Node *a,int value,int rank,bool isx)
{
	int r,j=0,k=0;
	r=rand()%value;
	Node *left,*right;
	left=new Node[value];
	right=new Node[value];
	if(isx)
	{
		for(int i=0;i<value;i++)
		{
			if((*(int *)&a[r].x)>(*(int *)&a[i].x))
			{
				left[j++]=a[i];
			}
			else if(r!=i)
			{
				right[k++]=a[i];
			}
		}
	}
	else
	{
		for(int i=0;i<value;i++)
		{
			if((*(int *)&a[r].y)>(*(int *)&a[i].y))
			{
				left[j++]=a[i];
			}
			else if(r!=i)
			{
				right[k++]=a[i];
			}
		}
	}
	if(j==rank-1)
	{
		return a[r];
	}
	else if(j<rank-1)
	{
		return quickSelect(right,k,rank-j-1,isx);
	}
	return quickSelect(left,j,rank,isx);
}

// class for Kdnode
class Kdnode
{
public:
	Kdnode *left,*right;
	Node *l;
	float split;
	Kdnode()
	{
		left=NULL;
		right=NULL;
		split=0;
		l=NULL;
	}
};
/*
Simple kd tree algorithm using quickselect to build a network
*/
Kdnode * buildKd(Node *a,int no,bool isx)
{
	if(no==0)
	{
		return NULL;
	}
	if(no==1)
	{
		Kdnode *kd=new Kdnode;
		kd->l=a;
		kd->split=isx ? a->x : a->y;
		return kd;
	}
	int j=0,k=0;
	Node *left,*right;
	left=new Node[no];
	right=new Node[no];
	Node median;
	median=quickSelect(a,no,(no+1)/2,isx);
	if(isx)
	{
		for(int i=0;i<no;i++)
		{
			if((*(int *)&a[i].x)<=(*(int *)&median.x))
			{
				left[j++]=a[i];
			}
			else
			{
				right[k++]=a[i];
			}
		}
	}
	else
	{
		for(int i=0;i<no;i++)
		{
			if((*(int *)&a[i].y)<=(*(int *)&median.y))
			{
				left[j++]=a[i];
			}
			else
			{
				right[k++]=a[i];
			}
		}
	}
	Kdnode *kd1=new Kdnode;
	kd1->left=buildKd(left,j,!isx);
	kd1->right=buildKd(right,k,!isx);
	kd1->l=a;
	kd1->split=isx ? median.x : median.y;
	return kd1;
}
/*
queryKd gives result of range query i.e. which nodes are in the range of the current node.
*/
vector <Node *> queryKd(Kdnode *tree,float x1,float y1,float x2, float y2,bool isx)
{
	if(tree==NULL)
	{
		vector <Node *> v;
		return v;
	}
	else if(tree->left==NULL&&tree->right==NULL)
	{
		if(tree->l==NULL)
		{
			cout<<"display some error"<<endl;   //TODO:remove this condition
		}
		vector <Node *> v;
		if((isgreater(tree->l->x,x1)||isequal(tree->l->x,x1)) && (islesser(tree->l->x,x2)||isequal(tree->l->x,x2)) &&
		   (isgreater(tree->l->y,y1)||isequal(tree->l->y,y1)) && (islesser(tree->l->y,y2)||isequal(tree->l->y,y2)))
		{
			v.push_back(tree->l);
		}
		return v;
	}
	if(isx)
	{
		if(islesser(x2,tree->split)||isequal(x2,tree->split))
		{
			return queryKd(tree->left, x1, y1, x2, y2, !isx);
		}
		else if(isgreater(x1,tree->split))
		{
			return queryKd(tree->right, x1, y1, x2, y2, !isx);
		}
		else
		{
			vector <Node *> v1=queryKd(tree->left, x1, y1, x2, y2, !isx);
			vector <Node *> v2=queryKd(tree->right, x1, y1, x2, y2, !isx);
			v1.insert(v1.end(),v2.begin(),v2.end());
			return v1;
		}
	}
	else
	{
		if(islesser(y2,tree->split)||isequal(y2,tree->split))
		{
			return queryKd(tree->left, x1, y1, x2, y2, !isx);
		}
		else if(isgreater(y1,tree->split))
		{
			return queryKd(tree->right, x1, y1, x2, y2, !isx);
		}
		else
		{
			vector <Node *> v1=queryKd(tree->left, x1, y1, x2, y2, !isx);
			vector <Node *> v2=queryKd(tree->right, x1, y1, x2, y2, !isx);
			v1.insert(v1.end(),v2.begin(),v2.end());
			return v1;
		}
	}
}

float distance(float x1,float y1,float x2,float y2)
{
	return (sqrt(pow(x2-x1,2)+pow(y2-y1,2)));
}

/*
graph is created using buildKd and queryKd. radius is the connectivity radius of each node
*/
void createGraph(Node *vertex,int n)
{

	float radius=0.5;
	Kdnode *kdtree= buildKd(vertex,n,true);
	for(int i=0;i<n;i++)
	{
		vector <Node *> v=queryKd(kdtree,vertex[i].x-radius,vertex[i].y-radius,vertex[i].x+radius,vertex[i].y+radius,true);

		for(vector <Node *>::iterator j=v.begin();j!=v.end();j++)
		{
			float d=distance(vertex[i].x,vertex[i].y,(*j)->x,(*j)->y);
			if(islesser(d,radius)||isequal(d,radius))
			{
				vertex[i].addNode(*j);
			}
		}
	}
}

/*
bfs algorithm to check connectivity network. Also hop_dist is updated here.
*/
bool bfs(Node *vertex,const int n,int source)
{
	bool *visited=new bool[n];
	int count=0;
	for(int i=0;i<n;i++)
	{
		visited[i]=false;
	}
	queue <pair<int,int> > q;						//to set queue
	q.push(make_pair(vertex[source].id,0));			// sink and hop distance of 0

	while(!q.empty())
	{
		pair<int,int> pr=q.front();					//from pair we are getting first element and storing in current
		int current=pr.first;						//sink
		q.pop();


		if(!visited[current])
		{
			vertex[current].hop_dist=pr.second;
			visited[current]=true;
			count++;
			for(vector <Node *>::iterator j=vertex[current].adj.begin();j!=vertex[current].adj.end();j++)
			{
				if(!visited[(*j)->id])
				{
					q.push(make_pair((*j)->id,pr.second+1));
				}
			}
		}
	}
	return count==n;
}

void displayTree(Kdnode *kdtree,int depth=0)
{
	if(kdtree != NULL)
	{
		displayTree(kdtree->left,depth+1);
		cout<<kdtree->split<<"      "<<depth<<endl;
		displayTree(kdtree->right,depth+1);
	}
}
struct HOP 													// for average buffer size for hop distance graph
{
	float avg;
	int count;
};

/*
This is the main function where simulation is performed
*/
bool simulate(Node *vertex,int value,int sink_id,double beta)
{
    ofstream f2;
    f2.open ("rgg_transition.txt");
    //vector<double> flag_t;//to keep check of first time of round generation
    vector<bool> alohaP(NODES,0);//stores whether a node will transmit or not(depends on both p and buffer)
    vector<int> selected_neighbor(NODES,0);//stores the neighbor to which the node will transmit data
    vector<int> selected_packet(NODES,0);//stores the neighbor to which the node will transmit data
    vector<double> start_t;//time when first packet of given round was generated
    vector<double> end_t;//time when last packet of given round was sunk
    vector<double> total_t;//latency for each round
    vector<double> pkt_arrived;//packets of each round arrived
    vector<double> queue;//to keep track of max queue of each round sunk
    vector<double> queue_diff;//queue size difference between sunk rounds
    double max_avg=0.0;
    int packets_sunk=0;
    double packet_avg=0.0;
    double packet_avg_max=0.0;

    for(int i=0;i<value;i++)
    {
        vertex[i].buffer.clear();
        vertex[i].rdbuffer.clear();
        vertex[i].temp_buf.clear();
        vertex[i].temp_rdbuf.clear();
        vertex[i].time_gen.clear();
        vertex[i].temp_time.clear();
        vertex[i].round.clear();
    }
    for (int i=0;i<INT_MAX;i++)//initialise all to 0
	{
	//flag_t.push_back(0);
	start_t.push_back(0);
	end_t.push_back(0);
	total_t.push_back(0);
	pkt_arrived.push_back(0);
	queue.push_back(0);
	queue_diff.push_back(0);
    }
    int max2=0;
    double rate=beta;
    int pkt_allowed=value-1;
    double gmax,gid,gx,gy=0;//global max an dits id and coordinates
    int max_round=0;//maximum rounds generated
    int max_round_sunk=0;//maximum rounds sunk
    double check=0;//variable to check packets genrated
    double packets=0;//variable to keep count of packets collected at sink
   /* vector<double> data;
    ifstream file("output_rgg.txt");// file path to read random no's

    	while(!file.eof())
	     {
		double d;
	        file >> d;
		data.push_back(d);
	     }*/

	int count=1;										// count for number of timeslot
	int find_temp = 0;
	// cout<<"SINK IS\t\t\t"<<sink_id<<endl;					//print sink id
	set <int> s;											// sink set. All the packets that sink recovers is maintained in this
	s.insert(vertex[sink_id].id);		// inserting sink id in list of packets recovered

	vector<vector<double> > round_vertex(RMAX,vector<double>(value,0));//to find vertex id's of rounds sunk
	vector<vector<double> > trans(value,vector<double>(value,0));//to hold transition matrix
	struct timeval now;
    gettimeofday(&now, NULL);

    //cout<<"Transfering Node"<<"\t\t\t"<<"Symbol Transfered"<<"\t\t\t"<<"Transfered To"<<endl;
	//while(count!=100000)	// while no. of time steps is 100000
	while(count<=16700000)
	{
	    alohaP[0]=0;
		for(int i=0;i<value-1;i++)//for packet generation at rate beta
		{
		    if(!sources[i+1])
                continue;

            double number;
			//double v=count*(value-1)+i;
			//number=data[v];//cout<<"v:"<<number<<" ";
			number=((double)rand()/(double)RAND_MAX);
			if(number<rate)  //(no.<\beta)
			{
				int x=vertex[i+1].round.size();//if(x==0)check++;//cout<<"x:"<<x<<" " <<"s_x:"<<start_t[x]<<" ";
				if(start_t[x]==0)
				{
					start_t[x]=count;//start_t.insert(start_t.begin() + x,count);
					//flag_t.insert(flag_t.begin() + x,1);
				}
				vertex[i+1].buffer.push_back(x+1);	check++;// adding data packet in every nodes buffer
				vertex[i+1].time_gen.push_back(count);
			 	vertex[i+1].rdbuffer.push_back(vertex[i+1].id);
				vertex[i+1].round.push_back(x+1);
				vertex[i+1].round.size()>max_round?max_round=vertex[i+1].round.size():max_round=max_round;
			}
	       }

        for(int i=1;i<value;i++) //fill values in alohaP, sink never transmits
        {
            double number=((double)rand()/(double)RAND_MAX);
            if(number<prob && !vertex[i].buffer.empty())//transmits with probability p if the buffer is not empty
            {
                alohaP[i]=1;
            }
            else
                alohaP[i]=0;
        }
		/*for(int i=0;i<value-1;i++)
		{ 	//srand((unsigned)time(NULL));
			double number=((double)rand()/(double)RAND_MAX);
			//number = rand() % 100 + 1;		// generating random probability will act as tossing coin by neighbor
			//number = number/100;
			cout<<number<< "  ";
			if(number<rate)  //(no.<\beta)
			{
				int x=vertex[i+1].round.size();//if(x==0)check++;//cout<<"x:"<<x<<" " <<"s_x:"<<start_t[x]<<" ";
				if(start_t[x]==0)
				{
					start_t[x]=count;//start_t.insert(start_t.begin() + x,count);
					//flag_t.insert(flag_t.begin() + x,1);
				}
				vertex[i+1].buffer.push_back(x+1);	check++;// adding data packet in every nodes buffer
			 	vertex[i+1].round.push_back(x+1);
				vertex[i+1].round.size()>max_round?max_round=vertex[i+1].round.size():max_round=max_round;
			}
	       }*/

	       //Precompute which node is going to select which packet and to which neighbor,and if the node is selecting itself for transmission
        //then it is same as setting alohaP=0

        for(int i=1;i<value;i++)
        {   if(!vertex[i].buffer.empty())
            {
            //cout<<"i:"<<i<<endl;
                int adj_size=vertex[i].adj1.size(); //degree of node i
                //cout<<vertex[i].buffer.size();
                int pos=rand()%(vertex[i].buffer.size()); // selecting random packet from buffer
                selected_packet[i]=pos;
				int node_select; // variable to select random node from neighbors
				do
				{
					//adj list has itself as neighbor so selecting random node which is other than itself
					node_select=rand()%adj_size;
				}while(vertex[i].id==vertex[i].adj1[node_select]);
                selected_neighbor[i]=node_select;

                //check if it is a self loop then it will not cause interference
                //if(vertex[i].adj1[node_select]!=vertex[sink_id].id)
                //{
                    int degree_ngb = vertex[vertex[i].adj1[node_select]].adj1.size();   // degree of neighbor node
                    double deg_ratio = (double)(adj_size)/(degree_ngb);     // probability of sending
                    double accept_prob;
                    deg_ratio < 1 ? accept_prob = deg_ratio:accept_prob = 1;
                    if(accept_prob != 1)
                    {
                        double number=((double)rand()/(double)RAND_MAX);
                        if(number > accept_prob)
                        {
                            alohaP[i]=0;
                        }
                    }

                //}
            }
        }

		for(int i=0; i<value; i++)//for packet movement in the network
		{
			//if not sink or buffer is non empty then node will transfer packet
			if((vertex[i].id!=sink_id)&&(!vertex[i].buffer.empty()) && (alohaP[i]==1))
			{
				int adj_size=vertex[i].adj1.size(); //degree of node i
				//int pos=rand()%(vertex[i].buffer.size()); // selecting random packet from buffer
				int pos;
				pos=selected_packet[i];
				int node_select; // variable to select random node from neighbors

				/*do
				{
					//adj list has itself as neighbor so selecting random node which is other than itself
					node_select=rand()%adj_size;
				}while(vertex[i].id==vertex[i].adj[node_select]);*/
                node_select=selected_neighbor[i];
                if (canReceive(vertex,vertex[i].adj1[node_select],vertex[i].id,alohaP))
                {
                    if(vertex[i].adj1[node_select]==vertex[sink_id].id)
                        {
                            packets_sunk++;
                            packet_avg=(double)((packet_avg*(packets_sunk-1)+count-vertex[i].time_gen[pos])/(double)packets_sunk);

                            //if(vertex[i].rdbuffer[pos]==158)cout<<"yes"<<" ";//cout<<vertex[i].id<<"-"<<vertex[i].buffer[pos]<<" ";
                            s.insert(vertex[i].buffer[pos]);		// inserting packet in sink set. i.e. packet is recovered by sink
                            packets++;
                            int x= vertex[i].buffer.at(pos);//pos has actual round id
                            pkt_arrived[x-1]+=1;//cout<<vertex[i].buffer.at(pos)<<" "<<x-1<<" " <<" "<<pkt_arrived[0]<<" "<<pkt_arrived[1]<<endl;
                            //round_vertex[x-1][(vertex[i].rdbuffer[pos])]=1;//set up flag if packet from given vertex is received
                            /*if (pkt_arrived[x-1]==375.0)
                            {
                                rounds_sunk++;
                            }*/
                            if (pkt_arrived[x-1]==pkt_allowed)//all packets of given round collected
                                {
                                    end_t[x-1]=count;//cout<<end_t[x-1]<<" ";///end_t.insert(end_t.begin() + (x-1),count);
                                    max_round_sunk++;max2=1;//if(max_round_sunk==38)cout<<"latency:"<<count<<endl;
                                    //cout<<max_round_sunk<<" "<<count<<endl;
                                    //rounds_sunk++;

                                }
                            vertex[i].buffer.erase(vertex[i].buffer.begin()+pos);
                            vertex[i].time_gen.erase(vertex[i].time_gen.begin()+pos);
                            vertex[i].rdbuffer.erase(vertex[i].rdbuffer.begin()+pos);
                        }
                    else
                        {
                            /*int degree_ngb = vertex[vertex[i].adj[node_select]].adj.size();   // degree of neighbor node
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
                            }*/
                            vertex[vertex[i].adj1[node_select]].temp_buf.push_back(vertex[i].buffer.at(pos));  // pushing packet in temp buffer
                            vertex[vertex[i].adj1[node_select]].temp_time.push_back(vertex[i].time_gen.at(pos));
                            vertex[i].buffer.erase(vertex[i].buffer.begin()+pos); //erase pkt from sending nodes buffer
                            vertex[i].time_gen.erase(vertex[i].time_gen.begin()+pos);
                            vertex[vertex[i].adj1[node_select]].temp_rdbuf.push_back(vertex[i].rdbuffer.at(pos));
                            vertex[i].rdbuffer.erase(vertex[i].rdbuffer.begin()+pos); //erase pkt from sending nodes buffer
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
			vertex[i].time_gen.insert(vertex[i].time_gen.end(),vertex[i].temp_time.begin(),vertex[i].temp_time.end());
			vertex[i].temp_buf.clear();
			vertex[i].temp_time.clear();
			vertex[i].rdbuffer.insert(vertex[i].rdbuffer.end(),vertex[i].temp_rdbuf.begin(),vertex[i].temp_rdbuf.end());
			vertex[i].temp_rdbuf.clear();
			//vert_buffer[i] = vertex[i].buffer.size();
			int(vertex[i].buffer.size())>max?max=int(vertex[i].buffer.size()):max=max;
			if(max>gmax)
			{
				gmax=max;
				gid=vertex[i].id;gx=vertex[i].x;gy=vertex[i].y;
			}
		}
        double avg_queue=0.0;
        for(int i=1;i<value;i++)
        {
            avg_queue=avg_queue+(double)vertex[i].buffer.size();//Gives number of unsunk packets
        }
        avg_queue=avg_queue/99.0;
        if(avg_queue>max_avg)
        {
             max_avg=avg_queue;
             if(count>=15700000)
                return false;
        }
        if(packet_avg>packet_avg_max)
        {
            packet_avg_max=packet_avg;
            if(count>=15700000)
                return false;
        }


        //cout<<packet_avg<<"   ";
        //cout<<max_avg<<" "<<count<<endl;
        if(count%100000==0)
            output<<packet_avg<<" "<<max_avg<<" "<<count<<" "<<packets_sunk<<" "<<beta<<endl;
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
    return true;
	double mean_latency,mean_queue,std_latency,std_queue=0;
	double temp1,temp2,temp3,temp4=0;
	/*for(int i=0;i<=max_round;i++)//print latency
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
//cout<<pkt_arrived[1]<<" ";//for(int i=0;i<value;i++)cout<<round_vertex[5][i]<<" ";
/*for(int i=0;i<5;i++)
	{
	for(int j=0;j<value;j++)
		{
			if(round_vertex[i][j]!=1)
			cout<<endl<<"round:"<<i<<" "<<"vertex:"<<j<<"x:"<<vertex[j].x<<" "<<"y:"<<vertex[j].y<<"size:"<<vertex[j].adj.size();
		}
		cout<<endl;
	}*/
//for(int i=0;i<6;i++)cout<<"i:"<<i<<" "<<vertex[7].adj[i]->id<<" "<<endl;
//for(int i=0;i<8;i++)cout<<"i:"<<i<<" "<<vertex[158].adj[i]->id<<" "<<vertex[vertex[158].adj[i]->id].adj.size()<<endl;
//cout<<"our"<<" "<<rate<<" "<< prob<<" "<<packets<<" "<<max_round<<" "<<max_round_sunk<<endl;

/*Metropolis chain based transition matrix and sink probability*/
//to create transition matrix
for(int i=0;i<value;i++)
{double ownprob=0;
	for(int j=0;j<vertex[i].adj.size();j++)
	{
		int node=vertex[i].adj[j]->id;
		double size1=vertex[i].adj.size();
		double size2=vertex[node].adj.size();
		double size;
		size1>size2?size=size1:size=size2;
		if(node!=i)
		{
			double prob;
			prob=double(1/(size-1));
			trans[i][node]=prob;
			ownprob+=prob;
		}
	}
trans[i][i]=1 - ownprob;//self loop
}
for(int i=0;i<value;i++)//print transition matrix to file
{
	for(int j=0;j<value;j++)
	{
		f2<<trans[i][j]<<" ";
	}
	f2<<endl;
}
double p,q=0;//sum of prob of txn from all neighbors of sink
for(int i=0;i<value;i++)
p+=trans[0][i];
q=p-trans[0][0];
cout<<p<<" "<<q<<endl;

/*Normal Random walk based transition matrix and sink probability*/
/*for(int i=0;i<vertex[0].adj.size();i++)
{
	int node=vertex[0].adj[i]->id;
	if(node!=0)
	{
		double size=vertex[node].adj.size();
		double prob;
		prob=double(1/(size));
		tprob+=prob;
	}
}
cout<<"tprob:"<<tprob<<endl;
//to create transition matrix
for(int i=0;i<value;i++)
{
	for(int j=0;j<vertex[i].adj.size();j++)
	{
		int node=vertex[i].adj[j]->id;
		double size=vertex[i].adj.size();
		if(node!=i)
		trans[i][node]=double(1/(size-1));
	}
}
for(int i=0;i<value;i++)//print transition matrix to file
{
	for(int j=0;j<value;j++)
	{
		f2<<trans[i][j]<<" ";
	}
	f2<<endl;
}double p=0;
for(int i=0;i<value;i++)
p+=trans[0][i];
cout<<p<<endl;*/
}
/*
to run the code you have to give seed for randomness
./a.out < 10000.txt 10 or ./a.out 10 < 10000.txt -----in both cases argv[0]=./a, argv[1]=10 and because of symbol < 1000.txt is input
*/
int main(int argc, char *argv[])
{
    /* parameters needed:
    nodes, fraction of nodes acting as source, probability
    source must be at the center
    select sources randomly.
    */
	/*if (argc < 3)
	{
		cout<<"seed missing"<<endl;
		return 1;
	}*/
	//srand (atoi(argv[1]));
	//srand (time(NULL));
    ifstream file("100.txt");
    ofstream dout("100_output.txt");
    //Also change line aound 676 with change in number of nodes												//to generate random seed w.r.t time
	double rate;//=atof(argv[1]);
	//prob=atof(argv[2]);
	int value;
	file>>value;
	Node *nod,i;
	nod=new Node[value];

	for(int i=0;i<value;i++)
	{
		//scanf("%f,%f",&nod[i].x,&nod[i].y);
		double d;
		file>>d;
		nod[i].x=d;
		double e;
		file>>e;
		nod[i].y=e;
		nod[i].id=i;
	}


    output.open ("temp.txt");
	createGraph(nod,value);
	int sink_id=0;

	// cout<<"-----------------------------------------"<<sink_id<<endl;
	if(bfs(nod,value,sink_id))				// to check connectivity of network. If connected then simulate
		{
                //cout<<"connected";
                //exit(0);
                //clear all global vectors
                prob=0.1;
                double src_frac=0.1;
                int src_nodes;
                for(int i=1;i<=91;i++)
                {
                    src_nodes=-1;
                    rate=0.001;
                    bool prev_rate = false;
                    for(int j=1;j<=91;j++)
                    {
                    	prev_rate = false;
                        int temp=src_frac*(NODES-1);
                        if(temp==src_nodes)
                            continue;
                        src_nodes=temp;
                        srand(10);
                        getrand(src_nodes,NODES-1);
                        //loop which change in rate as per the return value of simulation
                        double left_lim=-0.1;
                        double right_lim=-0.1;
                        //rate=0.001;
                        while(true)
                        {
                            srand(10);
                            //testing
                            cout<<"left_lim: "<<left_lim<<" right_lim: "<<right_lim<<" rate: "<<rate<<" prob: "<<prob<<" src: "<<src_frac<<endl;
                            bool check=simulate(nod,value,sink_id,rate);
                            if(check)
                            {
                                left_lim=rate;
                                if(right_lim<0)
                                {
                                    rate=rate*2.0;
                                }
                                else
                                {
                                    rate=(left_lim+right_lim)/2.0;
                                    if(rate>=0.001)
                                    {
                                        rate=(trunc(rate*10000.0))/10000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.0001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.0001;
                                    }
                                    else if(rate>=0.0001)
                                    {
                                        rate=(trunc(rate*100000.0))/100000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.00001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.00001;
                                    }
                                   else if(rate>=0.00001)
                                    {
                                        rate=(trunc(rate*1000000.0))/1000000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.000001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.000001;
                                    }
                                   else if(rate>=0.000001)
                                    {
                                        rate=(trunc(rate*10000000.0))/10000000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.0000001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.0000001;
                                    }
                                    if ((rate>=0.00001 && rate<0.0001 && (right_lim-left_lim-0.000001)<=numeric_limits<double>::epsilon())||(rate>=0.0001 && rate<0.001 && (right_lim-left_lim-0.00001)<=numeric_limits<double>::epsilon())||(rate>=0.001 && (right_lim-left_lim-0.0001)<=numeric_limits<double>::epsilon()))
                                    {
                                        //write output and break
                                        dout<<NODES<<" "<<src_frac*(NODES-1)<<" "<<prob<<" "<<left_lim<<endl;
                                        rate = left_lim;
                                        prev_rate = true;
                                        cout<<NODES<<" "<<src_frac*(NODES-1)<<" "<<prob<<" "<<left_lim<<endl;
                                        cout<<"0"<<endl;
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                right_lim=rate;
                                if(left_lim<0)
                                {
                                    rate=rate/2.0;
                                }
                                else
                                {
                                    rate=(left_lim+right_lim)/2.0;
                                    if(rate>=0.001)
                                    {
                                        rate=(trunc(rate*10000.0))/10000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.0001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.0001;
                                    }
                                    else if(rate>=0.0001)
                                    {
                                        rate=(trunc(rate*100000.0))/100000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.00001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.00001;
                                    }
                                   else if(rate>=0.00001)
                                    {
                                        rate=(trunc(rate*1000000.0))/1000000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.000001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.000001;
                                    }
                                   else if(rate>=0.000001)
                                    {
                                        rate=(trunc(rate*10000000.0))/10000000.0;
                                        if(rate-left_lim<=numeric_limits<double>::epsilon())
                                            rate=rate+0.0000001;
                                        else if(right_lim-rate<=numeric_limits<double>::epsilon())
                                            rate=rate-0.0000001;
                                    }
                                    if ((rate>=0.00001 && rate<0.0001 && (right_lim-left_lim-0.000001)<=numeric_limits<double>::epsilon())||(rate>=0.0001 && rate<0.001 && (right_lim-left_lim-0.00001)<=numeric_limits<double>::epsilon())||(rate>=0.001 && (right_lim-left_lim-0.0001)<=numeric_limits<double>::epsilon()))
                                    {
                                        //write output and break
                                        dout<<NODES<<" "<<src_frac*(NODES-1)<<" "<<prob<<" "<<left_lim<<endl;
                                        rate = left_lim;
                                        prev_rate = true;

                                        cout<<NODES<<" "<<src_frac*(NODES-1)<<" "<<prob<<" "<<left_lim<<endl;
                                        cout<<"1"<<endl;
                                        break;
                                    }
                                }
                            }
                            if(rate<0.00001)
                            {
                                dout<<NODES<<" "<<src_frac*(NODES-1)<<" "<<prob<<" "<<"less than 0.00001";

                                cout<<"2"<<endl;
                                break;
                            }
                        }

                        src_frac+=0.01;
                    }
                    if(!prev_rate)
                    {
                    	rate = 0.001;
                    }
                    prob+=0.01;
                }



		}

	else
		cout<<"not connceted"<<endl;
	return 0;
}
