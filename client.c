//Jeff Saracco
//a test program to connect to the FLU program to test remote infections

#include "inet.h"

struct flu *new_flu_strain();


int sockfd;
struct sockaddr_in serv_addr;
char sendline[MAXLINE], recvline[MAXLINE];

struct flu
{
	int H, N; /* H and N proteins */
	int strain; /*strain, variable length*/
	int old_strain;
	int infect_day; /* day flu infected */
	int vaccine;
	struct flu *next;
};

struct flu *firstFlu;
void initialize(char **argv);// initialize the network
void rewire();
int isANeighbor(int host, int randHost);
void delNeighbor(int host, int delHost);
void clustercoeff();


int main()
{
	int host;
	        
	memset( (char *) &serv_addr, 0, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = inet_addr(SERV_HOST_ADDR);
	serv_addr.sin_port = htons(SERV_TCP_PORT);
	
	/*
	* open a TCP socket - internet stream socket
	*/
   
	if ( (sockfd = socket(AF_INET, SOCK_STREAM,0)) < 0 ) 
	{
		perror("client: can't open stream socket");
		exit(1);
	}
		
	/*
	* connect to the server
	*/

	if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr) ) < 0) 
	{	
		perror("client: can't connect to server");
		exit(2);
	}
	
	firstFlu = new_flu_strain();

	for(;;)
	{
		printf("Enter a host to infect: \n");
		scanf("%d", &host);
		
		if( send(sockfd, &host, sizeof(int), 0) == -1 )
		{
			printf("SEND ERROR!!\n");
			return 0;
		}
		send(sockfd, &firstFlu->H, sizeof(int), 0);
		send(sockfd, &firstFlu->N, sizeof(int), 0);
		send(sockfd, &firstFlu->strain, sizeof(int), 0);
		send(sockfd, &firstFlu->old_strain, sizeof(int), 0);
		send(sockfd, &firstFlu->vaccine, sizeof(int), 0);  
	
	}

}//main


struct flu *new_flu_strain()
{
	struct flu *f;

	f = (struct flu *) malloc(sizeof(struct flu));
	if (!f) 
	{
		printf("malloc of f in new_flu_strain failed!\n");
		exit(0);
	}
	f->H = 0; // make -1 = unassigned 
	f->N = 0; // make -1 = unassigned 
	f->strain = 0; // make -1 = unassigned 
	f->old_strain = 0;	
	f->vaccine = 1;
	f->infect_day = -1;
	f->next = NULL;
	return (f);
}// new_flu_strain



void initialize(char **argv)// initialize the network
{
	int max, offset;
	int i, j, k;
	int additional = 1;
	neighbors = original_number_neighbors;
	
	unit = nodes / 10;				// unit size set to tenth of total nodes
	
	for(i = 0; i < nodes; i ++)
		IDs[i] = i;
	
	if ((unit / 2) > 2999)
		unit = unit / 2;

	if ((nodes % unit) != 0)
		max = (nodes / unit) + 1; 
	else
		max = (nodes / unit); 
		
		
	level = (struct point *) malloc (sizeof(struct point) * max);
	for(i = 0; i < max; i++)\
	{
		level[i].id = (struct people *) malloc (sizeof(struct people) * unit);
		for(j = 0; j < unit; j++)
		{
			level[i].id[j].neighbor = (int *) malloc (sizeof(int) * original_number_neighbors*6);
			level[i].id[j].myNumNeighbors = 0;
			level[i].id[j].flu_strains_infected = NULL; 
			level[i].id[j].flu_strains_recovered = NULL; 
			level[i].id[j].clustco = 0;
		}
	}
	max = nodes;
	for(i = 0; i < nodes; i++) //intital  wiring loop
	{
		k = level[i/unit].id[i%unit].myNumNeighbors;
		
		for(j = 1; j <= (neighbors); j++) //set right neighbors
		{
			offset = i + j;
			if(offset >= max)
			{
				offset = (offset - (max-1)) - 1;
			}
			level[i/unit].id[i%unit].neighbor[k++] = offset;
			level[i/unit].id[i%unit].myNumNeighbors++;
		}
		
		for(j = (neighbors); j > 0; j--) //set left neighbors
		{
			offset = i - j;
			if(offset < 0)
			{
				offset += max;
			}
			level[i/unit].id[i%unit].neighbor[k++] = offset;
			level[i/unit].id[i%unit].myNumNeighbors++;
		}

	}
	
	rewire();
	clustercoeff();
	if (makepajek)
		make_pajek(argv);
	if (makedegreedist)
		make_degree_distribution(argv);
	if (tree)
		make_tree_file(argv);
	if (print_strains)
		make_strains_file(argv);

	initInfect();
	if (print_hosts)
		fprintf(host_output,"Day, hostID,infected_day,H,N,strain\n");

}//initialize



void rewire()
{
	int i,j,k,l,m, temp;
	float randnumb;
	int randHost;

	for(j = 0; j < nodes; j++)
	{
		for(k = 0; k < level[j/unit].id[j%unit].myNumNeighbors; k++)
		{
			randnumb = (rand()/(float)RAND_MAX); //random number generation
			
			if((randnumb < swnP) && ((level[j/unit].id[j%unit].neighbor[k] - j > 0) || (j >= nodes - neighbors)) && level[j/unit].id[j%unit].neighbor[k]>0)
			{
				//	printf("node %d 's %d neighbor must be rewired\n", j, level[j/unit].id[j%unit].neighbor[k]);
					randHost = nodes*(rand()/(float)RAND_MAX);				//get a random host
					while(randHost == j || isANeighbor(j,randHost))			//if it is us or one of our neighbors....
						randHost = nodes*(rand()/(float)RAND_MAX);			//	then get get a different random host
				//	printf("it will be rewired to %d\n",randHost);
					delNeighbor(level[j/unit].id[j%unit].neighbor[k],j);	//delete the host from the original neighbors array
					level[j/unit].id[j%unit].neighbor[k]=0-randHost;		//set the neighbor to the negative value of the random host, negative so that we know it has been rewired already

					//add the host to the end of the neighbor hosts array
					level[randHost/unit].id[randHost%unit].neighbor[level[randHost/unit].id[randHost%unit].myNumNeighbors++]=0-j; 
					
			}
		}
	}
	
	
}//rewire

int isANeighbor(int host, int randHost)
{
	int i;
	for(i = 0; i < level[host/unit].id[host%unit].myNumNeighbors; i++)
	{
		if(abs(level[host/unit].id[host%unit].neighbor[i]) == randHost)
			return 1;
	}
	
	return 0;
}//isANeighbor

void delNeighbor(int host, int delHost)
{
	int i,j;
	
	for(i = 0; i < level[host/unit].id[host%unit].myNumNeighbors; i++)
	{
		if(abs(level[host/unit].id[host%unit].neighbor[i]) == delHost)
		{
			for(j = i; j < level[host/unit].id[host%unit].myNumNeighbors-1 ; j++)
				level[host/unit].id[host%unit].neighbor[j] = level[host/unit].id[host%unit].neighbor[j+1];
			break;
		}
	}
		
	level[host/unit].id[host%unit].myNumNeighbors--;

}//delNeighbor

void clustercoeff()
{
	int k, l, m, n;
	int nb;
	int numpal, count;
	double clust;
	
	for(k = 0; k < nodes; k++) 
	{
		numpal = level[k/unit].id[k%unit].myNumNeighbors;
		count = 0;
		// go through all of hosts neighbors
		// go through hosts' neighbors' neighbors each
		// determine if hosts nieghbors are neighbors of each other
		
		for(l = 0; l < level[k/unit].id[k%unit].myNumNeighbors; l++)
		{ /* counts the number of connections between each hosts' neighbors*/
			nb = abs(level[k/unit].id[k%unit].neighbor[l]);
		
			for(m = 0; m < level[nb/unit].id[nb%unit].myNumNeighbors; m++)
			{
				for(n = 0; n < level[k/unit].id[k%unit].myNumNeighbors; n++)
				{
					if(abs(level[nb/unit].id[nb%unit].neighbor[m]) == abs(level[k/unit].id[k%unit].neighbor[n])) 
					{
						count++; /*number of connections between that hosts' neighbors*/
					}
				}
			}
		}
		if(numpal*(numpal-1) != 0) 
		{
			/* calculates clustering coefficient for that host*/
			clust = count/((double)(numpal*(numpal-1)));
			level[k/unit].id[k%unit].clustco = clust;
		}
	}
}
