/* 
A Distributed Multiple City Simulation of Small-World Networks
	each "city" resides on its own node in the cluster, the IP 
	addresses of the nodes are stored in the header file inet.h
	
	
	Originally Written by: Jeff Saracco, March 2006


*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "inet.h"



#define ver "3.0" /* Network Influenza Model (NIM) Version # */
#define first_host 0

int main( int argc, char **argv);
void runProg(char **argv);
void freeNetwork();
void listen_for_connects(int i);
void connect_cities(int i);
void get_commands(int i);
void rewire();
int isANeighbor(int host, int randHost);
void delNeighbor(int host, int delHost);
int isAForNeigh(int host, int foreigner, int city);
void delForNeigh(int host, int foreigner, int city);
void initialize(char **argv);
void readvars();
void global_init();
void global_rewire();
void print_people();
void makeGraphicFile();

struct people 
{
	struct neighbor_info *neighbor;
	int myNumNeighbors;
	int myNumForNeigh;
	int cross;
	struct flu *flu_strains_infected; /* current strains in host */
	struct flu *flu_strains_recovered; /* past strains of flu */
	double clustco;
};

struct neighbor_info
{
	int ID;
	int city;
};

struct point
{
	struct people *id;
} *level;

struct infectious
{
	int h;
	struct infectious *next;
	struct infectious *prev;
};

struct flu
{
	int H, N; /* H and N proteins */
	int strain; /*strain, variable length*/
	int old_strain;
	int infect_day; /* day flu infected */
	int vaccine;
	struct flu *next;
};

struct strain_list
{
   int first_day;
   int last_day_inf;
   int strain;
   int old_strain;
   int instance;
   int hamming_distance;
   struct strain_list *next;
};

int unit, nodes;
struct strain_list *first_strain;
struct flu *first_flu;
int neighbors;
char printinfo, vaccinate;
int day,vacc_day,net_type,num_days,num_runs,numbits_Hsubtypes,numbits_Nsubtypes, num_net_types, strat_comb, pc_comb;
int half_pop_inf,num_days_vacc, num_vacc_strat, vacc_strats[36], vacc_strat, print_average, tree, count_strains, currentvacc;
int strategy,percentage,original_number_neighbors,net_types[2], duration, peak_day,countrun = 0;
int vacc_days[100],num_pc_vacc,num_pop_size,pop_sizes[100],num_NCR,num_mutate,print_strains,tot_count_strains;
int num_days_infectious,latent_period,num_reps;
int print_hosts,makepajek,makedegreedist,num_infect,infect_H[20],infect_N[20];
int vaccine_subtype_H[3],vaccine_subtype_N[3],IDs[1000000], cross_cuts[1000000];
int num_vaccinate,num_flu_vaccine,originalstrain[20],vacstrain[20],num_p_values;
long numbits_strain;
float NCRs[100],p_values[100],mutation_rates[100],mut_rate,swnP,NCR,num_infected,peak_num = 0;
float totalinf[10000],averageinf[10000],totalrec[10000],averagerec[10000],numstrain[10000],averagetotstrain[10000];
float totalstrain[10000],averagestrain[10000],per_vaccinate,pc_vacc[100];
float totalinfections[10000],averageinfections[10000];
FILE *input, *average_flu_days, *host_output, *flu_sum_output, *tree_files, *strain_files;

int CITY_NUMBER;
int GO_AHEAD = 0, go2 = 0;

/* declarations for sockets */
int *city_array;
int num_Cities = 4;
int this_city;
int SERV_TCP_PORT = 5529;
int num_connected = 0;
int n_cross_city_edges = 10;
int num_for_neigh;
int *foreign_array;
float swnP_city = 0.5;

struct timeval tv;
fd_set readfds;


int main(int argc, char **argv)
{

	char s[120];
	int i;
	pthread_t doWork, connects[num_Cities], listens[num_Cities], commands[num_Cities];
	
	if (argc != 3) 
	{
		printf ("Try: %s RandNumberSeed CityNumber\n",argv[0]);
		exit(0);
	}
	srand(25 * atoi(argv[1]));
	
	CITY_NUMBER = atoi(argv[2]);
	
	//srand((unsigned)time(NULL));
	strcpy(s,"sum_file");
	strcat(s,argv[1]);
	strcat(s," city ");
	strcat(s,argv[2]);
	strcat(s,".csv");
	if ((flu_sum_output = fopen(s, "w")) == NULL) 
	{
		printf("failed to open flu summary output file: be sure it's closed in Excel\n");
		exit(0);
	}
	
	pthread_create(&doWork, NULL, (void *) runProg, (void *) &argv);
	
	for(i = CITY_NUMBER; i <= num_Cities; i++) //create listening threads
	{
		if( i != CITY_NUMBER )
			pthread_create(&listens[i-1], NULL, (void *) listen_for_connects, i);
	}
	sleep(5);
	SERV_TCP_PORT = 5529;
	for(i = CITY_NUMBER-1; i > 0; i--) //create talking threads
	{
		if( i != CITY_NUMBER )
			pthread_create(&connects[i-1], NULL, (void *) connect_cities, i);
	}
	
	while(num_connected < num_Cities-1)
	{}
		
	for(i = 1; i <= num_Cities; i++) //create commands listening threads
	{
		if( i != CITY_NUMBER )
			pthread_create(&commands[i-1], NULL, (void *) get_commands, i);	
	}
	
	for(i = 1; i <= num_Cities; i++) //clean up listening threads
	{
		if( i != CITY_NUMBER )
			pthread_join(listens[i-1],NULL);
	}
	
	for(i = 1; i <= num_Cities; i++) //clean up talking threads
	{
		if( i != CITY_NUMBER )
			pthread_join(connects[i-1],NULL);
	}
	
	pthread_join(doWork,NULL);
	
}

void listen_for_connects(int i)//listen_for_connects
{
	int host,numbytes;
	int city = i;
	int city_sock;
	int clilen,childpid, n;
	struct sockaddr_in cli_addr, serv_addr, serv_addr2;
	char recvline[MAXLINE];

	tv.tv_sec = 2;
	tv.tv_usec = 500000;
	
	struct flu *firstFlu = (struct flu *) malloc(sizeof(struct flu));
	struct flu *f = (struct flu *) malloc(sizeof(struct flu));
	
	FD_ZERO(&readfds);
	
	/*
    * open a TCP socket - Internet stream socket
    */
   if ( (this_city = socket(AF_INET, SOCK_STREAM,0)) < 0 ) 
   {
      perror("server: can't open stream socket");
      exit(1);
   }
   
   
	/*
    * bind our local address so that the client can send us
    */
	memset((char *) &serv_addr, 0, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
	serv_addr.sin_port = htons(SERV_TCP_PORT++);

	if ( bind(this_city, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
	{
		perror("server: can't bind local address");
		exit(1);
	}

	printf("waiting for City %d to connect.....\n", i);
	listen(this_city, 15);
  while(1) 
   {
		/*
		* wait for a connection from a client process
		* - concurrent server -
		*/
		clilen = sizeof(cli_addr);
		city_sock = accept(this_city, (struct sockaddr *) &cli_addr, &clilen); 
      
		if (city_sock < 0 ) 
		{
			if (errno == EINTR)
				continue; /*try again*/
			perror("server: accept error");
			exit(1);
		}
	
		city_array[city-1] = city_sock;
		FD_SET(city_array[city-1], &readfds);

		
		printf("city %d connected\n",city);
		num_connected++;
		
	}//while 1
	
}//listen_for_connects


void connect_cities(int city)
{
	int i = 1, toReceive;
	
	int clilen,childpid, n;
	struct sockaddr_in cli_addr, serv_addr, serv_addr2;
	char recvline[MAXLINE];
	
	struct flu *f = first_flu;

	memset( (char *) &serv_addr, 0, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	
	printf("Trying to connect to city %d...\n",city);
	
	switch(city)
	{
		case 1:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_1);
		break;
		
		case 2:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_2);
		break;
		
		case 3:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_3);
		break;
		
		case 4:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_4);
		break;
		
		case 5:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_5);
		break;
		
		case 6:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_6);
		break;
		
		case 7:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_7);
		break;
		
		case 8:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_8);
		break;
		
		case 9:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_9);
		break;
		
		case 10:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_10);
		break;
		
		case 11:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_11);
		break;
		
		case 12:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_12);
		break;
		
		case 13:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_13);
		break;
		
		case 14:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_14);
		break;
		
		case 15:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_15);
		break;
		
		case 16:
			serv_addr.sin_addr.s_addr = inet_addr(CITY_16);
		break;
		
		default:
		break;
	}//switch
	
	serv_addr.sin_port = htons(SERV_TCP_PORT++);
	
	/*
	* open a TCP socket - internet stream socket
	*/
   
	if ( (city_array[city-1] = socket(AF_INET, SOCK_STREAM,0)) < 0 ) 
	{
		perror("client: can't open stream socket");
		exit(1);
	}
		
	/*
	* connect to the server
	*/

	if (connect(city_array[city-1], (struct sockaddr *) &serv_addr, sizeof(serv_addr) ) < 0) 
	{	
		perror("client: can't connect to server");
		exit(2);
	}

	printf("CONNECTED to city %d\n",city);
	num_connected++;

}//connect_cities

void get_commands(int i)
{
	int city = i;
	int toReceive, temp;
	int foreigner, host, host_city;
	
	while( day < num_days )
	{
		if(recv(city_array[city-1], &toReceive, sizeof(int), 0) < 0)
			break;
		
	//	recv(city_array[city-1], &toReceive, sizeof(int), 0);
		
		switch(toReceive)
		{
			case -12:
				go2 = 1;
				break;
		
			case -11:
				break;
				
			case -10:
				GO_AHEAD = 1;
				break;
			
			case -9:// global rewiring code
				recv(city_array[city-1], &foreigner, sizeof(host), 0);
				recv(city_array[city-1], &host, sizeof(host), 0);
				recv(city_array[city-1], &host_city, sizeof(host), 0);
				printf("received host: %d to connect to %d from city %d\n",host,foreigner,host_city);
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors].ID = foreigner;
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors++].city = host_city;
				level[host/unit].id[host%unit].myNumForNeigh++;
				if(level[host/unit].id[host%unit].myNumForNeigh == 1)
					foreign_array[num_for_neigh++] = host;
				num_connected++;
				toReceive = -11;
				break;
			
			case -6:
				recv(city_array[city-1], &foreigner, sizeof(host), 0);
				recv(city_array[city-1], &host, sizeof(host), 0);
				recv(city_array[city-1], &host_city, sizeof(host), 0);
		//		printf("received deletion of host %d's neighbor %d, from city %d\n",host, foreigner, city);
				delForNeigh(host, foreigner, city);
				toReceive = -11;				
				break;
				
			default:
		//		printf("ERROR, Unknown Function %d from city %d\n",toReceive,city);
				toReceive = -11;
				break;		
		}
		
		temp=toReceive;
	}
	
	close(city_array[city-1]);

}//get_commands

void runProg(char **argv)
{
	int a, b, c, f, l, k, j, g, y;
	char average_flu_day[120],s[120];
	int si;
	int numinfhosts, numrechosts, numinfections;
	struct flu *recovereds;
	int toSend = 512;
	
	readvars();
	city_array = (int *) malloc(sizeof(int) * num_Cities);
	foreign_array = (int *) malloc(sizeof(int) * n_cross_city_edges * num_Cities); //worst case scenario
	fprintf(flu_sum_output,"\nNet_type,swnP,vacc_strats,per_vaccs,popsize,NCR,mut_rate,num_inf,tot_inf,duration,peak_inf,peak_day,half_pop_inf_day,# strains\n");
	if(!vaccinate) {
		printf("WARNING: Not vaccinating but summary output will report vacc_strat!\n");
		printf("\tAlthough slower, consider vaccinating but setting proportion \n\tof pop to vacc = 0.0%\n\n");
	}
	
	swnP = p_values[0];
	nodes = pop_sizes[0];
	//	printf("Population: %d\n",nodes);
	while(num_connected < num_Cities-1)
	{}
	initialize(argv);  /* set up model and do statistics */
	printf("done setting up network and rewiring\n");
	printf("NORMAL EXECUTION\n");
	makeGraphicFile();
	exit(0);
}

void freeNetwork()
{
	int i, j;
	struct flu *temp;
	
	for(i = 0; i < nodes/unit; i ++)
	{
		for(j = 0; j < unit; j ++)
		{
			free(level[i].id[j].neighbor);
			temp = level[i].id[j].flu_strains_infected;
			while(temp)
			{
				temp = temp->next;
				free(level[i].id[j].flu_strains_infected);
				level[i].id[j].flu_strains_infected = temp;
			}
			temp = level[i].id[j].flu_strains_recovered;
			while(temp)
			{
				temp = temp->next;
				free(level[i].id[j].flu_strains_recovered);
				level[i].id[j].flu_strains_recovered= temp;
			}
			free(temp);
		}
		free(level[i].id);
	}
	free(level);
}

void initialize(char **argv)// initialize the network
{
	int max, offset;
	int i, j, k;
	int additional = 1;
	neighbors = original_number_neighbors;
	int totForNeigh = 0, indFN = 0;
	
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
			level[i].id[j].neighbor = (struct neighbor_info *) malloc (sizeof(struct neighbor_info) * original_number_neighbors*6);
			level[i].id[j].myNumNeighbors = 0;
			level[i].id[j].myNumForNeigh = 0;
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
			level[i/unit].id[i%unit].neighbor[k].ID = offset;
			level[i/unit].id[i%unit].neighbor[k++].city = CITY_NUMBER;
			level[i/unit].id[i%unit].myNumNeighbors++;
		}
		
		for(j = (neighbors); j > 0; j--) //set left neighbors
		{
			offset = i - j;
			if(offset < 0)
			{
				offset += max;
			}
			level[i/unit].id[i%unit].neighbor[k].ID = offset;
			level[i/unit].id[i%unit].neighbor[k++].city = CITY_NUMBER;
			level[i/unit].id[i%unit].myNumNeighbors++;
		}

	}
	sleep(1);
	rewire();
	GO_AHEAD = -10;
	if(CITY_NUMBER == num_Cities)
		send(city_array[0],&GO_AHEAD,sizeof(int),0);	

	while(GO_AHEAD != 1)
	{}

	if(CITY_NUMBER == 1)
	{	
		GO_AHEAD = -10;
		for(i = 1; i <= 3; i++)
			send(city_array[i],&GO_AHEAD,sizeof(int),0);
		printf("SENT GO AHEADS TO ALL CITIES\n");
	}
	
	num_connected = 0;
	if(CITY_NUMBER%2 != 0)
		global_init();
	printf("OUT OF GLOBAL INIT\n\n\n");
	
	GO_AHEAD = -10;
	if(CITY_NUMBER == num_Cities)
		send(city_array[0],&GO_AHEAD,sizeof(int),0);	

	while(GO_AHEAD != 1)
	{}
	printf("ALL CITIES GLOBALLY INIT'D\n");
	if(CITY_NUMBER == 1)
	{	
		GO_AHEAD = -10;
		for(i = 1; i <= 3; i++)
			send(city_array[i],&GO_AHEAD,sizeof(int),0);
		printf("SENT GO AHEADS TO ALL CITIES\n");
	}
	
	
	if(CITY_NUMBER == 1)
		go2 = 1;
	while(go2 != 1)
	{}
		
	global_rewire();	
	printf("OUT OF GLOBAL REWIRE\n\n\n");
	go2 = -12;
	if(CITY_NUMBER != num_Cities)
		send(city_array[CITY_NUMBER], &go2, sizeof(int),0);
	
	GO_AHEAD = -10;
	if(CITY_NUMBER == num_Cities)
		send(city_array[0],&GO_AHEAD,sizeof(int),0);

	/*if(CITY_NUMBER = num_Cities - 1)
	{	
		sleep(1);
		send(city_array[CITY_NUMBER],&GO_AHEAD,sizeof(int),0);
	}*/
	
	while(GO_AHEAD != 1)
	{}
	printf("ALL CITIES GLOBALLY REWIRED\n");
	if(CITY_NUMBER == 1)
	{	
		GO_AHEAD = -10;
		for(i = 1; i <= 3; i++)
			send(city_array[i],&GO_AHEAD,sizeof(int),0);
		printf("SENT GO AHEADS TO ALL CITIES\n");
	}

//	print_people();
	
	printf("CITY %d has %d foreign neighbors %d\n",CITY_NUMBER,totForNeigh,indFN);
	//sleep(10);
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
			
			if((randnumb < swnP) && ((level[j/unit].id[j%unit].neighbor[k].ID - j > 0) || (j >= nodes - neighbors)) && level[j/unit].id[j%unit].neighbor[k].ID>0)
			{
				//	printf("node %d 's %d neighbor must be rewired\n", j, level[j/unit].id[j%unit].neighbor[k]);
					randHost = nodes*(rand()/(float)RAND_MAX);				//get a random host from this city
					while(randHost == j || isANeighbor(j,randHost))			//if it is us or one of our neighbors....
						randHost = nodes*(rand()/(float)RAND_MAX);			//	then get get a different random host
				//	printf("it will be rewired to %d\n",randHost);
					delNeighbor(level[j/unit].id[j%unit].neighbor[k].ID,j);	//delete the host from the original neighbors array
					level[j/unit].id[j%unit].neighbor[k].ID=0-randHost;		//set the neighbor to the negative value of the random host, negative so that we know it has been rewired already
					level[j/unit].id[j%unit].neighbor[k].city = CITY_NUMBER;
					//add the host to the end of the neighbor hosts array
					level[randHost/unit].id[randHost%unit].neighbor[level[randHost/unit].id[randHost%unit].myNumNeighbors].ID=0-j; 
					level[randHost/unit].id[randHost%unit].neighbor[level[randHost/unit].id[randHost%unit].myNumNeighbors++].city = CITY_NUMBER;
					
			}
		}
	}
}//rewire

void global_init()
{
	int host, i, j, k;
	int foreigner;
	int host_city;
	int foreign_city;
	int code = -9;
	int num_crosstown_neighbors = n_cross_city_edges/2;
	
	for(k = 1; k <= num_Cities; k++)
	{
	
		if(CITY_NUMBER == k)
		{
			host = nodes*(rand()/(float)RAND_MAX);
			foreigner = nodes*(rand()/(float)RAND_MAX);
			host_city = CITY_NUMBER; //city sending below
			
			foreign_city = CITY_NUMBER;
			if(foreign_city == num_Cities)
				foreign_city = 0;
					
			for(i = 0; i < num_crosstown_neighbors; i++)
			{
				send(city_array[foreign_city],&code,sizeof(int),0);
				send(city_array[foreign_city], &host, sizeof(int), 0);
				send(city_array[foreign_city], &foreigner, sizeof(int), 0);
				send(city_array[foreign_city], &host_city, sizeof(int), 0);
		//		printf("sent host: %d to connect to %d from city %d\n",host,foreigner,foreign_city+1);
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors].ID = foreigner;
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors++].city = foreign_city+1;
				level[host/unit].id[host%unit].myNumForNeigh++;
				if(level[host/unit].id[host%unit].myNumForNeigh == 1)
					foreign_array[num_for_neigh++] = host; 
				host = nodes*(rand()/(float)RAND_MAX);
				foreigner = nodes*(rand()/(float)RAND_MAX);
			}
			
			foreign_city = CITY_NUMBER - 2;
			
			if(foreign_city == -1)
				foreign_city = num_Cities - 1;

		
			for(i = 0; i < num_crosstown_neighbors; i++)
			{
				send(city_array[foreign_city],&code,sizeof(int),0);
				send(city_array[foreign_city], &host, sizeof(int), 0);
				send(city_array[foreign_city], &foreigner, sizeof(int), 0);
				send(city_array[foreign_city], &host_city, sizeof(int), 0);
		//		printf("sent host: %d to connect to %d from city %d\n",host,foreigner,foreign_city+1);
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors].ID = foreigner;
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors++].city = foreign_city+1;
				level[host/unit].id[host%unit].myNumForNeigh++;
				if(level[host/unit].id[host%unit].myNumForNeigh == 1)
					foreign_array[num_for_neigh++] = host;
				host = nodes*(rand()/(float)RAND_MAX);
				foreigner = nodes*(rand()/(float)RAND_MAX);
			}
		}
		
	}
	
	/*GO_AHEAD = -10;
	if(CITY_NUMBER != num_Cities)
		send(city_array[CITY_NUMBER],&GO_AHEAD,sizeof(int),0);
	else
		send(city_array[0],&GO_AHEAD,sizeof(int),0);*/
}//global_init

void global_rewire()
{
	int i,j,k,l,m, temp, host;
	float randnumb;
	int randHost;
	int randCity;
	int code = -6;
	int foreign_city;

	for(j = 0; j < num_for_neigh; j++)
	{
		host = foreign_array[j];
		for(k = 0; k < level[host/unit].id[host%unit].myNumNeighbors; k++)
		{
			randnumb = (rand()/(float)RAND_MAX); //random number generation
					
			if((randnumb < swnP_city) && level[host/unit].id[host%unit].neighbor[k].ID > 0 && level[host/unit].id[host%unit].neighbor[k].city != CITY_NUMBER)
			{
				printf("node %d 's %d neighbor from city %d must be rewired\n", host, level[host/unit].id[host%unit].neighbor[k].ID, level[host/unit].id[host%unit].neighbor[k].city);
				randHost = nodes*(rand()/(float)RAND_MAX);				//get a random host from another city
				randCity = num_Cities*(rand()/(float)RAND_MAX);			//choose the random city
				randCity++;
				while(randCity == CITY_NUMBER)
				{
					randCity = num_Cities*(rand()/(float)RAND_MAX);
					randCity++;
				}
				while(isAForNeigh(j,randHost,randCity))			//if it is us or one of our neighbors....
					randHost = nodes*(rand()/(float)RAND_MAX);			//	then get get a different random host
				printf("it will be rewired to %d in city %d\n",randHost,randCity);
				//delForNeigh(host,level[j/unit].id[j%unit].neighbor[k].ID,level[j/unit].id[j%unit].neighbor[k].city);	//delete the host from the original neighbors array
				
				code = -6;
				foreign_city = level[host/unit].id[host%unit].neighbor[k].city-1;
				
				send(city_array[foreign_city],&code,sizeof(int),0);
				send(city_array[foreign_city], &host, sizeof(int), 0);
				send(city_array[foreign_city], &level[host/unit].id[host%unit].neighbor[k].ID, sizeof(int), 0);
				send(city_array[foreign_city], &CITY_NUMBER, sizeof(int), 0);
				
				level[host/unit].id[host%unit].neighbor[k].ID=0-randHost;		//set the neighbor to the negative value of the random host, negative so that we know it has been rewired already
				level[host/unit].id[host%unit].neighbor[k].city = randCity;
				//add the host to the end of the neighbor hosts array
				code = -9;
				foreign_city = randCity-1;
				send(city_array[foreign_city],&code,sizeof(int),0);
				send(city_array[foreign_city], &host, sizeof(int), 0);
				send(city_array[foreign_city], &randHost, sizeof(int), 0);
				send(city_array[foreign_city], &CITY_NUMBER, sizeof(int), 0);
			}
			
		}
	}
}//global_rewire

int isAForNeigh(int host, int foreigner, int city)
{
	int i;
	
	for(i = 0; i < level[host/unit].id[host%unit].myNumNeighbors; i++)
	{
		if(abs(level[host/unit].id[host%unit].neighbor[i].ID) == foreigner && level[host/unit].id[host%unit].neighbor[i].city == city)
			return 1;
	}
	
	return 0;
}//isAForNeigh

void delForNeigh(int host, int foreigner, int city)
{
	int i,j, k, found = 0;
	
	for(i = 0; i < level[host/unit].id[host%unit].myNumNeighbors; i++)
	{
		if((abs(level[host/unit].id[host%unit].neighbor[i].ID) == foreigner) && (level[host/unit].id[host%unit].neighbor[i].city == city) )
		{
			for(j = i; j < level[host/unit].id[host%unit].myNumNeighbors-1 ; j++)
				level[host/unit].id[host%unit].neighbor[j] = level[host/unit].id[host%unit].neighbor[j+1];
			found = 1;
			break;
		}
	}
	
		
	level[host/unit].id[host%unit].myNumNeighbors--;
	level[host/unit].id[host%unit].myNumForNeigh--;
	
	if(level[host/unit].id[host%unit].myNumForNeigh == 0)
	{
		for(k = 0; k < num_for_neigh-1; k++)
		{
			foreign_array[k] = foreign_array[k+1];
		}
		num_for_neigh--;
	}
	
	if(found == 0)
		printf("\n\n\n\n\n%d's NEIGHBOR %d from city %d NOT FOUND!!\n\n\n\n\n\n\n",host,foreigner,city);
	
}//dellForNeigh

int isANeighbor(int host, int randHost)
{
	int i;
	for(i = 0; i < level[host/unit].id[host%unit].myNumNeighbors; i++)
	{
		if(abs(level[host/unit].id[host%unit].neighbor[i].ID) == randHost)
			return 1;
	}
	
	return 0;
}//isANeighbor

void delNeighbor(int host, int delHost)
{
	int i,j;
	
	for(i = 0; i < level[host/unit].id[host%unit].myNumNeighbors; i++)
	{
		if(abs(level[host/unit].id[host%unit].neighbor[i].ID) == delHost)
		{
			for(j = i; j < level[host/unit].id[host%unit].myNumNeighbors-1 ; j++)
				level[host/unit].id[host%unit].neighbor[j] = level[host/unit].id[host%unit].neighbor[j+1];
			break;
		}
	}
		
	level[host/unit].id[host%unit].myNumNeighbors--;

}//delNeighbor

void print_people()
{
	int j,k,host;
	
	for(j = 0; j < num_for_neigh; j++)
	{
		host = foreign_array[j];
		printf("\nHost: %d neighbors: (%d of them)\n",host, level[host/unit].id[host%unit].myNumNeighbors);
		for(k = 0; k < level[host/unit].id[host%unit].myNumNeighbors; k++)
		{
			printf("\t%d city: %d\n",abs(level[host/unit].id[host%unit].neighbor[k].ID),level[host/unit].id[host%unit].neighbor[k].city);
		}
	}

}//print_people

void readvars()
{
	int i,j;
	char c;

	if (printinfo) 
			printf("in readvars\n");
	if ((input = fopen("fluvars.in", "r")) == NULL) 
	{
		printf("failed to open input file");
		exit(0);
	}
	fscanf(input, "%d %*s", &printinfo);
	//	printf("printinfo:%d\n",printinfo);
	fscanf(input, "%d %*s", &num_days);
	//	printf("num_days:%d\n",num_days);
	fscanf(input, "%d %*s", &num_runs);
	//	printf("num_runs:%d\n",num_runs);
	fscanf(input, "%d %*s", &num_reps);
	//	printf("num_reps:%d\n",num_reps);
	fscanf(input, "%d %*s", &num_net_types);
	//	printf("num_net_types:%d\n",num_net_types);
   
	if (num_net_types < 0 || num_net_types > 2) 
	{
		printf("only 1 or 2 net types allowed\n");
		exit(0);
	}
   
	for(i=0;i<num_net_types;i++)
	{
		fscanf(input, "%d", &net_types[i]);
	//	printf("net_types[%d]:%d\n",i,net_types[i]);
	}
   
	fscanf(input, "%d %*s", &num_p_values);
		//printf("num_p_values:%d\n",num_p_values);
   
	for(i=0;i < num_p_values;i++)
	{
		fscanf(input, "%f", &p_values[i]);
//		printf("p_values[%d]:%f\n",i,p_values[i]);
		if (p_values[i] < 0 || p_values[i] > 1.) 
		{
			printf("swnP = %f is not allowed\n",p_values[i]);
			exit(0);
		}
	}
	
	fscanf(input, "%s %*s", &c);
	vaccinate = atoi(&c);
//		printf("vaccinate:%d\n",vaccinate);
	fscanf(input, "%d %*s", &num_days_vacc);
		//printf("num_days_vacc:%d\n",num_days_vacc);
	
	for(i=0;i<num_days_vacc;i++)
	{
		fscanf(input, "%d", &vacc_days[i]);
		//printf("vacc_days[%d]:%d\n",i,vacc_days[i]);
		
		if(vacc_days[i] < 0 || vacc_days[i] > num_days)
		{
	//		printf("vacc_day = %f is not allowed\n", vacc_days[i]);
			exit(0);
		}
	}

	fscanf(input, "%d %*s", &num_vacc_strat);
	//	printf("num_vacc_strat:%d\n",num_vacc_strat);
	if (num_vacc_strat < 0) 
	{
		printf("cannot have a vacc_strat less than 0\n");
		exit(0);
	}
	for(i=0;i<num_vacc_strat;i++)
	{
		for(j=0;j<num_days_vacc;j++)
		{
			fscanf(input, "%d", &vacc_strats[i*num_days_vacc+j]);
	//		printf("vacc_strats[%d]:%d\n",i*num_days_vacc+j,vacc_strats[i*num_days_vacc+j]);
		}
	}
	fscanf(input, "%d %*s", &num_pc_vacc);
	//	printf("num_pc_vacc:%d\n",num_pc_vacc);
	
	for(i=0;i<num_pc_vacc;i++)
	{
		for(j=0;j<num_days_vacc;j++)
		{
			fscanf(input, "%f", &pc_vacc[i*num_days_vacc+j]);
	//		printf("pc_vacc[%d]:%f\n",i*num_days_vacc+j, pc_vacc[i*num_days_vacc+j]);
			if (pc_vacc[i*num_days_vacc+j] < 0. || pc_vacc[i*num_days_vacc+j] > 1.) 
			{
				printf("pc_vacc = %f is not allowed\n",pc_vacc[i]);
				exit(0);
			}
		}
	}
	fscanf(input, "%d %*s", &num_pop_size);
	//	printf("num_pop_size:%d\n",num_pop_size);
	for(i=0;i<num_pop_size;i++)
	{
		fscanf(input, "%d", &pop_sizes[i]);
	//	printf("pop_sizes[%d]:%d\n",i,pop_sizes[i]);
	}
	fscanf(input, "%d %*s", &num_NCR);
	//	printf("num_NCR:%d\n",num_NCR);
	for(i=0;i<num_NCR;i++)
	{
		fscanf(input, "%f", &NCRs[i]);
	//	printf("NCRs[%d]:%f\n",i,NCRs[i]);
		if (NCRs[i] < 0. || NCRs[i] > 1.) {
	//		printf("NCR = %f is not allowed\n",NCRs[i]);
			exit(0);
		}
	}
	fscanf(input, "%d %*s", &num_mutate);
	//	printf("num_mutate:%d\n",num_mutate);
	for(i=0;i<num_mutate;i++)
	{
		fscanf(input, "%f", &mutation_rates[i]);
	//	printf("mutation_rates[%d]:%f\n",i,mutation_rates[i]);
		if (mutation_rates[i] < 0. || mutation_rates[i] > 1.) 
		{
			printf("mutation_rate = %f is not allowed\n",mutation_rates[i]);
			exit(0);
		}
	}	
	fscanf(input, "%d %*s", &original_number_neighbors);
	//	printf("original_number_neighbors:%d\n",original_number_neighbors);
	fscanf(input, "%d %*s", &print_hosts);
	//	printf("print_hosts:%d\n",print_hosts);
	fscanf(input, "%d %*s", &tree);
	//	printf("tree:%d\n",tree);
	fscanf(input, "%d %*s", &print_strains);	
	//	printf("print_strains:%d\n",print_strains);
	fscanf(input, "%d %*s", &makepajek);
	//	printf("makepajek:%d\n",makepajek);
	if (makepajek && nodes > 1000)
		printf("WARNING:  Making Pajek file with %d hosts\n",nodes);
	fscanf(input, "%d %*s", &print_average);
	//	printf("print_average:%d\n",print_average);
	fscanf(input, "%d %*s", &makedegreedist);
	//	printf("makedegreedist:%d\n",makedegreedist);
	fscanf(input, "%d %*s", &numbits_Hsubtypes);
	//	printf("numbits_Hsubtypes:%d\n",numbits_Hsubtypes);
	fscanf(input, "%d %*s", &numbits_Nsubtypes);
	//	printf("numbits_Nsubtypes:%d\n",numbits_Nsubtypes);
	fscanf(input, "%d %*s", &numbits_strain);
	//	printf("numbits_strain:%d\n",numbits_strain);
		if (numbits_Hsubtypes > 32 || numbits_Nsubtypes > 32 || numbits_strain > 32) 
		{
			printf("Currently not supporting > 32-bit subtypes or strains\n");
			exit(0);
		}
	fscanf(input, "%d %*s", &latent_period);
	//	printf("latent_period:%d\n",latent_period);
	fscanf(input, "%d %*s", &num_days_infectious);
	//	printf("num_days_infectious:%d\n", num_days_infectious);
	fscanf(input, "%d %*s", &num_flu_vaccine);
	//	printf("num_flu_vaccine:%d\n",num_flu_vaccine);
	if (num_flu_vaccine > 3) 
	{
		printf("only allowing up to 3 strains in vaccine\n");
		exit(0);
	}
	for (i=0;i<num_flu_vaccine;i++) 
	{
		fscanf(input, "%d %d %d %*s", &vaccine_subtype_H[i],&vaccine_subtype_N[i],&vacstrain[i]);
	//	printf("vaccine_subtype_H[%d]:%d\tvaccine_subtype_N[%d]:%d\tvacstrain[%d]:%d\n",i,vaccine_subtype_H[i],i,vaccine_subtype_N[i],i,vacstrain[i]);
	}
	fscanf(input, "%d %*s", &num_infect);
	//	printf("num_infect:%d\n",num_infect);
	if (num_infect > 20) 
	{
		printf("only allowing up to 20 strains to infect\n");
		exit(0);
	}
	for (i=0;i<num_infect;i++) 
	{
		fscanf(input, "%d %d %d %*s", &infect_H[i],&infect_N[i],&originalstrain[i]);
	//	printf("infect_H[%d]:%d\tinfect_N[%d]:%d\toriginalstrain[%d]:%d\n",i,infect_H[i],i,infect_N[i],i,originalstrain[i]);
	}

	// print input file to flu daily output file header 
	
	fprintf(flu_sum_output,"NIM = %s\n",ver); // version number 
  	fprintf(flu_sum_output,"numruns = %d\n", num_runs);
  	fprintf(flu_sum_output,"num_days = %d\n", num_days);
  	fprintf(flu_sum_output,"original_number_neighbors = %d\n", original_number_neighbors);
	if (print_hosts)
		fprintf(flu_sum_output,"printing hosts\n");
	else
		fprintf(flu_sum_output,"not printing hosts\n");
	if (makepajek)
		fprintf(flu_sum_output,"making pajek file\n");
	else
		fprintf(flu_sum_output,"not making pajek file\n");
	if (print_average)
		fprintf(flu_sum_output,"making average file\n");
	else
		fprintf(flu_sum_output,"not making average file\n");
	if (makedegreedist)
		fprintf(flu_sum_output,"making degree distribution file\n");
	else
		fprintf(flu_sum_output,"not making degree distribution file\n");
		fprintf(flu_sum_output,"nettype ");
	for(i=0;i<num_net_types;i++)
	{
		fprintf(flu_sum_output,"%d ",net_types[i]);
	}
	fprintf(flu_sum_output,"\n");
  	fprintf(flu_sum_output,"swnP ");
	
	for(i=0;i<num_p_values;i++)
	{
		fprintf(flu_sum_output,"%f ",p_values[i]);
	}
	fprintf(flu_sum_output,"\n");
	fprintf(flu_sum_output,"vaccination strategy(ies) = ");
	
	for(i=0;i<num_vacc_strat;i++)
	{
		fprintf(flu_sum_output, "%d ", vacc_strats[i]);
	}
	fprintf(flu_sum_output,"|| where 0 = random; 1 = hubs; 2 = low clustering; 3 = high clustering; 4 = cross-cut");
	fprintf(flu_sum_output,"\ndays_of_vaccination");
	
	for(i=0;i<num_days_vacc;i++)
	{
		fprintf(flu_sum_output,"  %d",vacc_days[i]);
	}
	fprintf(flu_sum_output,"\npercent_vaccination");
	
	for(i=0;i<num_pc_vacc*num_days_vacc;i++)
	{
		fprintf(flu_sum_output, " %f", pc_vacc[i]);
	}
   
	fprintf(flu_sum_output,"\npopulation size");
	for(i=0;i<num_pop_size;i++)
	{
		fprintf(flu_sum_output, "%d ", pop_sizes[i]);
	}	
	fprintf(flu_sum_output,"\nNCR ");
	for(i=0;i<num_NCR;i++)
	{
		fprintf(flu_sum_output,"%f ", NCRs[i]);
	}
	fprintf(flu_sum_output,"\nmutation rates");
	for(i=0;i<num_mutate;i++)
	{
		fprintf(flu_sum_output,"%f ", mutation_rates[i]);
	}
	fprintf(flu_sum_output,"\n");
	fprintf(flu_sum_output,"numbits_Hsubtypes = %d\n", numbits_Hsubtypes);
  	fprintf(flu_sum_output,"numbits_Nsubtypes = %d\n", numbits_Nsubtypes);
	fprintf(flu_sum_output,"numbits_strain = %d\n", numbits_strain);
  	fprintf(flu_sum_output,"latent_period = %d\n", latent_period);
	fprintf(flu_sum_output,"num_days_infectious = %d\n", num_days_infectious);
	fprintf(flu_sum_output,"num strains in vaccine = %d\n",num_flu_vaccine);
	fprintf(flu_sum_output,"strains in vaccine: \n");
	for (i=0;i<num_flu_vaccine;i++) 
	{
		fprintf(flu_sum_output,"H = %d, N = %d, strain = %d\n", vaccine_subtype_H[i],vaccine_subtype_N[i],vacstrain[i]);
	}
	if (num_infect > 0) 
	{
		fprintf(flu_sum_output, "initially infect %d hosts with the following strain(s):\n", num_infect);
		for (i=0;i<num_infect;i++)
			fprintf(flu_sum_output, "H = %d, N = %d, strain = %d\n", infect_H[i],infect_N[i],originalstrain[i]);
	}
	else 
	{
  		fprintf(flu_sum_output, "no hosts are to be infected with flu\n");
	}
	fclose(input);
}//readvars


void makeGraphicFile()
{
	FILE *network;
	int q, r, i;
	GO_AHEAD = -10;
	int code = -10;
	
	if(CITY_NUMBER == 1)
	{
		network = fopen("network.dat","w");
		fprintf(network,"%d\n",num_Cities);
		fclose(network);
	}
	else
	{
		while(GO_AHEAD != 1)
		{}
	}
	
	network = fopen("network.dat","a");
	fprintf(network,"%d\n",nodes);
	for(q = 0; q < nodes; q++)
	{
		fprintf(network, "%d\t",level[q/unit].id[q%unit].myNumNeighbors);
		for(r = 0; r < level[q/unit].id[q%unit].myNumNeighbors; r++)
		{
			fprintf(network, "%d\t%d\t",abs(level[q/unit].id[q%unit].neighbor[r].ID),level[q/unit].id[q%unit].neighbor[r].city-1);
		}
		fprintf(network, "\n");
	}
	fclose(network);
	GO_AHEAD = -10;
	
	if(CITY_NUMBER == num_Cities)
		send(city_array[0],&code,sizeof(int),0);
	else
		send(city_array[CITY_NUMBER],&code,sizeof(int),0);
	
	while(GO_AHEAD != 1)
	{}
	
	if(CITY_NUMBER == 1)
	{	
		GO_AHEAD = -10;
		for(i = 1; i <= 3; i++)
			send(city_array[i],&GO_AHEAD,sizeof(int),0);
//		printf("SENT GO AHEADS TO ALL CITIES\n");
	}
	printf("ALL CITIES wrote to graphics file\n");
}//makeGraphicFile
