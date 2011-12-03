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



#define ver "2.0" /* Network Influenza Model (NIM) Version # */
#define first_host 0

int main( int argc, char **argv);
void runProg(char **argv);
void freeNetwork();
void acceptFromOtherCity();
void connectToCity();
void rewire();
int isANeighbor(int host, int randHost);
void delNeighbor(int host, int delHost);
void initInfect();
struct strain_list *new_strain_list();
struct flu *new_flu_strain();
void infect(int host, struct flu *f, int d);
void initialize(char **argv);
void vaccinate_hosts();
void readvars();
void vacc_high_clust();
void siftDownHighClust(int root, int bottom);
void vacc_clust();
void vacc_cross_cuts();
void vacc_hubs();
void cross_cut();
void degreedist();
void siftDownDegree(int root, int bottom);
void sort_cross_cuts();
void siftDownCrossCuts(int root, int bottom);
void only_one_across(int hostID);
void vacc_low_clust();
void siftDownLowClust(int root, int bottom);
void vaccinate_node(int IDnum);
void spread_flu_from_host(int host);
void spreadflu(void);
struct infectious *get_infectious_host_list(void);//changed struct to int
int is_infectious(int h);
float strains_matchN(int h);
float N_prop_match(int h);
struct infectious *new_infected_host(void);
struct flu *mutate(struct flu *f);
float percent_match(int a, int b, int nbits);
int hamming_distance(int a, int b, int nbits);
void make_pajek(char *argv[]);
void make_degree_distribution(char *argv[]);
void make_tree_file(char *argv[]);
void make_strains_file(char *argv[]);
void clustercoeff(void);
void recovery(void);
void statistics(void);
void summarystats(void);
void suminfected (int d);
void average(void);
void free_strain_list(void);
void print_tree();

struct people 
{
	int *neighbor;
	int myNumNeighbors;
	int cross;
	struct flu *flu_strains_infected; /* current strains in host */
	struct flu *flu_strains_recovered; /* past strains of flu */
	double clustco;
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

/* declarations for sockets */
int sockfd,newsockfd,clilen,childpid, n;
struct sockaddr_in cli_addr, serv_addr, serv_addr2;
char recvline[MAXLINE];

struct timeval tv;
fd_set readfds;


int main(int argc, char **argv)
{

	char s[120];
	
	pthread_t doWork, hear, speak;
	
	if (argc != 3) 
	{
		printf ("Try: argv[0] RandNumberSeed CityNumber\n");
		exit(0);
	}
	srand(25 * atoi(argv[1]));
	
	CITY_NUMBER = atoi(argv[2]);
	
	//srand((unsigned)time(NULL));
	strcpy(s,"sum_file");
	strcat(s,argv[1]);
	strcat(s," ");
	strcat(s,argv[2]);
	strcat(s,".csv");
	if ((flu_sum_output = fopen(s, "w")) == NULL) 
	{
		printf("failed to open flu summary output file: be sure it's closed in Excel\n");
		exit(0);
	}
	
	pthread_create(&doWork, NULL, runProg, argv);
	
	if(CITY_NUMBER == 1)
		pthread_create(&hear, NULL, acceptFromOtherCity, NULL);
	else if (CITY_NUMBER == 2)
		pthread_create(&speak, NULL, connectToCity, NULL);
	
	pthread_join(doWork, NULL);
	pthread_join(hear,NULL);
	pthread_join(speak,NULL);
}

void acceptFromOtherCity()
{
	int host,numbytes;

	tv.tv_sec = 2;
	tv.tv_usec = 500000;
	
	struct flu *firstFlu = (struct flu *) malloc(sizeof(struct flu));
	struct flu *f = (struct flu *) malloc(sizeof(struct flu));
	
	FD_ZERO(&readfds);
	
	/*
    * open a TCP socket - Internet stream socket
    */
   if ( (sockfd = socket(AF_INET, SOCK_STREAM,0)) < 0 ) 
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
	serv_addr.sin_port = htons(SERV_TCP_PORT);

	if ( bind(sockfd, (struct sockaddr *)&serv_addr,sizeof(serv_addr)) < 0) 
	{
		perror("server: can't bind local address");
		exit(1);
	}

	printf("waiting for City to connect.....\n");
	listen(sockfd,5);
  while(1) 
   {
		/*
		* wait for a connection from a client process
		* - concurrent server -
		*/
		clilen = sizeof(cli_addr);
		newsockfd = accept(sockfd, (struct sockaddr *) &cli_addr, &clilen); 
      
		if (newsockfd < 0 ) 
		{
			if (errno == EINTR)
				continue; /*try again*/
			perror("server: accept error");
			exit(1);
		}
		
		FD_SET(newsockfd, &readfds);
		
			printf("city connected\n");
		//	close(sockfd); 

			while(day < num_days)
			{
				recv(newsockfd, &host, sizeof(host), 0);
				if(host < nodes)
				printf("received %d host, gonna infect on %d\n",host,day);
				recv(newsockfd, &firstFlu->H, sizeof(int), 0);
				recv(newsockfd, &firstFlu->N, sizeof(int), 0);
				recv(newsockfd, &firstFlu->strain, sizeof(int), 0);
				recv(newsockfd, &firstFlu->old_strain, sizeof(int), 0);
				recv(newsockfd, &firstFlu->vaccine, sizeof(int), 0); 
					
				if(host < nodes)
					infect(host,firstFlu,day);
				host = nodes+1;
					
			}
		//}
		
	}//while 1
	
	close(newsockfd);
}//acceptfromothercity

void connectToCity()
{
	int hostToInfect = 670, i;
	
	struct flu *f = first_flu;

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
	

	printf("CONNECTED\n");
	
	while(day < num_days)// && hostToInfect < nodes)
	{
	/*	if( send(sockfd, &hostToInfect, sizeof(int), 0) == -1 )
		{
			printf("SEND ERROR!!\n");
			exit(1);
		}
		send(sockfd, &f->H, sizeof(int), 0);
		send(sockfd, &f->N, sizeof(int), 0);
		send(sockfd, &f->strain, sizeof(int), 0);
		send(sockfd, &f->old_strain, sizeof(int), 0);
		send(sockfd, &f->vaccine, sizeof(int), 0);
		printf("host %d sent\n",hostToInfect);
		hostToInfect++;
//		usleep(500000);
	*/
	}
	
	close(sockfd);
	
}//connectToCity


void runProg(char **argv)
{
	int a, b, c, f, l, k, j, g, y;
	char average_flu_day[120],s[120];
	int si;
	int numinfhosts, numrechosts, numinfections;
	struct flu *recovereds;
	int toSend = 512;
	
	readvars();
	fprintf(flu_sum_output,"\nNet_type,swnP,vacc_strats,per_vaccs,popsize,NCR,mut_rate,num_inf,tot_inf,duration,peak_inf,peak_day,half_pop_inf_day,# strains\n");
	if(!vaccinate) {
		printf("WARNING: Not vaccinating but summary output will report vacc_strat!\n");
		printf("\tAlthough slower, consider vaccinating but setting proportion \n\tof pop to vacc = 0.0%\n\n");
	}
	
	swnP = p_values[0];
	nodes = pop_sizes[0];
	printf("Population: %d\n",nodes);
	strategy = 0;
	vacc_strat = vacc_strats[strategy];
	percentage = 0;
	per_vaccinate = pc_vacc[percentage];
	NCR = NCRs[0];
	mut_rate = mutation_rates[0];
//	vaccinate = 1;
	
	for (y = 0; y < num_days; y++) 
	{
		totalinf[y]=0;
		totalrec[y]=0;
		totalinfections[y]=0;
		totalstrain[y]=0;
		numstrain[y]=0;
		averagestrain[y]=0;
		averagetotstrain[y]=0;
	}
	for (g = 1; g<=num_runs; g++)
	{
		initialize(argv);  /* set up model and do statistics */
		printf("done setting up network and rewiring\n");
		duration = 0;
		peak_num = 0;
		peak_day = 0;
		currentvacc = 0;
		vacc_day = vacc_days[currentvacc++];
		for (day = 0; day < num_days; day++) 
		{
			//usleep(20000);
			if(day <= 120)
				sleep(1);
			if(day == 10)
			{
					if(CITY_NUMBER == 2)
					{	
						if( send(sockfd, &toSend, sizeof(int), 0) == -1 )
						{
							printf("SEND ERROR!!\n");
							exit(1);
						}
					}
			}
			printf("%d\n",day);
			//printf("vaccinate %d\n",vaccinate);
			if(vaccinate==1 && day==vacc_day)
			{
//				printf("vaccinating on day %d\n",day);
				vaccinate_hosts(); /* vaccinates if its a vacc_day */
				vacc_day = vacc_days[currentvacc];
			}
			recovery(); /* recovery check of all hosts */
			spreadflu();
			statistics();

			numinfhosts = 0;
			numrechosts = 0;
			numinfections = 0;
			for (si=0; si < nodes; si++) 
			{
				if (level[si/unit].id[si%unit].flu_strains_infected) 
				{
					numinfhosts++;              /* count up infected hosts */
				}
				if (level[si/unit].id[si%unit].flu_strains_recovered) 
				{
						if (level[si/unit].id[si%unit].flu_strains_recovered->vaccine == 0) 
						{
							numrechosts++;              /* count up recovered hosts, but not vaccinated hosts */
						}
						recovereds = level[si/unit].id[si%unit].flu_strains_recovered;
						while(recovereds) 
						{
							if (recovereds->vaccine == 0) 
							{
								numinfections++;       /* counts all previous infections, but not vaccinations */
							}
							recovereds = recovereds->next;
						}
				}
			}
			totalinf[day]= numinfhosts;  /* add infections from current run to total */
			totalrec[day]= numrechosts;  /* add recovereds from current run to total */
			totalinfections[day]= numinfections;  /* add count all previous infections from current run to total */
			numstrain[day]= count_strains;  /* add count of strains/day from current run to total */
			totalstrain[day]=  tot_count_strains;  /* add cummulative count of strains from current run to total */
			}
			free_strain_list();
			freeNetwork();
		}
		average();
		summarystats();
		/******** Print Summary File ******/
		if (net_type == 0)
			fprintf(flu_sum_output,"SWN,");
		else
			fprintf(flu_sum_output,"SFN,");
		fprintf(flu_sum_output,"%4.2f,",swnP);
		for(j=0;j<num_days_vacc;j++)
		{
			fprintf(flu_sum_output," %d",vacc_strats[j]);
		}
		fprintf(flu_sum_output,",");
		for(j=0;j<num_days_vacc;j++)
		{
			fprintf(flu_sum_output," %4.2f",pc_vacc[j]);
		}
		fprintf(flu_sum_output,",%d,%4.2f,%4.4f,%4.2f,%4.2f,%d,%4.2f,%d,%d,%d\n",nodes,NCR,mut_rate,num_infected,averageinfections[num_days-1],duration,peak_num,peak_day,half_pop_inf,count_strains);
		fflush(flu_sum_output);
		if(print_average)
		{
			fclose(average_flu_days);

		}
		printf("NORMAL EXECUTION\n");	
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


void initInfect()
{
	struct flu *f;
	struct strain_list *s;
	int i, j, host;

   	// begin strain list, add all strains to initially infect with 
	first_strain = new_strain_list();
	first_strain->strain = originalstrain[0];
	first_strain->first_day = 0;
	first_strain->last_day_inf = 0;
	first_strain->old_strain = originalstrain[0];
	s = first_strain;
	
	for(j = 1; j < num_infect; j++) 
	{
		s->next = new_strain_list();
		s = s->next;
		s->strain = originalstrain[j];
		s->first_day = 0;
		s->last_day_inf = 0;
		s->old_strain = originalstrain[j];
	}
	
	
	for (i = 0; i < num_infect; i++) 
	{
		host = nodes*(rand()/(float)RAND_MAX); // get random host
		// get uninfected host and one that hasn't already been vaccinated 
		while (level[host/unit].id[host%unit].flu_strains_infected || level[host/unit].id[host%unit].flu_strains_recovered)
			host= nodes*(rand()/(float)RAND_MAX);
		f = new_flu_strain();
		f->H = infect_H[i];
		f->N = infect_N[i];
		f->strain = originalstrain[i];
		first_flu = f;
		
		
		infect(host,f,0);  //day = 0 
	}

}//initInfect

struct strain_list *new_strain_list()
{
	struct strain_list *l;

	l = (struct strain_list *) malloc(sizeof(struct strain_list));
	if (!l) 
	{
		printf("malloc of l in new_strain_list failed!\n");
		exit(0);
	}
	l->first_day = -1;
	l->last_day_inf = -1;
	l->strain = -1;
	l->old_strain = -1;
	l->instance = 0;
	l->hamming_distance = -1;
	l->next = NULL;
	return (l);
}//new_strain_list


struct flu *new_flu_strain()
{
	struct flu *f;

	f = (struct flu *) malloc(sizeof(struct flu));
	if (!f) 
	{
		printf("malloc of f in new_flu_strain failed!\n");
		exit(0);
	}
	f->H = -1; // make -1 = unassigned 
	f->N = -1; // make -1 = unassigned 
	f->strain = -1; // make -1 = unassigned 
	f->old_strain = 0;	
	f->vaccine = 0;
	f->infect_day = -1;
	f->next = NULL;
	return (f);
}// new_flu_strain

void infect(int host, struct flu *f, int d) // infects host "h" with flu "f" on day "d" 
{
	if (level[host/unit].id[host%unit].flu_strains_infected) // infect if not already infected 
	{ 
		printf("host isn't supposed to be infected\n");
		//exit(0);
	}
	else 
	{
//		printf("host:%d is getting infected\n",host);
		level[host/unit].id[host%unit].flu_strains_infected = f;
		level[host/unit].id[host%unit].flu_strains_infected->infect_day = d;
	}
}//infect


void vaccinate_hosts()
{
	struct flu *f_new, *f;
	int i, j, host;
	num_vaccinate = (int)(per_vaccinate*nodes);
	switch (vacc_strat) 
	{
		case 0:
	//			printf("vaccinating randomly\n");
			// randomly choose unvaccinated host to vaccinate 
				for (i = 0; i < num_vaccinate;i++) 
				{
					host = nodes*(rand()/(float)RAND_MAX); // get random host
					while(level[host/unit].id[host%unit].flu_strains_recovered || level[host/unit].id[host%unit].flu_strains_infected)  // already recovered or infected, get another 
					{
						host = nodes*(rand()/(float)RAND_MAX); // get random host
						printf("in here %d \n",host);
						
					}
					for (j = 0; j < num_flu_vaccine; j++)  // vaccinate with all strains 
					{
						f_new = new_flu_strain(); // flu strain from vaccine to give host 
						f_new->H = vaccine_subtype_H[j];
						f_new->N = vaccine_subtype_N[j];
						f_new->strain = vacstrain[j];
						f_new->vaccine = 1;
						f = level[host/unit].id[host%unit].flu_strains_recovered; // point to recovered list

						if (!f) // first vaccine flu strain in recovered list 
						{ 
							level[host/unit].id[host%unit].flu_strains_recovered = f_new;
						}
						else 
						{
							while (f->next) 
							{
								f = f->next;
							}
							f->next = f_new;
						}
					}	
				}	
	//			printf("done vaccinating randomly\n");
			break;
			
			case 1 :
			//	printf("vaccinating hubs\n");
				// randomly choose most connected inds to vaccinate 
				vacc_hubs();
			//	printf("done vaccinating hubs\n");
			break;
			
			case 2 :
			//	printf("vaccinating low clustering\n");
                 // choose nodes with lowest clustering coefficient
				vacc_low_clust();
			//	printf("done vaccinating low clust\n");
			break;
			case 3 :
			//	printf("high clust\n");
                 // choose nodes with highest clustering coefficient
				vacc_high_clust();
			//	printf("done with high clust\n");
			break;
			case 4 :
			//	printf("cross cuts\n");
				// choose nodes with cross-cut edges
				vacc_cross_cuts();
			//	printf("done with cross cuts\n");
			break;
			default :
				printf("problem in vaccinate_hosts()\n");
				exit(0);
		}
		strategy = strategy+1;
		percentage = percentage+1;
		if(vacc_days[currentvacc+1] && (strategy < num_vacc_strat))
		{
			currentvacc = currentvacc+1;
			vacc_day = vacc_days[currentvacc];
			vacc_strat = vacc_strats[strategy];
			per_vaccinate = pc_vacc[percentage];
		}
}//vaccinate_hosts


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

void vacc_clust()
{
	double clust_un_ob;
	int i, ra, shots_left, pot_shots = 0, pos = 0, s, counter;
	int h;
	
	clust_un_ob = level[pos/unit].id[pos%unit].clustco; /* array returns appropriately sorted */
	shots_left = num_vaccinate;
	while(shots_left > 0) 
	{
		for (i = pos; i < nodes; i++) 
		{
			if(clust_un_ob == level[IDs[i]/unit].id[IDs[i]%unit].clustco) 
			{
				pot_shots++;
			}
			else
				break;
		}
		if (pot_shots < shots_left)
		{
            //vaccinates all groups of the same degree until it comes
			//	to one with more indiviuals than vaccines left
			for (s = 0; s < pot_shots; s++) 
			{
				h = IDs[pos+s];
				if(level[h/unit].id[h%unit].flu_strains_recovered ||level[h/unit].id[h%unit].flu_strains_infected)
				{
					//do not vaccinate
				}
				else
				{
					vaccinate_node(IDs[pos+s]);
					shots_left--;
				//	printf("host:%d, vaccinated!!!!\n",h);
				}
			}
			pos += pot_shots;
			pot_shots = 0;
			clust_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].clustco;
		}
		else 
		{
            // vaccinates the last group randomly
			counter = pot_shots;
			
			for(i = pos; i < pos + pot_shots; i++)
			{
				h = IDs[i];
				if(level[h/unit].id[h%unit].flu_strains_recovered || level[h/unit].id[h%unit].flu_strains_infected)
				{
					IDs[i] = IDs[pos];
					IDs[pos] = h;
					pos++;
					counter--;
				}
			}
			
			while(counter>0 && shots_left > 0)
			{

					ra = ((int)((rand()/(float)RAND_MAX)*counter));
					counter--;
					h = IDs[pos+ra];
					vaccinate_node(h);
					shots_left--;
					IDs[pos+ra] = IDs[pos];
					IDs[pos] = h;
					pos++;
				
			}
			// if more vaccines than eligible nodes, go to next bin
			pot_shots = 0;
			clust_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].clustco;
		}
	}
}//vacc_clust


void vacc_hubs()
{
	int deg_un_ob, i, ra, shots_left, pot_shots = 0, pos = 0, s, counter;
	int h;
	
	degreedist();
	deg_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].myNumNeighbors;
	shots_left = num_vaccinate;
	pos = 0;
	while(shots_left > 0 && deg_un_ob > 0) 
	{
		for (i = pos; i < nodes; i++) 
		{
			if(deg_un_ob == level[IDs[i]/unit].id[IDs[i]%unit].myNumNeighbors) 
			{
				pot_shots++;
			}
			else
				break;
		}
		if (pot_shots < shots_left)
		{
            //vaccinates all groups of the same degree until it comes
			//	to one with more indiviuals than vaccines left
			for (s = 0; s < pot_shots; s++) 
			{
				h = IDs[pos+s];
				if(level[h/unit].id[h%unit].flu_strains_recovered ||level[h/unit].id[h%unit].flu_strains_infected)
				{
					//do not vaccinate
				}
				else
				{
					vaccinate_node(IDs[pos+s]);
			//		printf("(1)host:%d, vaccinated!!!!\n",IDs[pos+s]);
					shots_left--;
				}
			}
			pos += pot_shots;
			pot_shots = 0;
			deg_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].myNumNeighbors;
		}
		else 
		{
            // vaccinates the last group randomly
			counter = pot_shots;
			for(i = pos; i < pos + pot_shots; i++)
			{
				h = IDs[i];
				if(level[h/unit].id[h%unit].flu_strains_recovered || level[h/unit].id[h%unit].flu_strains_infected)
				{
					IDs[i] = IDs[pos];
					IDs[pos] = h;
					pos++;
					counter--;
				}
			}
			
			while(shots_left > 0 && counter > 0)
			{
				ra = ((int)((rand()/(float)RAND_MAX)*counter));
				counter--;
				h = IDs[pos+ra];
//				printf("swnP:%f\n",swnP);
//				for(i = 0; i < nodes; i ++)
//					printf("IDs[%d]:%d:myNumNeigh%d\n",i, IDs[i], level[IDs[i]/unit].id[IDs[i]%unit].myNumNeighbors);
//				printf("\n");
				vaccinate_node(h);
				shots_left--;
					
				IDs[pos+ra] = IDs[pos];
				IDs[pos] = h;
				pos++;
			}
			// if more vaccines than eligible nodes, go to next bin
			pot_shots = 0;
			deg_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].myNumNeighbors;
		}
	}
}//vacc_hubs

void vacc_cross_cuts()
{
   int deg_un_ob, i, ra, shots_left, pot_shots = 0, pos = 0, s, j, counter;
   int rn, lastr, firstr, temp, P1, P2;
   int h;
   cross_cut();
   for(j=0; j<nodes; j++) 
   {
      if(level[IDs[j]/unit].id[IDs[j]%unit].cross > neighbors) 
	  {
         only_one_across(IDs[j]);
      }
   }
   sort_cross_cuts();
   deg_un_ob = level[pos/unit].id[pos%unit].cross;
   shots_left = num_vaccinate;

	while(shots_left > 0 && deg_un_ob > 0) 
	{
		for (i = pos; i < nodes; i++) 
		{
			if(deg_un_ob == level[IDs[i]/unit].id[IDs[i]%unit].cross) 
			{
				pot_shots++;
			}
			else
				break;
		}
		if (pot_shots < shots_left)
		{
            //vaccinates all groups of the same degree until it comes
			//	to one with more indiviuals than vaccines left
			for (s = 0; s < pot_shots; s++) 
			{
				h = IDs[pos+s];
				if(level[h/unit].id[h%unit].flu_strains_recovered ||level[h/unit].id[h%unit].flu_strains_infected)
				{
					//do not vaccinate
				}
				else
				{
					vaccinate_node(IDs[pos+s]);
	//				printf("host:%d, vaccinated!!!!\n",IDs[pos+ra]);
					shots_left--;
				}
			}
			pos += pot_shots;
			pot_shots = 0;
			deg_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].cross;
		}
		else 
		{
            // vaccinates the last group randomly
			counter = pot_shots;
			
			for(i = pos; i < pos + pot_shots; i++)
			{
				h = IDs[i];
				if(level[h/unit].id[h%unit].flu_strains_recovered || level[h/unit].id[h%unit].flu_strains_infected)
				{
					IDs[i] = IDs[pos];
					IDs[pos] = h;
					pos++;
					counter--;
				}
			}
			
				
			while(counter > 0 && shots_left > 0)
			{
				ra = ((int)((rand()/(float)RAND_MAX)*counter));
				counter--;
				h = IDs[pos+ra];
				vaccinate_node(h);
//					printf("host:%d, vaccinated!!!!\n",IDs[pos+ra]);
				shots_left--;
				
				IDs[pos+ra] = IDs[pos];
				IDs[pos] = h;
				pos++;
			}
			// if more vaccines than eligible nodes, go to next bin
			pot_shots = 0;
			deg_un_ob = level[IDs[pos]/unit].id[IDs[pos]%unit].cross;
		}
	}

}//vacc_cross_cuts

void only_one_across(int hostID)
{
   int h;
   int nb;
   int i, neighID, diff, cross, big, small, far_neigh;

   h = hostID;
   cross = 0;
   nb = 0;
   while (nb < level[h/unit].id[h%unit].myNumNeighbors)
   { /* goes through all the neighbors to find neighbors far away */
		neighID = level[h/unit].id[h%unit].neighbor[nb];
		if(neighID > hostID) 
		{
			big = neighID;
			small = hostID;
			diff = neighID-hostID;
		}
		else 
		{
			big = hostID;
			small = neighID;
			diff = hostID-neighID;
		}
      
		if(diff>nodes/2) 
		{ /* hosts with very high ID numbers are close to those with very low ID numbers */
			diff = nodes-big+small;
		}
		
		if(diff>cross) 
		{ /* finds the farthest neighbor away */
			cross = diff;
			far_neigh = level[h/unit].id[h%unit].neighbor[nb];
		}
		nb++;
   }
   
   i = 0;
   
   while(i >= 0 && i < nodes) 
   { /* makes sure that when a cross-cut edge is
                                       found, only one of the nodes will be
                                       counted in the cross-cut array, so only
                                       one will be vaccinated */
		if(IDs[i] == far_neigh) 
		{
			level[IDs[i]/unit].id[IDs[i]%unit].cross = -1;
			i = -1;
		}
		else 
		{
			i++;
		}
   }
}//only one across


void cross_cut()
{
   int h, nb;
   int i, neighID, hostID, diff, cross, big, small;

   h = first_host;
   for (i = 0; i < nodes; i ++) 
   {
      hostID = i;
      cross = 0;
      nb = 0;
		while (nb < level[i/unit].id[i%unit].myNumNeighbors) 
		{ /* goes through all the neighbors to find neighbors far away*/
			neighID = level[i/unit].id[i%unit].neighbor[nb];
			if(neighID > hostID) 
			{
				big = neighID;
				small = hostID;
				diff = neighID-hostID;
			}
			else 
			{
				big = hostID;
				small = neighID;
				diff = hostID-neighID;
			}
         
			if(diff>nodes/2) 
			{
				diff = nodes - big + small;
			}
			if(diff > cross) 
			{
				cross = diff;
			}
			nb++;
		}
		
		level[i/unit].id[i%unit].cross = cross; /*  creates value of distance between farthest neighbor */
   }
   sort_cross_cuts();
}//cross_cut

void degreedist()
{
  int i, temp;

  for (i = (nodes / 2)-1; i >= 0; i--)
    siftDownDegree(i, nodes);

  for (i = nodes-1; i >= 1; i--)
  {
    temp = IDs[0];
    IDs[0] = IDs[i];
    IDs[i] = temp;
    siftDownDegree(0 , i-1);
  }
}//degreedist


void siftDownDegree(int root, int bottom)
{
  int done, maxChild, temp;

  done = 0;
  while ((root*2 <= bottom) && (!done))
  {
    if (root*2 == bottom)
      maxChild = root * 2;
    else if (level[IDs[root * 2]/unit].id[IDs[root * 2]%unit].myNumNeighbors < level[IDs[root * 2 + 1]/unit].id[IDs[root * 2 + 1]%unit].myNumNeighbors)
      maxChild = root * 2;
    else
      maxChild = root * 2 + 1;

    if (level[IDs[root]/unit].id[IDs[root]%unit].myNumNeighbors > level[IDs[maxChild]/unit].id[IDs[maxChild]%unit].myNumNeighbors)
    {
      temp = IDs[root];
      IDs[root] = IDs[maxChild];
      IDs[maxChild] = temp;
      root = maxChild;
    }
    else
      done = 1;
  }
}//siftDownDegree


void sort_cross_cuts()
{
  int i, temp;

  for (i = (nodes / 2)-1; i >= 0; i--)
    siftDownCrossCuts(i, nodes - 1);

  for (i = nodes-1; i >= 1; i--)
  {
    temp = IDs[0];
    IDs[0] = IDs[i];
    IDs[i] = temp;
    siftDownCrossCuts(0 , i-1);
  }
}//sort_cross_cuts


void siftDownCrossCuts(int root, int bottom)
{
  int done, maxCross, temp;

  done = 0;
  while ((root*2 <= bottom) && (!done))
  {
    if (root*2 == bottom)
      maxCross = root * 2;
    else if (level[IDs[root * 2]/unit].id[IDs[root * 2]%unit].cross < level[IDs[root * 2 + 1]/unit].id[IDs[root * 2 + 1]%unit].cross)
      maxCross = root * 2;
    else
      maxCross = root * 2 + 1;

    if (level[IDs[root]/unit].id[IDs[root]%unit].cross > level[IDs[maxCross]/unit].id[IDs[maxCross]%unit].cross)
    {
      temp = IDs[root];
      IDs[root] = IDs[maxCross];
      IDs[maxCross] = temp;
      root = maxCross;
    }
    else
      done = 1;
  }
}//siftDownCrossCuts



void vacc_low_clust()
{
  int i, temp;

  for (i = (nodes / 2)-1; i >= 0; i--)
    siftDownLowClust(i, nodes);

  for (i = nodes-1; i >= 1; i--)
  {
    temp = IDs[0];
    IDs[0] = IDs[i];
    IDs[i] = temp;
    siftDownLowClust(0 , i-1);
  }
  vacc_clust();
}//vacc_low_clust


void siftDownLowClust(int root, int bottom)
{
  int done, maxClust, temp;

  done = 0;
  while ((root*2 <= bottom) && (!done))
  {
    if (root*2 == bottom)
      maxClust = root * 2;
    else if (level[IDs[root * 2]/unit].id[IDs[root * 2]%unit].clustco > level[IDs[root * 2 + 1]/unit].id[IDs[root * 2 + 1]%unit].clustco)
      maxClust = root * 2;
    else
      maxClust = root * 2 + 1;

    if (level[IDs[root]/unit].id[IDs[root]%unit].clustco < level[IDs[maxClust]/unit].id[IDs[maxClust]%unit].clustco)
    {
      temp = IDs[root];
      IDs[root] = IDs[maxClust];
      IDs[maxClust] = temp;
      root = maxClust;
    }
    else
      done = 1;
  }
}//siftDownLowClust

void vacc_high_clust()
{
  int i, temp;

  for (i = (nodes / 2)-1; i >= 0; i--)
    siftDownHighClust(i, nodes);

  for (i = nodes-1; i >= 1; i--)
  {
    temp = IDs[0];
    IDs[0] = IDs[i];
    IDs[i] = temp;
    siftDownHighClust(0 , i-1);
  }
  vacc_clust(); 
}//vacc_high_clust


void siftDownHighClust(int root, int bottom)
{
  int done, maxClust, temp;

  done = 0;
  while ((root*2 <= bottom) && (!done))
  {
    if (root*2 == bottom)
      maxClust = root * 2;
    else if (level[IDs[root * 2]/unit].id[IDs[root * 2]%unit].clustco < level[IDs[root * 2 + 1]/unit].id[IDs[root * 2 + 1]%unit].clustco)
      maxClust = root * 2;
    else
      maxClust = root * 2 + 1;

    if (level[IDs[root]/unit].id[IDs[root]%unit].clustco > level[IDs[maxClust]/unit].id[IDs[maxClust]%unit].clustco)
    {
      temp = IDs[root];
      IDs[root] = IDs[maxClust];
      IDs[maxClust] = temp;
      root = maxClust;
    }
    else
      done = 1;
  }
}//siftDownHighClust

void vaccinate_node(int IDnum)
{
	int h;
	struct flu *f_new, *f;
	int j;

	h = IDnum;
	for (j = 0; j < num_flu_vaccine; j++) 
	{
		f_new = new_flu_strain();
		f_new->H = vaccine_subtype_H[j];
		f_new->N = vaccine_subtype_N[j];
		f_new->strain = vacstrain[j];
		f_new->vaccine = 1;
		f = level[h/unit].id[h%unit].flu_strains_recovered;
		if (!f) 
		{
			level[h/unit].id[h%unit].flu_strains_recovered = f_new;
			f_new->next = NULL;
		}
		else 
		{
			while (f->next) 
			{
				f = f->next;
			}
			f->next = f_new;
		}
	}
}//vaccinate_node

void spread_flu_from_host(int h)
	 /* spread flu from host *h: success based on NCR * [1-H_prop_match()]
	  *  where NCR = neighbor contact rate. If = 1. then all neighbors come into
	  *      contact with sick host and can get flu based only on H */
{
	int nb, count=0;
	struct flu *f, *recovereds, *recovstrains;
	float Hmatch, temp_match, strainmatch;
	int did_infection;
	int numspreadto = 0;
	
	float random = 0;

	did_infection = 0;
	if (!level[h/unit].id[h%unit].flu_strains_infected) 
	{
		printf("problem in spread_flu_from_host() function - host should be infected\n");
		exit(0);
	}
	/* first check if neighbor can get this flu strain */
	while (count < level[h/unit].id[h%unit].myNumNeighbors) 
	{ /* while h has a neighbor (nb) in list */
		nb = abs(level[h/unit].id[h%unit].neighbor[count]);
		
		f = new_flu_strain(); /* ASSUMING THIS FLU STRAIN WILL BE ADDED - MIGHT NOT BE */
		f->H = level[h/unit].id[h%unit].flu_strains_infected->H;
		f->N = level[h/unit].id[h%unit].flu_strains_infected->N;
		f->strain = level[h/unit].id[h%unit].flu_strains_infected->strain;
		Hmatch = 0.;
		
		if (!level[nb/unit].id[nb%unit].flu_strains_infected) 
		{ 
			/* if not already infected */
			/* Test if neighbors can get sick by matching H proteins */
			recovereds = level[nb/unit].id[nb%unit].flu_strains_recovered; /* beginning list of recovered strains */
			while (recovereds) 
			{
				/* check matching of H against each recovered strain */
				temp_match = percent_match(f->H, recovereds->H, numbits_Hsubtypes);
				if (temp_match > Hmatch)
					Hmatch = temp_match;
				recovereds = recovereds->next;
			}
			/****** ATTEMPT TO INFECT NEIGHBOR BASED ON H protein ********/
			/* this involves matching H proteins, probability of contact with
			 * *           neighbor */
			 random = (rand()/(float)RAND_MAX);
			if (Hmatch < 1 && random < NCR) 
			{
				/* if Hmatch < 1 (not perfect match) and randnum < NCR (neighbor contact rate)
				 * *                 then infect. */
				if (mut_rate > 0.) 
				{
					f->old_strain = f->strain;
					f = mutate(f);
				}
				infect(nb, f, day);
	//			printf("random number = %f\n",random);
				did_infection = 1;
			}
			if (Hmatch == 1) 
			{ 
				/*if Hs match perfectly, look at strains*/
				recovstrains = level[nb/unit].id[nb%unit].flu_strains_recovered; /* beginning list of recovered strains */
				strainmatch = 0;
				while (recovstrains) 
				{
					/* check matching of strain against each recovered strain */
					temp_match = percent_match(f->strain, recovstrains->strain, numbits_strain);
					if (temp_match > strainmatch) 
					{
						strainmatch = temp_match;
					}
					recovstrains = recovstrains->next;
				}
				/****** ATTEMPT TO INFECT NEIGHBOR BASED ON CROSS-IMMUNITY ********/
				/* this involves matching strains proteins, probability of contact with
				 * *                 neighbor */
				strainmatch = 1 - strainmatch;  /*  now strainmatch = probability of spreading
								    to neighbor (perfect match = 0.) */
				if ((rand()/(float)RAND_MAX) < strainmatch && (rand()/(float)RAND_MAX) < NCR) 
				{
					/* if randnum < NCR (neighbor contact rate) then infect. */
					if (mut_rate > 0.) 
					{
						f->old_strain = f->strain;
						f = mutate(f); 
					}
					infect(nb, f, day);
					did_infection = 1;
				}
			}
		}
		if (did_infection == 0) 
		{ /* didn't infect a host with f - need to free it */
			free(f);
		}
		count++;
	}
}//spread flu from host


void spreadflu(void)
{
	      /* 1. get list of infectious hosts.
	       * 2. randomly choose infectious host.
	       * 3. try to spread flu to each neighbor, then free infected host from list.
	       * 4. if not last host go to 2.
		   */

	struct infectious *infected_list, *first_infected_list;
	int countinfected, rand_infected;

	if (printinfo) printf("in spreadflu()\n");
	first_infected_list = get_infectious_host_list();
	countinfected = 0;
	infected_list = first_infected_list;
	while (infected_list) 
	{ /* first count number of infectious */
		countinfected++;
		infected_list = infected_list->next;
	}
	while (countinfected) 
	{ 
		/* do for all in infected_list */
		rand_infected = (int)(countinfected*(rand()/(float)RAND_MAX));
		infected_list = first_infected_list;
		while (rand_infected--) 
		{
			infected_list = infected_list->next;
		}
		spread_flu_from_host(infected_list->h);
		/* free memory for host pointed to by infected_list */
		if (infected_list->prev) 
		{
			if (infected_list->next) 
			{
				(infected_list->prev)->next = infected_list->next;
				(infected_list->next)->prev = infected_list->prev;
			}
			else 
			{
				(infected_list->prev)->next = NULL;
			}
		}
		else {
			if (infected_list->next) {
				first_infected_list = infected_list->next;
				first_infected_list->prev = NULL;
			}
		}
		free(infected_list);
		countinfected--;
	}
}


//struct infectious *get_infectious_host_list(void)
struct infectious *get_infectious_host_list(void)
	/* returns pointer to list of infectious hosts */
{
	int h;
	struct infectious *infected_list, *first_infected_host, *tempinfected;
	int infecteds, i;

	if (printinfo) printf("in get_infectious_host_list()\n");
	infecteds = 0;
	first_infected_host = NULL;
	for (i = 0; i < nodes; i++) 
	{
	
		if(level[i/unit].id[i%unit].flu_strains_infected != NULL)
		{
			if(day >= (level[i/unit].id[i%unit].flu_strains_infected->infect_day + latent_period) && 
				day < (level[i/unit].id[i%unit].flu_strains_infected->infect_day + latent_period + num_days_infectious))
			{
				if (is_infectious(i)) 
				{	
					if (infecteds == 0) 
					{
						first_infected_host = new_infected_host();
						infected_list = first_infected_host;
						infected_list->h = i; /* first host in infected_list points to h */
						infecteds=1; /* now there is at least one infected host */
					}
					else 
					{ /* add host to end of infected_list */
						tempinfected = new_infected_host();
						tempinfected->h = i;
						tempinfected->prev = infected_list; /* point to previous infected in list */
						infected_list->next = tempinfected; /* previous infected points to current infected */
						infected_list = tempinfected;
					}
				}
			}
		}
	}
	return (first_infected_host);
}


/* test if host h is infected:
*
*   h->flu_strains_infected
*
*   and
*
*      whether day (today) is in infectious period for host:
*
*         flu->infect_day + latent_period <= day <
*               flu->infect_day + latent_period + num_days_infectious
*
*                    and
*
*                          if N doesn't preclude production of virus particles N_prop_match() < 1
*                                   or 1 - N_prop_match() > 0.  A perfect match (1.0) would then yield
*                                               0 and this is used as the probability that spreading takes place. */

int is_infectious(int h)
{
//	printf("in is_infectious:%d\n", h);
	if ( 1 - N_prop_match(h) > 0.) 
	{
		return 1; /* host is infectious */
	}
	else if (1 - N_prop_match(h) == 0.) 
	{
		if ((rand()/(float)RAND_MAX) < 1 - strains_matchN(h)) 
		{ /* if strains are perfect match, not infectious */
			return 1;  /*host is infectious */
		}
		else 
		{
			return 0;  /* host is not infectious */
		}
	}
	else 
	{
		return 0;
	}
}

float strains_matchN(int h)
{
	struct flu *f;
	float best_match,temp_match;

	best_match = 0.;
	if (level[h/unit].id[h%unit].flu_strains_infected) 
	{
		if ((day - level[h/unit].id[h%unit].flu_strains_infected->infect_day) < (num_days_infectious + latent_period) &&
				day >= (level[h/unit].id[h%unit].flu_strains_infected->infect_day + latent_period)) 
		{
			f = level[h/unit].id[h%unit].flu_strains_recovered;
			while (f) 
			{
				temp_match = percent_match(level[h/unit].id[h%unit].flu_strains_infected->strain, f->strain, numbits_strain);
				if (temp_match > best_match) 
				{
					best_match = temp_match;
				}
				f = f->next;
			}
		}
	}
	return (best_match);
}

float N_prop_match(int h)
{
	   /* returns highest prop. match between host h's flu_strains_infected->N and
	    *       h's flu_strains_recovered-N */

	struct flu *f;
	float best_match,temp_match;
	
	best_match = 0.;
	if (level[h/unit].id[h%unit].flu_strains_infected) 
	{ 
		/* host is infected */
		/* test if current "day" is within infectious period */
		if ((day - level[h/unit].id[h%unit].flu_strains_infected->infect_day) < (num_days_infectious + latent_period) &&
				day >= (level[h/unit].id[h%unit].flu_strains_infected->infect_day + latent_period)) 
		{
			/* host is infectious and in infectious period. test if particle
			 * *             can be produced based on N proteins from previously seen strains. */
			f = level[h/unit].id[h%unit].flu_strains_recovered;
			while (f) 
			{
				temp_match = percent_match(level[h/unit].id[h%unit].flu_strains_infected->N,f->N,numbits_Nsubtypes);
				if (temp_match > best_match)
					best_match = temp_match;
				f = f->next;
			}
		}
	}
	return (best_match);
}

struct infectious *new_infected_host(void)
{
	struct infectious *s;

	s = (struct infectious *) malloc(sizeof(struct infectious));
	if (!s) 
	{
		printf("malloc of s in new_infectious failed!\n");
		exit(0);
	}
	s->h = -1;
	s->prev = NULL;
	s->next = NULL;
	return (s);
}
struct flu *mutate(struct flu *f)
{
	long mutationsite, site;
	/*** equal mutation probability per 'nucleotide' ***/
	
	for (site = 0; site < numbits_strain; site++)
	{
		if ((rand()/(float)RAND_MAX) < mut_rate) {
			mutationsite = 1 << site;
			f->strain = (f->strain & mutationsite) ?
				f->strain & (~mutationsite) : f->strain | mutationsite;
		}
	}
	return (f);
}

float percent_match(int a, int b, int nbits)
{
	int concord, i, mask;
	
	concord = 0;
	for (i = 0; i < nbits; i++) {
		mask = 1 << i;
		if ((a & mask) == (b & mask))
			concord++;
	}
	return ((float)concord/nbits);
}

int hamming_distance(int a, int b, int nbits)
{
	int hamming_dist, concord, i, mask;

	concord = 0;
	for (i = 0; i < nbits; i++) 
	{
		mask = 1 << i;
		if ((a & mask) == (b & mask))
			concord++;
	}
	hamming_dist = nbits - concord;
	return (hamming_dist);
}

void make_pajek(char *argv[])
{
	int nb;
	int i, j;
	char s[30];
	FILE *pajek;

	strcpy(s,"flu");
	strcat(s,argv[1]);
	strcat(s,".net");
	if ((pajek = fopen(s, "w")) == NULL) 
	{
		printf("failed to open pajek output file\n");
		exit(0);
	}
	fprintf(pajek,"*Vertices %d\n",nodes);
	for (i = 0; i < nodes; i++) 
	{
		fprintf(pajek,"%d \"%d\"\n",i+1,i);
	}
	fprintf(pajek,"*Edges\n");
	for (i = 0; i < nodes; i++) 
	{
		for (j = 0; j < level[i/unit].id[i%unit].myNumNeighbors; j++)
		{
			fprintf(pajek,"%d %d\n",i,level[i/unit].id[i%unit].neighbor[j]);
		}
	}
	fclose(pajek);
}

void make_degree_distribution(char *argv[])
{
	int h, nb;
	int i, numneigh;
	double lognumneigh, lognode, lnode, lnumneigh;
	char s[30];
	FILE *degree;

	/* srand(25*atoi(argv[1])); WHAT IN THE WORLD WAS THIS FOR??? COMMENTED 7-14-04 */
	strcpy(s,"degree");
	strcat(s,argv[1]);
	strcat(s,".csv");
	if ((degree = fopen(s, "w")) == NULL) {
		printf("failed to open degree distribution output file\n");

		exit(0);
	}
	fprintf(degree,"node ID, degree, log ID, log degree\n");
	h = first_host;
	for (i = 0; i < nodes; i++) {
		numneigh = level[i/unit].id[i%unit].myNumNeighbors;
		lnumneigh = numneigh;
		lnode = i;
		lognode = log10(lnode);
		lognumneigh = log10(lnumneigh);
		fprintf(degree,"%i, %i, %f, %f\n", i+1, numneigh, lognode, lognumneigh);
	}	
}

void make_tree_file(char *argv[])
{
	char tree_file[80], s[20];

	countrun++;
	strcpy(tree_file,"tree file");
	strcat(tree_file,argv[1]);
	strcat(tree_file,"_reps");
	sprintf(s,"%d",countrun);
	strcat(tree_file,s);
	strcat(tree_file,".csv");
	if ((tree_files = fopen(tree_file, "w")) == NULL) {
		printf("failed to open average_flu_day (%s) output file\n",tree_files);
		exit(0);
	}
	fprintf(tree_files,"newbit,oldbit,Strain,Oldstrain,hamdist,1st day,tot day,last day\n");
}

void make_strains_file(char *argv[])
{
	char strain_file[80], s[20];

	if (tree == 0)
		countrun++;
	strcpy(strain_file,"strains");
	strcat(strain_file,argv[1]);
	strcat(strain_file,"_reps");
	sprintf(s,"%d",countrun);
	strcat(strain_file,s);
	strcat(strain_file,"_swnP");
	sprintf(s,"%f",swnP);
	strcat(strain_file,s);
	strcat(strain_file,".csv");
	if ((strain_files = fopen(strain_file, "w")) == NULL) 
	{
		printf("failed to open average_flu_day (%s) output file\n",strain_files);
			exit(0);
	}
	fprintf(strain_files,"day,# strains,total #\n");
}

void clustercoeff(void)
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

void recovery(void)
{
	   struct flu *f;
	   int i;

	   for (i = 0; i < nodes; i++) 
	   {
		   if (level[i/unit].id[i%unit].flu_strains_infected &&    /* if infected and */
				   day == level[i/unit].id[i%unit].flu_strains_infected->infect_day + num_days_infectious + latent_period) 
		   {
			   /* it's the day of recovery, this individual recovers! */

			   f = level[i/unit].id[i%unit].flu_strains_recovered; /* point to recovered list */
			   if (!f) 
			   { 
				   /* no recovered strains */
				   level[i/unit].id[i%unit].flu_strains_recovered = level[i/unit].id[i%unit].flu_strains_infected;
				   level[i/unit].id[i%unit].flu_strains_recovered->next = NULL;
				   level[i/unit].id[i%unit].flu_strains_infected = NULL; /* host no longer infected */
//			   	   printf("host: %d recovered\n", i);
			   }
			   else 
			   { /* zoom to last recovered strain */
				   while (f->next) 
				   {
					   f = f->next;
				   }
				   f->next = level[i/unit].id[i%unit].flu_strains_infected; /* point to flu to be recovered */
				   level[i/unit].id[i%unit].flu_strains_infected = NULL; /* host no longer infected */
//			   	   printf("host: %d recovered\n", i);
			   }
		   }
	   }
}

void statistics(void)
{
	int i, count_s;
	struct strain_list *s, *l, *temp_s;

	if (printinfo) printf("in statistics()\n");
	
	for (i = 0; i < nodes; i++) 
	{
		count_s = 0;
		if (level[i/unit].id[i%unit].flu_strains_infected && level[i/unit].id[i%unit].flu_strains_infected->infect_day < day + latent_period) 
		{
			s = first_strain;
			while(s) 
			{
				
				if(s->strain == level[i/unit].id[i%unit].flu_strains_infected->strain) 
				{  
					/* if strain is already in list, check to see if it died out, add to list if it has */
					count_s = 1;
					if (s->last_day_inf <= 1 || day - s->last_day_inf < 5)  /* if this is the first day of infection or the strain is still alive, increment last day */
						s->last_day_inf = day;
					if (day - s->last_day_inf >= 5 && s->instance < 2) 
					{   
						/* if strain has died out add it to the list, and increment instance */
						l = first_strain;
						while(l) 
						{
							/* go to end of list and add strain */
							if(l->next == NULL)
								temp_s = l;
							l=l->next;
						}
						l = new_strain_list();
						temp_s->next = l;
						l->strain = level[i/unit].id[i%unit].flu_strains_infected->strain;
						l->old_strain = level[i/unit].id[i%unit].flu_strains_infected->old_strain;
						l->hamming_distance = hamming_distance(l->strain,l->old_strain,numbits_strain);
						l->first_day = day;
						l->last_day_inf = 0;
						l->instance = 1;
						l->next = NULL;
						s->instance = 2;
					}
				}
				if(s->next == NULL) 
				{
					temp_s = s;
				}
				s = s->next;
			}
			
			if(count_s == 0) 
			{ 
				/* if strain isn't in list yet, add it */
				s = new_strain_list();
				temp_s->next = s;
				s->strain = level[i/unit].id[i%unit].flu_strains_infected->strain;
				s->old_strain = level[i/unit].id[i%unit].flu_strains_infected->old_strain;
				s->hamming_distance = hamming_distance(s->strain,s->old_strain,numbits_strain);
				s->first_day = day;
				s->last_day_inf = 0;
				s->instance = 1;
				s->next = NULL;
			}
		}
		if (print_hosts) 
		{ 
			/* print all infected hosts to a file */
			if (level[i/unit].id[i%unit].flu_strains_infected) 
			{ 
				/* infected host */
				fprintf(host_output,"%d,%d,%d,%d,%d,%d\n",
						day,i,level[i/unit].id[i%unit].flu_strains_infected->infect_day,
						level[i/unit].id[i%unit].flu_strains_infected->H,level[i/unit].id[i%unit].flu_strains_infected->N, level[i/unit].id[i%unit].flu_strains_infected->strain);
			}
		}
	}

	count_strains = 0;
	tot_count_strains = 0;
	temp_s = first_strain;
	while(temp_s) 
	{
		if(day >= temp_s->first_day && day <= temp_s->last_day_inf)
			count_strains++;
		tot_count_strains++;
		temp_s = temp_s->next;
	}
	if(print_strains)
		fprintf(strain_files,"%d,%d,%d\n",day,count_strains,tot_count_strains);
}

void summarystats(void)
{
	int i;

	for(i = 0; i < num_days; i++) 
	{
		if(averageinf[i] > peak_num) 
		{
			peak_num = averageinf[i];       /* records peak of infections for stats file */
			peak_day = i;
		}
		if(averageinf[i] == 0 && averageinf[i-1] != 0) 
		{
			duration = i;                    /* records duration for stats file */
		}
		if(averageinf[num_days - 1] > 0) 
		{
			duration = num_days -1;                   /* if duration exceeds run time set duration to last day */
		}
		if(averagerec[i-1]<nodes/2 && averagerec[i]>=nodes/2)
		{
			half_pop_inf = i;
		}
	}
	num_infected = averagerec[num_days - 1];  /* amount of people that got sick is taken from the count of recovereds on the last day */
	count_strains = averagetotstrain[num_days - 1];  /* strains reported is cummulative strains from last day */
}

void suminfected (int d)
{
	/* adds the total number of nodes infected each run to the previous total*/
	   int i;
	   int numinfhosts, numrechosts, numinfections;
	   struct flu *recovereds;

	   numinfhosts = 0;
	   numrechosts = 0;
	   numinfections = 0;
	   for (i=0; i < nodes; i++) 
	   {
		   if (level[i/unit].id[i%unit].flu_strains_infected) 
		   {
			   numinfhosts++;              /* count up infected hosts */
		   }
		   if (level[i/unit].id[i%unit].flu_strains_recovered) 
		   {
			   if (level[i/unit].id[i%unit].flu_strains_recovered->vaccine == 0) 
			   {
				   numrechosts++;              /* count up recovered hosts, but not vaccinated hosts */
			   }
			   recovereds = level[i/unit].id[i%unit].flu_strains_recovered;
			   while(recovereds) 
			   {
				   if (recovereds->vaccine == 0) 
				   {
					   numinfections++;       /* counts all previous infections, but not vaccinations */
				   }
				   recovereds = recovereds->next;
			   }
		   }
	   }
	   totalinf[d]= numinfhosts;  /* add infections from current run to total */
	   totalrec[d]= numrechosts;  /* add recovereds from current run to total */
	   totalinfections[d]= numinfections;  /* add count all previous infections from current run to total */
	   numstrain[d]= count_strains;  /* add count of strains/day from current run to total */
	   totalstrain[d]=  tot_count_strains;  /* add cummulative count of strains from current run to total */
}

void average(void)
{
	/* averages the number of infected individuals for each day*/
	int y;

	if(print_average)
		fprintf(average_flu_days,"Day, aveNinfected, aveNrecovered, aveinfections, # strains, totstrains\n");
	for (y = 0; y < num_days; y++) 
	{
		averageinf[y] = totalinf[y]/num_runs;     /* average infected */
		averagerec[y] = totalrec[y]/num_runs;     /* average recovered */
		averageinfections[y] = totalinfections[y]/num_runs;  /* average total infections */
		averagestrain[y] = numstrain[y]/num_runs;  /* average strains present per day */
		averagetotstrain[y] = totalstrain[y]/num_runs;  /* average cummulative strains */
		if(print_average)
			fprintf(average_flu_days,"%i,%f,%f,%f,%f,%f\n",y,averageinf[y],averagerec[y],averageinfections[y],averagestrain[y],averagetotstrain[y]);
	}
}

void free_strain_list(void)
{
	struct strain_list *s, *temp_s;

	s = first_strain;
	while(s) 
	{
		temp_s = s;
		s = s->next;
		free(temp_s);
	}
}

void print_tree(void)
{
	struct strain_list *s;
	int new_bit[33], old_bit[33], numbits, i, new_strain, old_strain;

	s = first_strain;
	while(s) {
		new_strain = s->strain;
		old_strain = s->old_strain;
		for(i = 0; i < numbits_strain; i++) {
			numbits = numbits_strain - i;
			new_bit[i] = ((int)new_strain/(pow(2,numbits-1)));
			old_bit[i] = ((int)old_strain/(pow(2,numbits-1)));
			if(new_bit[i] == 1)
				new_strain = new_strain - (pow(2,numbits-1));
			if(old_bit[i] == 1)
				old_strain = old_strain - (pow(2,numbits-1));
			fprintf(tree_files,"%d",new_bit[i]);
		}
		fprintf(tree_files,",");
		for(i = 0; i < numbits_strain; i++) {
			fprintf(tree_files,"%d",old_bit[i]);
		}
		fprintf(tree_files,",");
		fprintf(tree_files,"%d,%d,%d,%d,%d,%d\n",s->strain,s->old_strain,s->hamming_distance,s->first_day,s->last_day_inf - s->first_day + 1,s->last_day_inf + 1);
		s=s->next;
	}
}
