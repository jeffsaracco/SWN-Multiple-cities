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
#include <stdint.h>
#include "inet.h"



#define ver "3.0" /* Network Influenza Model (NIM) Version # */
#define first_host 0

struct flu;
int main( int argc, char **argv);
void runProg(char **argv);
void freeNetwork();
void read_cities(const char *path);
int connect_to_city(int city);
void connect_cities();
void send_ints(int city, const int *values, int n);
void send_code(int city, int code);
void send_cmd(int city, int code, int host, int foreigner, int host_city);
int recv_int(int sock, int *value);
void *get_commands(void *arg);
void barrier();
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
void send_infection(int city, int host, struct flu *f);
void apply_imported_infections();
void report_day();
void write_epidemic_meta();
int command_arg_count(int code);
void initInfect();
struct strain_list *new_strain_list();
struct flu *new_flu_strain();
void infect(int host, struct flu *f, int d);
void vaccinate_hosts();
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
int try_infect(int nb, struct flu *f);
void spreadflu(void);
struct infectious *get_infectious_host_list(void);
int is_infectious(int h);
float strains_matchN(int h);
float N_prop_match(int h);
struct infectious *new_infected_host(void);
struct flu *mutate(struct flu *f);
float percent_match(int a, int b, int nbits);
int hamming_distance(int a, int b, int nbits);
void clustercoeff(void);
void recovery(void);
void statistics(void);
void summarystats(void);
void average(void);
void free_strain_list(void);

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
int printinfo, vaccinate;
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

/* Commands sent between cities; command_arg_count() says how many ints follow each.
   Add/delete: the sender's host, the receiver's host, the sender's city number.
   Infect: day, the receiver's host, and the strain's H, N and strain bits.
   Day stats (to city 1): day, infected, recovered, vaccinated, then the number
   of infections imported today from each city. */
#define CMD_DEL_FOREIGN_NEIGHBOR  -6
#define CMD_ADD_FOREIGN_NEIGHBOR  -9
#define CMD_BARRIER              -10
#define CMD_INFECT               -13
#define CMD_DAY_STATS            -14
#define DAY_STATS_FIELDS           4 /* day, infected, recovered, vaccinated */

#define MAX_CITIES 64

/* declarations for sockets */
int num_Cities;
char city_host[MAX_CITIES][256];	/* address of each city, indexed by city-1 */
int city_port[MAX_CITIES];
char city_name[MAX_CITIES][64];		/* where the visualization puts each city */
float city_lat[MAX_CITIES], city_lon[MAX_CITIES];
int city_array[MAX_CITIES];			/* socket to each other city, indexed by city-1 */
int n_cross_city_edges = 10;
int num_for_neigh;
int *foreign_array;
float swnP_city = 0.5;
int seed_city = 1;					/* city that gets the initial infections (0 = every city) */
FILE *epidemic_output;				/* per-day counts for every city, written by city 1 */
int imports_today[MAX_CITIES];		/* infections that took here today, by source city */

/* net_lock guards the network and the fields below, which the get_commands
   threads update as commands arrive from other cities. */
pthread_mutex_t net_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t net_cond = PTHREAD_COND_INITIALIZER;
int barriers_from[MAX_CITIES];		/* barriers received from each city */
int barriers_sent = 0;
struct imported_infection
{
	int day, host, H, N, strain;
	struct imported_infection *next;
};
struct imported_infection *imported_head[MAX_CITIES], *imported_tail[MAX_CITIES]; /* by source city */
int city_stats[MAX_CITIES][DAY_STATS_FIELDS + MAX_CITIES];	/* city 1 only: latest day stats from each city */


int main(int argc, char **argv)
{
	char s[120];
	int i;
	pthread_t commands[MAX_CITIES];

	if (argc != 3 && argc != 4)
	{
		printf ("Try: %s RandNumberSeed CityNumber [CitiesFile]\n",argv[0]);
		exit(0);
	}
	setvbuf(stdout, NULL, _IOLBF, 0);
	CITY_NUMBER = atoi(argv[2]);
	srand(25 * atoi(argv[1]) + CITY_NUMBER - 1); // city 1 matches the original seeding
	read_cities(argc == 4 ? argv[3] : "cities.txt");
	if (CITY_NUMBER < 1 || CITY_NUMBER > num_Cities)
	{
		printf("CityNumber must be between 1 and %d\n", num_Cities);
		exit(1);
	}

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

	connect_cities();

	for(i = 1; i <= num_Cities; i++) //create commands listening threads
	{
		if( i != CITY_NUMBER )
			pthread_create(&commands[i-1], NULL, get_commands, (void *)(intptr_t) i);
	}

	runProg(argv); // exits when the run is done
	return 0;
}

/* Reads one "host port" line per city; line N is city N. */
void read_cities(const char *path)
{
	FILE *f;
	char line[512];
	int fields, name_start, i;
	int placed[MAX_CITIES] = { 0 };

	if ((f = fopen(path, "r")) == NULL)
	{
		perror(path);
		exit(1);
	}
	num_Cities = 0;
	while (fgets(line, sizeof(line), f))
	{
		if (line[0] == '#' || line[0] == '\n')
			continue;
		if (num_Cities == MAX_CITIES)
		{
			printf("too many cities in %s (max %d)\n", path, MAX_CITIES);
			exit(1);
		}
		fields = sscanf(line, "%255s %d %f %f %n", city_host[num_Cities], &city_port[num_Cities],
			&city_lat[num_Cities], &city_lon[num_Cities], &name_start);
		if (fields < 2)
		{
			printf("bad line in %s: %s", path, line);
			exit(1);
		}
		placed[num_Cities] = (fields == 4);
		city_name[num_Cities][0] = 0;
		if (placed[num_Cities])
		{
			line[strcspn(line, "\r\n")] = 0;
			snprintf(city_name[num_Cities], sizeof(city_name[0]), "%s", line + name_start);
		}
		if (city_name[num_Cities][0] == 0)
			snprintf(city_name[num_Cities], sizeof(city_name[0]), "City %d", num_Cities + 1);
		num_Cities++;
	}
	fclose(f);
	if (num_Cities < 2)
	{
		printf("%s must list at least 2 cities\n", path);
		exit(1);
	}
	for (i = 0; i < num_Cities; i++)
		if (!placed[i]) // no location given: spread cities around the equator
		{
			city_lat[i] = 0;
			city_lon[i] = -180 + 360.0 * i / num_Cities;
		}
}//read_cities

/* Opens a TCP connection to a city, retrying while it starts up. */
int connect_to_city(int city)
{
	struct addrinfo hints, *addr;
	char port[16];
	int sock, attempt;

	memset(&hints, 0, sizeof(hints));
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	snprintf(port, sizeof(port), "%d", city_port[city-1]);

	printf("Trying to connect to city %d...\n",city);
	for (attempt = 0; attempt < 60; attempt++, sleep(1))
	{
		if (getaddrinfo(city_host[city-1], port, &hints, &addr) != 0)
			continue; // name not resolvable yet (e.g. container still starting)
		sock = socket(AF_INET, SOCK_STREAM, 0);
		if (sock < 0)
		{
			perror("client: can't open stream socket");
			exit(1);
		}
		if (connect(sock, addr->ai_addr, addr->ai_addrlen) == 0)
		{
			freeaddrinfo(addr);
			return sock;
		}
		close(sock);
		freeaddrinfo(addr);
	}
	printf("client: can't connect to city %d at %s:%s\n", city, city_host[city-1], port);
	exit(2);
}//connect_to_city

/* Connects every pair of cities once: each city dials the lower-numbered
   cities and accepts connections from the higher-numbered ones.  The dialing
   city sends its number first so the acceptor knows who it is. */
void connect_cities()
{
	struct sockaddr_in serv_addr;
	int listener, sock, city, i, on = 1;

	if ( (listener = socket(AF_INET, SOCK_STREAM,0)) < 0 )
	{
		perror("server: can't open stream socket");
		exit(1);
	}
	setsockopt(listener, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));
	memset((char *) &serv_addr, 0, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
	serv_addr.sin_port = htons(city_port[CITY_NUMBER-1]);
	if ( bind(listener, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
	{
		perror("server: can't bind local address");
		exit(1);
	}
	listen(listener, MAX_CITIES);

	for (city = 1; city < CITY_NUMBER; city++)
	{
		city_array[city-1] = connect_to_city(city);
		send_ints(city, &CITY_NUMBER, 1);
		printf("CONNECTED to city %d\n",city);
	}

	for (i = CITY_NUMBER + 1; i <= num_Cities; i++)
	{
		printf("waiting for a city to connect.....\n");
		sock = accept(listener, NULL, NULL);
		if (sock < 0 )
		{
			perror("server: accept error");
			exit(1);
		}
		if (!recv_int(sock, &city) || city <= CITY_NUMBER || city > num_Cities)
		{
			printf("server: bad hello from a connecting city\n");
			exit(1);
		}
		city_array[city-1] = sock;
		printf("city %d connected\n",city);
	}
	close(listener);
}//connect_cities

void send_ints(int city, const int *values, int n)
{
	if (send(city_array[city-1], values, n * sizeof(int), 0) != n * (int) sizeof(int))
	{
		perror("SEND ERROR");
		exit(1);
	}
}

void send_code(int city, int code)
{
	send_ints(city, &code, 1);
}

void send_cmd(int city, int code, int host, int foreigner, int host_city)
{
	int msg[4] = { code, host, foreigner, host_city };
	send_ints(city, msg, 4);
}

int recv_int(int sock, int *value)
{
	return recv(sock, value, sizeof(int), MSG_WAITALL) == sizeof(int);
}

int command_arg_count(int code)
{
	switch (code)
	{
		case CMD_ADD_FOREIGN_NEIGHBOR:
		case CMD_DEL_FOREIGN_NEIGHBOR:
			return 3;
		case CMD_INFECT:
			return 5;
		case CMD_DAY_STATS:
			return DAY_STATS_FIELDS + num_Cities;
		default:
			return 0;
	}
}

/* One thread per other city: applies the commands that city sends us. */
void *get_commands(void *arg)
{
	int city = (int)(intptr_t) arg;
	int sock = city_array[city-1];
	int toReceive, i, n;
	int args[DAY_STATS_FIELDS + MAX_CITIES];
	int foreigner, host, host_city;
	struct imported_infection *inf;

	while (recv_int(sock, &toReceive))
	{
		n = command_arg_count(toReceive);
		for (i = 0; i < n; i++)
			if (!recv_int(sock, &args[i]))
				break;
		if (i < n)
			break;
		foreigner = args[0];
		host = args[1];
		host_city = args[2];

		pthread_mutex_lock(&net_lock);
		switch(toReceive)
		{
			case CMD_BARRIER:
				barriers_from[city-1]++;
				break;

			case CMD_ADD_FOREIGN_NEIGHBOR:// global rewiring code
				printf("received host: %d to connect to %d from city %d\n",host,foreigner,host_city);
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors].ID = foreigner;
				level[host/unit].id[host%unit].neighbor[level[host/unit].id[host%unit].myNumNeighbors++].city = host_city;
				level[host/unit].id[host%unit].myNumForNeigh++;
				if(level[host/unit].id[host%unit].myNumForNeigh == 1)
					foreign_array[num_for_neigh++] = host;
				break;

			case CMD_DEL_FOREIGN_NEIGHBOR:
		//		printf("received deletion of host %d's neighbor %d, from city %d\n",host, foreigner, city);
				delForNeigh(host, foreigner, city);
				break;

			case CMD_INFECT: // queued; the main thread tries it after that day's barrier
				inf = (struct imported_infection *) malloc(sizeof(struct imported_infection));
				inf->day = args[0];
				inf->host = args[1];
				inf->H = args[2];
				inf->N = args[3];
				inf->strain = args[4];
				inf->next = NULL;
				if (imported_tail[city-1])
					imported_tail[city-1]->next = inf;
				else
					imported_head[city-1] = inf;
				imported_tail[city-1] = inf;
				break;

			case CMD_DAY_STATS:
				memcpy(city_stats[city-1], args, sizeof(int) * n);
				break;

			default:
				printf("ERROR, Unknown Function %d from city %d\n",toReceive,city);
				break;
		}
		pthread_cond_broadcast(&net_cond);
		pthread_mutex_unlock(&net_lock);
	}

	close(sock);
	return NULL;
}//get_commands

/* Blocks until every city has reached this barrier.  Each connection's
   commands are applied in order, so once a city's barrier has arrived,
   everything it sent before the barrier has been applied here too. */
void barrier()
{
	int city, waiting;

	barriers_sent++;
	for (city = 1; city <= num_Cities; city++)
		if (city != CITY_NUMBER)
			send_code(city, CMD_BARRIER);

	pthread_mutex_lock(&net_lock);
	do
	{
		waiting = 0;
		for (city = 1; city <= num_Cities; city++)
			if (city != CITY_NUMBER && barriers_from[city-1] < barriers_sent)
				waiting = 1;
		if (waiting)
			pthread_cond_wait(&net_cond, &net_lock);
	} while (waiting);
	pthread_mutex_unlock(&net_lock);
}//barrier


void runProg(char **argv)
{
	int j;
	
	readvars();
	fprintf(flu_sum_output,"\nNet_type,swnP,vacc_strats,per_vaccs,popsize,NCR,mut_rate,num_inf,tot_inf,duration,peak_inf,peak_day,half_pop_inf_day,# strains\n");
	if(!vaccinate) {
		printf("WARNING: Not vaccinating but summary output will report vacc_strat!\n");
		printf("\tAlthough slower, consider vaccinating but setting proportion \n\tof pop to vacc = 0.0%%\n\n");
	}
	
	swnP = p_values[0];
	nodes = pop_sizes[0];
	strategy = 0;
	vacc_strat = vacc_strats[strategy];
	percentage = 0;
	per_vaccinate = pc_vacc[percentage];
	NCR = NCRs[0];
	mut_rate = mutation_rates[0];
	//	printf("Population: %d\n",nodes);
	foreign_array = (int *) malloc(sizeof(int) * nodes); //worst case: every host has a foreign neighbor
	if (CITY_NUMBER == 1)
		write_epidemic_meta();
	initialize(argv);  /* set up model and do statistics */
	printf("done setting up network and rewiring\n");

	duration = 0;
	peak_num = 0;
	peak_day = 0;
	currentvacc = 0;
	vacc_day = vacc_days[currentvacc++];
	for (day = 0; day < num_days; day++) 
	{
		if(vaccinate==1 && day==vacc_day)
		{
			vaccinate_hosts(); /* vaccinates if its a vacc_day */
			vacc_day = vacc_days[currentvacc];
		}
		recovery(); /* recovery check of all hosts */
		spreadflu(); /* infections of people in other cities are sent to those cities */
		barrier(); /* every city has sent today's infections */
		apply_imported_infections();
		statistics();
		report_day();
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

	printf("NORMAL EXECUTION\n");
	makeGraphicFile();
	exit(0);
}

void send_infection(int city, int host, struct flu *f)
{
	int msg[6] = { CMD_INFECT, day, host, f->H, f->N, f->strain };
	send_ints(city, msg, 6);
}

/* Tries the infections other cities sent us today, in city order so a run is
   repeatable no matter which city's messages arrived first. */
void apply_imported_infections()
{
	int city;
	struct imported_infection *inf;
	struct flu *f;

	for (city = 1; city <= num_Cities; city++)
	{
		while (1)
		{
			pthread_mutex_lock(&net_lock);
			inf = imported_head[city-1];
			if (inf && inf->day == day)
			{
				imported_head[city-1] = inf->next;
				if (!inf->next)
					imported_tail[city-1] = NULL;
			}
			else
				inf = NULL; // nothing left from today (a faster city may already be sending tomorrow's)
			pthread_mutex_unlock(&net_lock);
			if (!inf)
				break;

			f = new_flu_strain();
			f->H = inf->H;
			f->N = inf->N;
			f->strain = inf->strain;
			if (try_infect(inf->host, f))
				imports_today[city-1]++;
			free(inf);
		}
	}
}//apply_imported_infections

/* Counts today's infected/recovered/vaccinated hosts for the summary file, and
   gets them to city 1, which writes one line per day for all cities. */
void report_day()
{
	static const char *field_names[DAY_STATS_FIELDS] = { "day", "infected", "recovered", "vaccinated" };
	int si, city, from, field, waiting, any_imports = 0;
	int numinfhosts = 0, numrechosts = 0, numinfections = 0, numvaccinated = 0;
	int stats[1 + DAY_STATS_FIELDS + MAX_CITIES];
	struct flu *recovereds;

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
				for (recovereds = level[si/unit].id[si%unit].flu_strains_recovered; recovereds; recovereds = recovereds->next)
					if (recovereds->vaccine)
					{
						numvaccinated++;
						break;
					}
		}
	}
	totalinf[day]= numinfhosts;  /* add infections from current run to total */
	totalrec[day]= numrechosts;  /* add recovereds from current run to total */
	totalinfections[day]= numinfections;  /* add count all previous infections from current run to total */
	numstrain[day]= count_strains;  /* add count of strains/day from current run to total */
	totalstrain[day]=  tot_count_strains;  /* add cummulative count of strains from current run to total */

	stats[0] = CMD_DAY_STATS;
	stats[1] = day;
	stats[2] = numinfhosts;
	stats[3] = numrechosts;
	stats[4] = numvaccinated;
	memcpy(&stats[5], imports_today, sizeof(int) * num_Cities);
	memset(imports_today, 0, sizeof(imports_today));

	if (CITY_NUMBER != 1)
	{
		send_ints(1, stats, 1 + DAY_STATS_FIELDS + num_Cities);
		return;
	}

	pthread_mutex_lock(&net_lock);
	memcpy(city_stats[0], &stats[1], sizeof(int) * (DAY_STATS_FIELDS + num_Cities));
	do
	{
		waiting = 0;
		for (city = 2; city <= num_Cities; city++)
			if (city_stats[city-1][0] != day)
				waiting = 1;
		if (waiting)
			pthread_cond_wait(&net_cond, &net_lock);
	} while (waiting);

	fprintf(epidemic_output, "{\"type\":\"day\",\"day\":%d", day);
	for (field = 1; field < DAY_STATS_FIELDS; field++) // stats[0] is the day
	{
		fprintf(epidemic_output, ",\"%s\":[", field_names[field]);
		for (city = 0; city < num_Cities; city++)
			fprintf(epidemic_output, "%s%d", city ? "," : "", city_stats[city][field]);
		fprintf(epidemic_output, "]");
	}
	fprintf(epidemic_output, ",\"imports\":[");
	for (city = 0; city < num_Cities; city++)
		for (from = 0; from < num_Cities; from++)
			if (city_stats[city][DAY_STATS_FIELDS + from] > 0)
			{
				fprintf(epidemic_output, "%s[%d,%d,%d]", any_imports ? "," : "", from, city, city_stats[city][DAY_STATS_FIELDS + from]);
				any_imports = 1;
			}
	fprintf(epidemic_output, "]}\n");
	fflush(epidemic_output);
	pthread_mutex_unlock(&net_lock);
}//report_day

/* First line of epidemic.jsonl: what the visualization needs to lay out the cities. */
void write_epidemic_meta()
{
	int city;
	char *c;

	if ((epidemic_output = fopen("epidemic.jsonl", "w")) == NULL)
	{
		perror("epidemic.jsonl");
		exit(1);
	}
	for (city = 0; city < num_Cities; city++) // keep names JSON-safe
		for (c = city_name[city]; *c; c++)
			if (*c == '"' || *c == '\\' || (unsigned char) *c < ' ')
				*c = ' ';
	fprintf(epidemic_output, "{\"type\":\"meta\",\"version\":1,\"days\":%d,\"population\":%d,\"seedCity\":%d,\"cities\":[",
		num_days, nodes, seed_city - 1);
	for (city = 0; city < num_Cities; city++)
		fprintf(epidemic_output, "%s{\"name\":\"%s\",\"lat\":%.4f,\"lon\":%.4f}",
			city ? "," : "", city_name[city], city_lat[city], city_lon[city]);
	fprintf(epidemic_output, "]}\n");
	fflush(epidemic_output);
}//write_epidemic_meta

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
	rewire();
	barrier(); // every city finishes local rewiring before foreign edges arrive

	/* Cities add and then rewire their cross-city edges one at a time, so edges
	   always arrive in the same order and the same seeds build the same network. */
	for(i = 1; i <= num_Cities; i++)
	{
		if(i == CITY_NUMBER && CITY_NUMBER%2 != 0)
			global_init();
		barrier();
	}
	printf("ALL CITIES GLOBALLY INIT'D\n");

	for(i = 1; i <= num_Cities; i++)
	{
		if(i == CITY_NUMBER)
			global_rewire();
		barrier();
	}
	printf("ALL CITIES GLOBALLY REWIRED\n");

//	print_people();

	printf("CITY %d has %d hosts with foreign neighbors\n",CITY_NUMBER,num_for_neigh);
	//sleep(10);

	clustercoeff();
	initInfect();
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
	int num_crosstown_neighbors = n_cross_city_edges/2;
	
	pthread_mutex_lock(&net_lock);
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
				send_cmd(foreign_city+1, CMD_ADD_FOREIGN_NEIGHBOR, host, foreigner, host_city);
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
				send_cmd(foreign_city+1, CMD_ADD_FOREIGN_NEIGHBOR, host, foreigner, host_city);
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
	pthread_mutex_unlock(&net_lock);
}//global_init

void global_rewire()
{
	int i,j,k,l,m, temp, host;
	float randnumb;
	int randHost;
	int randCity;

	pthread_mutex_lock(&net_lock);
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
				while(isAForNeigh(host,randHost,randCity))			//if it is us or one of our neighbors....
					randHost = nodes*(rand()/(float)RAND_MAX);			//	then get get a different random host
				printf("it will be rewired to %d in city %d\n",randHost,randCity);
				//delForNeigh(host,level[j/unit].id[j%unit].neighbor[k].ID,level[j/unit].id[j%unit].neighbor[k].city);	//delete the host from the original neighbors array
				
				send_cmd(level[host/unit].id[host%unit].neighbor[k].city, CMD_DEL_FOREIGN_NEIGHBOR,
					host, level[host/unit].id[host%unit].neighbor[k].ID, CITY_NUMBER);
				
				level[host/unit].id[host%unit].neighbor[k].ID=0-randHost;		//set the neighbor to the negative value of the random host, negative so that we know it has been rewired already
				level[host/unit].id[host%unit].neighbor[k].city = randCity;
				//add the host to the end of the neighbor hosts array
				send_cmd(randCity, CMD_ADD_FOREIGN_NEIGHBOR, host, randHost, CITY_NUMBER);
			}
			
		}
	}
	pthread_mutex_unlock(&net_lock);
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
	
	if(found == 0)
	{
		printf("\n\n\n\n\n%d's NEIGHBOR %d from city %d NOT FOUND!!\n\n\n\n\n\n\n",host,foreigner,city);
		return;
	}
		
	level[host/unit].id[host%unit].myNumNeighbors--;
	level[host/unit].id[host%unit].myNumForNeigh--;
	
	if(level[host/unit].id[host%unit].myNumForNeigh == 0)
	{
		for(k = 0; foreign_array[k] != host; k++)
			;
		for(; k < num_for_neigh-1; k++)
			foreign_array[k] = foreign_array[k+1];
		num_for_neigh--;
	}
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
	fscanf(input, "%ld %*s", &numbits_strain);
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
	fprintf(flu_sum_output,"numbits_strain = %ld\n", numbits_strain);
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
	int q, r, city;

	/* cities take turns appending their part of the network to the shared file */
	for(city = 1; city <= num_Cities; city++)
	{
		if(city == CITY_NUMBER)
		{
			if ((network = fopen("network.dat", city == 1 ? "w" : "a")) == NULL)
			{
				perror("network.dat");
				exit(1);
			}
			if(city == 1)
				fprintf(network,"%d\n",num_Cities);
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
		}
		barrier();
	}
	printf("ALL CITIES wrote to graphics file\n");
}//makeGraphicFile

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
	
	
	for (i = 0; i < num_infect && (seed_city == 0 || seed_city == CITY_NUMBER); i++) 
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
	int i, j, host, eligible = 0;
	num_vaccinate = (int)(per_vaccinate*nodes);
	for (i = 0; i < nodes; i++) // only never-infected, unvaccinated hosts can get a shot
		if (!level[i/unit].id[i%unit].flu_strains_recovered && !level[i/unit].id[i%unit].flu_strains_infected)
			eligible++;
	if (num_vaccinate > eligible)
		num_vaccinate = eligible; // otherwise the strategies below search forever for someone to vaccinate
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
   far_neigh = -1;
   while (nb < level[h/unit].id[h%unit].myNumNeighbors)
   { /* goes through all the neighbors to find neighbors far away */
		if (level[h/unit].id[h%unit].neighbor[nb].city != CITY_NUMBER)
		{
			nb++;
			continue;
		}
		neighID = abs(level[h/unit].id[h%unit].neighbor[nb].ID);
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
			far_neigh = neighID;
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
			if (level[i/unit].id[i%unit].neighbor[nb].city != CITY_NUMBER)
			{
				nb++;
				continue;
			}
			neighID = abs(level[i/unit].id[i%unit].neighbor[nb].ID);
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
	struct flu *f;

	if (!level[h/unit].id[h%unit].flu_strains_infected) 
	{
		printf("problem in spread_flu_from_host() function - host should be infected\n");
		exit(0);
	}
	/* first check if neighbor can get this flu strain */
	while (count < level[h/unit].id[h%unit].myNumNeighbors) 
	{ /* while h has a neighbor (nb) in list */
		nb = abs(level[h/unit].id[h%unit].neighbor[count].ID);
		
		f = new_flu_strain(); /* ASSUMING THIS FLU STRAIN WILL BE ADDED - MIGHT NOT BE */
		f->H = level[h/unit].id[h%unit].flu_strains_infected->H;
		f->N = level[h/unit].id[h%unit].flu_strains_infected->N;
		f->strain = level[h/unit].id[h%unit].flu_strains_infected->strain;

		if (level[h/unit].id[h%unit].neighbor[count].city == CITY_NUMBER)
			try_infect(nb, f);
		else
		{
			/* neighbor lives in another city: that city tries the infection after today's barrier */
			send_infection(level[h/unit].id[h%unit].neighbor[count].city, nb, f);
			free(f);
		}
		count++;
	}
}//spread flu from host

/* Tries to give host nb flu strain f, based on the neighbor contact rate and nb's
   immunity from the strains it has recovered from.  Frees f if nb doesn't get it. */
int try_infect(int nb, struct flu *f)
{
	struct flu *recovereds, *recovstrains;
	float Hmatch, temp_match, strainmatch;
	int did_infection = 0;
	float random = 0;

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
	return did_infection;
}//try_infect

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
			if (level[k/unit].id[k%unit].neighbor[l].city != CITY_NUMBER)
				continue;
			nb = abs(level[k/unit].id[k%unit].neighbor[l].ID);
		
			for(m = 0; m < level[nb/unit].id[nb%unit].myNumNeighbors; m++)
			{
				for(n = 0; n < level[k/unit].id[k%unit].myNumNeighbors; n++)
				{
					if(abs(level[nb/unit].id[nb%unit].neighbor[m].ID) == abs(level[k/unit].id[k%unit].neighbor[n].ID) &&
						level[nb/unit].id[nb%unit].neighbor[m].city == level[k/unit].id[k%unit].neighbor[n].city) 
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
