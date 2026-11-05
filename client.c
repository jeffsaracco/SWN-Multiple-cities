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
		if (scanf("%d", &host) != 1)
			break;

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

	close(sockfd);
	return 0;
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
