/* ------------------ inet.h -------------------------
 * Definitions for TCP and UDP client / server programs
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <errno.h>
#include <netdb.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

/* Used by the two-city swn.c and the client.c test harness.  The
   multi-city swn_cities.c reads city addresses from cities.txt instead. */
//#define SERV_HOST_ADDR  "137.238.7.32"
#define SERV_HOST_ADDR  "127.0.0.1"
#define SERV_TCP_PORT   5529


#define MAXLINE 256
