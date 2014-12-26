#include <stdio.h>
#include <string.h> //für strlen
#include <stdlib.h> //für EXIT_SUCCESS
#include <pthread.h>
#include <sys/resource.h>
#include <errno.h>
#include "timechain/timechain.h"

static float pi = 3.141592653;


int size = 18000;
float timestep = 0.01;
float J = -1.0;
float Delta = 1.5;
int timemax;
int k = 1000;
//int MaxId = 1000;
int MaxId = 10;
float zcouplingmax = 0.0;
float q;
int qnumber;
int runningThreads = 0;
pthread_mutex_t lock_runningThreads;
int max_Threads = 8;
int waitingtime = 10;

/* Read data of the Heisenbergchain
 */
void request_data (void)
{
	scanf("%d", &size);
	scanf("%f", &timestep);
	scanf("%f", &Delta);
	scanf("%f", &J);
	scanf("%d", &timemax);
	scanf("%d", &k);
	scanf("%f", &zcouplingmax);
}


void *calculate_spinchain(void *id)
{
	int i;
	struct spinchain *chain;
	printf("halo: %f\n", timestep);
	chain = (struct spinchain *)create_spinchain2(size, timestep, J, Delta, qnumber, (int*)id, zcouplingmax);
	if(chain != NULL){
		//print_chain(chain);
		//printmode_chain(chain, &q);
		//printforce_chain(chain);
		//plot_chain(chain);
		plotmode_chain(chain, q);
		for (i=0; i< timemax; i++){
			progress_rk(chain);
			//printmode_chain(chain, &q);
			//progress_eul(chain);
			if(i%100 == 0) plotmodecycle_chain(chain);
			//if(i%10 == 0) printforce_chain(chain);
			//if(i%3000 == 0) plot_chain(chain);
		}
		//print_chain(chain);
		//printmode_chain(chain, &q);
		//printforce_chain(chain);
		//plot_chain(chain);
		plotmodeend_chain(chain);
		free_spinchain(chain);
	}
	pthread_mutex_lock(&lock_runningThreads);
	runningThreads--;
	pthread_mutex_unlock(&lock_runningThreads);
	return 0;
}



int checkandincrease_runningThreads(int maximum)
{
	while (pthread_mutex_trylock(&lock_runningThreads) != 0) sleep(waitingtime);
	if (runningThreads > maximum){
		pthread_mutex_unlock(&lock_runningThreads);
		return -1;
	}
	runningThreads++;
	pthread_mutex_unlock(&lock_runningThreads);
	return 0;
}




void getlimit(void)
{
	struct rlimit info;
	if(getrlimit(RLIMIT_AS, &info) != 0){
		fputs("failed to load Limits\n", stderr);
		fputs(strerror(errno), stderr);
	}
	printf("Limits of memory: %d\n", (long)info.rlim_max);
}

/* Here we take a random derivation of spins in classical dynamic
 * Then we process them in time and look upon their result(fourrierAnalysis)
 * TODO: Randomnumbergenerator; improve process of time; save and plot
 */
int main (int argc, char **argv)
{
	int i;
	int id[MaxId];

	request_data();

	timemax = timemax / timestep;
	q = 2 * pi * k / size;
	qnumber = (int)size / k;

	pthread_t calcChain_Thread[MaxId];
	for(i=0;i < MaxId;i++)
	{
		id[i]=i;
		pthread_create (calcChain_Thread+i, NULL, calculate_spinchain, id+i);
		pthread_detach (*(calcChain_Thread+i));
		
		while (checkandincrease_runningThreads(max_Threads) != 0) sleep(waitingtime);
	}
	while (checkandincrease_runningThreads(0) != 0) sleep(waitingtime);
	return EXIT_SUCCESS;
}
