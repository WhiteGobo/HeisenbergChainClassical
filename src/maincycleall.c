#include <stdio.h>
#include <string.h> //für strlen
#include <stdlib.h> //für EXIT_SUCCESS
#include <pthread.h>

static float pi = 3.141592653;


int size = 18000;
float timestep = 0.01;
float J = -1.0;
float Delta = 1.5;
int timemax;
int k = 1000;
int MaxId = 200;
float zcouplingmax = 0.0;
float q;

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
	chain = (struct spinchain *)create_spinchain(&size, &timestep, &J, &Delta, &q, (int*)id, zcouplingmax);
	//print_chain(chain);
	//printmode_chain(chain, &q);
	//printforce_chain(chain);
	//plot_chain(chain);
	plotmode_chain(chain, &q);
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
	return 0;
}


/* Here we take a random derivation of spins in classical dynamic
 * Then we process them in time and look upon their result(fourrierAnalysis)
 * TODO: Randomnumbergenerator; improve process of time; save and plot
 */
int main (int argc, char **argv)
{
	int i;
	int id[MaxId];

	float q;

	request_data();

	timemax = 1000;
	timemax = timemax / timestep;
	q = 2 * pi * k / size;

	pthread_t calcChain_Thread[MaxId];
	for(i=0;i < MaxId;i++)
	{
		printf("Simulation: %d\n", i);
		id[i]=i;
		pthread_create (calcChain_Thread+i, NULL, calculate_spinchain, id+i);
		pthread_detach (*(calcChain_Thread+i));

	}
	return EXIT_SUCCESS;
}
