#include <stdio.h>
#include <string.h> //für strlen
#include <stdlib.h> //für EXIT_SUCCESS

static float pi = 3.141592653;


int size = 18000;
float timestep = 0.01;
float J = -1.0;
float Delta = 1.5;
int time, id;
int k = 1000;
int MaxId = 200;

/* Read data of the Heisenbergchain
 */
void request_data (void)
{
	scanf("%d", &size);
	scanf("%f", &timestep);
	scanf("%f", &Delta);
	scanf("%f", &J);
	scanf("%d", &time);
	scanf("%d", &k);
}


/* Here we take a random derivation of spins in classical dynamic
 * Then we process them in time and look upon their result(fourrierAnalysis)
 * TODO: Randomnumbergenerator; improve process of time; save and plot
 */
int main (int argc, char **argv)
{
	int i;
	float q;
	struct spinchain *chain;

	request_data();

	time = 1000;
	time = time / timestep;
	q = 2 * pi * k / size;

	id = 0;
	while(id < MaxId)
	{
		chain = (struct spinchain *)create_spinchain(&size, &timestep, &J, &Delta, &q, &id);
		//print_chain(chain);
		//printmode_chain(chain, &q);
		//printforce_chain(chain);
		//plot_chain(chain);
		plotmode_chain(chain, &q);
		for (i=0; i< time; i++){
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
		id++;
	}
	return EXIT_SUCCESS;
}

