#include <stdio.h>
#include <string.h> //für strlen
#include <stdlib.h> //für EXIT_SUCCESS

static float pi = 3.141592653;

int size, id, k, time;
float timestep, J, Delta, q, a, alpha;



/* Here we take a random derivation of spins in classical dynamic
 * Then we process them in time and look upon their result(fourrierAnalysis)
 * TODO: Randomnumbergenerator; improve process of time; save and plot
 */
int main (int argc, char **argv)
{
	size = 10000;
	timestep = 0.01;
	Delta = 1.5;
	J = 1;
	time = 1000;
	id = 0;
	k = 20;

	//printf("Wir berechnen eine Spin-Kette:\n");
	
	//int size, id, k;
	//float timestep, J, Delta, q, a, alpha;

	//printf("Bitte geben sie die Größe des Systems ein: ");
	//scanf("%d", &size);
	//printf("\n");

	//printf("Bitte geben sie die Zeitschrittlänge ein: ");
	//scanf("%f", &timestep);
	//printf("\n");

	//printf("Bitte geben sie die Inhomogenität ein: ");
	//scanf("%f", &Delta);
	//printf("\n");

	//printf("Bitte geben sie die Stärke der Dipolkopplung ein: ");
	//scanf("%f", &J);
	//printf("\n");

	k = 10.0;
	q = 2 * pi * k / size;
	//A = 5.0;
	a = 0.6;
	alpha = 1.1;

	//int time;
	//printf("Bitte geben sie die Dauer der Simulation ein: ");
	//scanf("%d", &time);
	//printf("\n");
	time = time / timestep;

	//printf("Bitte geben sie die Id für die Zufallszahlen ein: ");
	//scanf("%d", &id);
	//printf("%d\n", id);


	//printf("Bitte geben sie die Wellenlänge ein: ");
	//scanf("%d", &k);
	//printf("%d, %daaaaaaaa \n", id, k);
	//q = 2 * pi * k / size;

	struct spinchain *chain;
	chain = (struct spinchain *)create_spinchain(&size, &timestep, &J, &Delta, &q, &id);
	//print_chain(chain);
	//printmode_chain(chain, &q);
	//printforce_chain(chain);
	plot_chain(chain);
	//plotmode_chain(chain, &q);
	int i;
	for (i=0; i< time; i++){
		//progress_rk(chain);
		//printmode_chain(chain, &q);
		//progress_eul(chain);
		//if(i%10 == 0) plotmodecycle_chain(chain);
		//if(i%10 == 0) printforce_chain(chain);
		//if(i%50 == 0) plot_chain(chain);
		//if(i%5000 == 0) printmode_chain(chain, &q);
	}
	//print_chain(chain);
	//printmode_chain(chain, &q);
	//printforce_chain(chain);
	//plot_chain(chain);
	//plotmodeend_chain(chain);
	free_spinchain(chain);
	return EXIT_SUCCESS;
}

