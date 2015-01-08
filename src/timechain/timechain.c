#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "../reader/reader.h"
#include "../plotter/plotter.h"

static double pi = 3.141592653;

struct spinchain {
	double *spins;
	double *spinssecond;
	double *randomzcoupling;
	int size;
	double timestep;
	double time;
	double J; //dipolecoupling constant
	double Delta; //Inhomogenity
	int id;
	struct inhalt *plotter;
	double qanalysis;
	double q;
	int qnumber;
	double max_randomzcoupling;
};

void progress_rk(struct spinchain *chain);
void progress_eul(struct spinchain *chain);
void print_chain(struct spinchain *chain);
void printmode_chain(struct spinchain *chain, double *q);
void plotmodebegin_chain(struct spinchain *chain, double q);
void printforce_chain(struct spinchain *chain);
struct spinchain *create_spinchain(int size, double timestep, double J, double Delta, int qnumber, int *id, double zcouplingmax);
void free_spinchain(struct spinchain *chain);
double beginningx(double cosqr, double  a, double  alpha, double A, int r);
double beginningy(double cosqr, double  a, double  alpha, double A, int r);
double beginningz(double cosqr, double  a, double  alpha, double A, int r);
double timedex(double y, double z, double y2, double z2, double y3, double z3, double J, double Delta, double zcoupling);
double timedey(double x, double z, double x2, double z2, double x3, double z3, double J, double Delta, double zcoupling);
double timedez(double x, double y, double x2, double y2, double x3, double y3, double J, double Delta);
void plotmode_chain(struct spinchain *chain, double qanalysis);
void plotmodecycle_chain(struct spinchain *chain);
void plotmodeend_chain(struct spinchain *chain);


/*
 *One step further with Runge-Kutta of S_r
 * timede = Dt S_r = J^r_- - J^r_+
 * J^r_+- = +- J S_r x (s^x_r-1 , s^y_r-1, Delta s^z_r-1)  + J h_r x S_r
 * x = S^x_r     y = S^y_r     z = S^z_r
 * x2= S^x_r-1   y2= S^y_r-1
 * x3= S^x_r+1   
 */
void progress_rk(struct spinchain *chain)
{
	int i;
	double sigma[9];
	double erg[chain->size * 3];
	double x,y,z,xm,ym,zm,xp,yp,zp,J,Delta, h, zcoupling;
	double spinlaenge;

	J=chain->J;
	Delta=chain->Delta;
	h=chain->timestep;
	chain->time = chain->time + h;
	for (i=1; i< chain->size-1; i++)
	{
		x=*(chain->spins+i*3);
		y=*(chain->spins+i*3+1);
		z=*(chain->spins+i*3+2);
		xm=*(chain->spins+i*3-3);
		ym=*(chain->spins+i*3-2);
		zm=*(chain->spins+i*3-1);
		xp=*(chain->spins+i*3+3);
		yp=*(chain->spins+i*3+4);
		zp=*(chain->spins+i*3+5);
		zcoupling=*(chain->randomzcoupling+i);

		sigma[0] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
		sigma[1] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
		sigma[2] = timedez(x,y,xm,ym,xp,yp,J,Delta);
		x=x + 0.5*h*sigma[0];
		y=y + 0.5*h*sigma[1];
		z=z + 0.5*h*sigma[2];
		sigma[3] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
		sigma[4] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
		sigma[5] = timedez(x,y,xm,ym,xp,yp,J,Delta);
		x=x - 1.5*h*sigma[0] + 2*h*sigma[3];
		y=y - 1.5*h*sigma[1] + 2*h*sigma[4];
		z=z - 1.5*h*sigma[2] + 2*h*sigma[5];
		sigma[6] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
		sigma[7] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
		sigma[8] = timedez(x,y,xm,ym,xp,yp,J,Delta);
		erg[i*3]=*(chain->spins+i*3) + h*(sigma[0]+4*sigma[3]+sigma[6])/6;
		erg[i*3+1]=*(chain->spins+i*3+1) + h*(sigma[1]+4*sigma[4]+sigma[7])/6;
		erg[i*3+2]=*(chain->spins+i*3+2) + h*(sigma[2]+4*sigma[5]+sigma[8])/6;
		spinlaenge=sqrtf(erg[i*3+2]*erg[i*3+2] + erg[i*3+1]*erg[i*3+1] + erg[i*3]*erg[i*3]);
		*(chain->spinssecond+i*3  ) = erg[i*3  ] / spinlaenge;
		*(chain->spinssecond+i*3+1) = erg[i*3+1] / spinlaenge;
		*(chain->spinssecond+i*3+2) = erg[i*3+2] / spinlaenge;
	}

	//hier für den ersten und letzten spin berechnet
	//Es gelten nämlich periodische Randbedingungen
	x=*(chain->spins);
	y=*(chain->spins+1);
	z=*(chain->spins+2);
	xm=*(chain->spins+chain->size*3-3);
	ym=*(chain->spins+chain->size*3-2);
	zm=*(chain->spins+chain->size*3-1);
	xp=*(chain->spins+3);
	yp=*(chain->spins+4);
	zp=*(chain->spins+5);
	zcoupling=*(chain->randomzcoupling);
	
	sigma[0] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
	sigma[1] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
	sigma[2] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x + 0.5*h*sigma[0];
	y=y + 0.5*h*sigma[1];
	z=z + 0.5*h*sigma[2];
	sigma[3] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
	sigma[4] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
	sigma[5] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x - 1.5*h*sigma[0] + 2*h*sigma[3];
	y=y - 1.5*h*sigma[1] + 2*h*sigma[4];
	z=z - 1.5*h*sigma[2] + 2*h*sigma[5];
	sigma[6] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
	sigma[7] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
	sigma[8] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	erg[0]=*(chain->spins) + h*(sigma[0]+4*sigma[3]+sigma[6])/6;
	erg[1]=*(chain->spins+1) + h*(sigma[1]+4*sigma[4]+sigma[7])/6;
	erg[2]=*(chain->spins+2) + h*(sigma[2]+4*sigma[5]+sigma[8])/6;
	spinlaenge=sqrtf(erg[0]*erg[0] + erg[1]*erg[1] + erg[2]*erg[2]);
	*(chain->spinssecond ) = erg[0] / spinlaenge;
	*(chain->spinssecond+1) = erg[1] / spinlaenge;
	*(chain->spinssecond+2) = erg[2] / spinlaenge;

	//nur noch für den letzten spin berechnet
	//Es gelten nämlich periodische Randbedingungen
	x=*(chain->spins+chain->size*3-3);
	y=*(chain->spins+chain->size*3-2);
	z=*(chain->spins+chain->size*3-1);
	xm=*(chain->spins+chain->size*3-6);
	ym=*(chain->spins+chain->size*3-5);
	zm=*(chain->spins+chain->size*3-4);
	xp=*(chain->spins);
	yp=*(chain->spins+1);
	zp=*(chain->spins+2);
	zcoupling=*(chain->randomzcoupling+chain->size);

	sigma[0] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
	sigma[1] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
	sigma[2] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x + 0.5*h*sigma[0];
	y=y + 0.5*h*sigma[1];
	z=z + 0.5*h*sigma[2];
	sigma[3] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
	sigma[4] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
	sigma[5] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x - 1.5*h*sigma[0] + 2*h*sigma[3];
	y=y - 1.5*h*sigma[1] + 2*h*sigma[4];
	z=z - 1.5*h*sigma[2] + 2*h*sigma[5];
	sigma[6] = timedex(y,z,ym,zm,yp,zp,J,Delta,zcoupling);
	sigma[7] = timedey(x,z,xm,zm,xp,zp,J,Delta,zcoupling);
	sigma[8] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	erg[chain->size*3-3]=*(chain->spins+chain->size*3-3) + h*(sigma[0]+4*sigma[3]+sigma[6])/6;
	erg[chain->size*3-2]=*(chain->spins+chain->size*3-2) + h*(sigma[1]+4*sigma[4]+sigma[7])/6;
	erg[chain->size*3-1]=*(chain->spins+chain->size*3-1) + h*(sigma[2]+4*sigma[5]+sigma[8])/6;
	spinlaenge=sqrtf(erg[chain->size*3-3]*erg[chain->size*3-3] + erg[chain->size*3-2]*erg[chain->size*3-2] + erg[chain->size*3-1]*erg[chain->size*3-1]);
	*(chain->spinssecond+chain->size*3-3) = erg[chain->size*3-3] / spinlaenge;
	*(chain->spinssecond+chain->size*3-2) = erg[chain->size*3-2] / spinlaenge;
	*(chain->spinssecond+chain->size*3-1) = erg[chain->size*3-1] / spinlaenge;

	double *spinhelp = chain->spins;
	chain->spins = chain->spinssecond;
	chain->spinssecond = spinhelp;
}



/*
 * One step further with Euler of S_r
 * timede = Dt S_r = J^r_- - J^r_+
 * J^r_+- = +- J S_r x (s^x_r-1 , s^y_r-1, Delta s^z_r-1)  + J h_r x S_r
 * x = S^x_r     y = S^y_r     z = S^z_r
 * x2= S^x_r-1   y2= S^y_r-1
 * x3= S^x_r+1   
 */
void progress_eul(struct spinchain *chain)
{
	int i;
	double erg[chain->size * 3];
	double x,y,z,x2,y2,z2,x3,y3,z3,J,Delta, h, zcoupling;

	J=chain->J;
	Delta=chain->Delta;
	h=chain->timestep;
	chain->time = chain->time + h;
	for (i=1; i< chain->size-1; i++)
	{
		x=*(chain->spins+i*3);
		y=*(chain->spins+i*3+1);
		z=*(chain->spins+i*3+2);
		x2=*(chain->spins+i*3-3);
		y2=*(chain->spins+i*3-2);
		z2=*(chain->spins+i*3-1);
		x3=*(chain->spins+i*3+3);
		y3=*(chain->spins+i*3+4);
		z3=*(chain->spins+i*3+5);

		erg[i*3]=*(chain->spins+i*3) + h * timedex(y,z,y2,z2,y3,z3,J,Delta,zcoupling);
		erg[i*3+1]=*(chain->spins+i*3+1) + h*timedey(x,z,x2,z2,x3,z3,J,Delta,zcoupling);
		erg[i*3+2]=*(chain->spins+i*3+2) + h*timedez(x,y,x2,y2,x3,y3,J,Delta);
	}

	//hier für den ersten und letzten spin berechnet
	//Es gelten nämlich periodische Randbedingungen
	x=*(chain->spins);
	y=*(chain->spins+1);
	z=*(chain->spins+2);
	x2=*(chain->spins+chain->size*3-3);
	y2=*(chain->spins+chain->size*3-2);
	z2=*(chain->spins+chain->size*3-1);
	x3=*(chain->spins+3);
	y3=*(chain->spins+4);
	z3=*(chain->spins+5);
	
	erg[0]=*(chain->spins) + h*timedex(y,z,y2,z2,y3,z3,J,Delta,zcoupling);
	erg[1]=*(chain->spins+1) + h*timedey(x,z,x2,z2,x3,z3,J,Delta,zcoupling);
	erg[2]=*(chain->spins+2) + h*timedez(x,y,x2,y2,x3,y3,J,Delta);

	//nur noch für den letzten spin berechnet
	//Es gelten nämlich periodische Randbedingungen
	x=*(chain->spins+chain->size*3-3);
	y=*(chain->spins+chain->size*3-2);
	z=*(chain->spins+chain->size*3-1);
	x2=*(chain->spins+chain->size*3-6);
	y2=*(chain->spins+chain->size*3-5);
	z2=*(chain->spins+chain->size*3-4);
	x3=*(chain->spins);
	y3=*(chain->spins+1);
	z3=*(chain->spins+2);

	erg[chain->size*3-3]=*(chain->spins+i*3) + h*timedex(y,z,y2,z2,y3,z3,J,Delta,zcoupling);
	erg[chain->size*3-2]=*(chain->spins+i*3+1) + h*timedey(x,z,x2,z2,x3,z3,J,Delta,zcoupling);
	erg[chain->size*3-1]=*(chain->spins+i*3+2) + h*timedez(x,y,x2,y2,x3,y3,J,Delta);
	
	double *spin = chain->spins;
	//now lets copy our list
	for (i=0; i<chain->size*3; i++)
	{
		*(spin++) = erg[i];
	}
}



void print_chain(struct spinchain *chain)
{
	int r;
	double *spin = chain->spins;
	printf("Richtung der Spins: \n");
	for (r=0; r<chain->size; r++)
	{
		printf("%E ", *spin++);
		printf("%E ", *spin++);
		printf("%E ", *spin++);
		printf("\n");
	}
}


void plot_chain(struct spinchain *chain)
{
	struct inhalt *data;
	data = (struct inhalt *)init_inhalt();
	ins_numbers(data, chain->spins, 3, chain->size);

	char name[255];
	sprintf(name, "output/heismode%04.0f_%03d_%.3f_%.2f.dat", chain->time, chain->id, chain->q, chain->Delta);
	char headerstring[255];
	sprintf(headerstring,"# id:                 %d\n# size:               %d\n# time:               %E\n# coupling constant:  %E\n# inhomogenity Delta: %E\n# zcouplingmax:     %E\n", chain->id, chain->size, chain->time,chain->J, chain->Delta, chain->max_randomzcoupling);
	print_main(data, name, headerstring);
	close_inhalt(data);
}


/* analyse q_z-mode (fourier analysis)
 */
void printmode_chain(struct spinchain *chain, double *q)
{
	double sum = 0.0;
	int r;
	double *spin = chain->spins + 2;
	for (r=0; r< chain->size; r++)
	{
		sum= sum + (cosf((*q)*r))* (*spin);
		spin= spin+3;
	}
	printf("qmode: %f \n", sum);
	printf("q    : %f \n", *q);
}


/* create a plotter device for continuus plotting
 * uses plotmodecycle and plotmodeend
 */
void plotmode_chain(struct spinchain *chain, double qanalysis)
{
	struct inhalt *data;
	char name[255];
	sprintf(name, "output/heismode%.3f_%03d_%.3f_%.2f.dat", qanalysis, chain->id, chain->q, chain->Delta);
	//printf(name);
	//printf("\n");
	char headerstring[255];
	sprintf(headerstring,"# id:                      %d\n# size:                    %d\n# coupling constant:       %E\n# inhomogenity Delta:      %E\n# maximum randomzcoupling: ", chain->id, chain->size, chain->J, chain->Delta, chain->max_randomzcoupling);
	data = (struct inhalt *)init_plotter(name, headerstring);
	chain->plotter = data;
	chain->qanalysis = qanalysis;
	//print_main(data, name, string);
	//close_inhalt(data);
}


void plotmodecycle_chain(struct spinchain *chain)
{
	double qmode = 0.0;
	int r;
	double *spin = chain->spins + 2;
	for (r=0; r< chain->size; r++)
	{
		qmode= qmode + cosf((double)(2*pi*r%(chain->qnumber) / chain->qnumber))* (*spin);
		spin= spin+3;
	}
	cycle_plotter(chain->plotter, &(chain->time), &(qmode));
}


void plotmodeend_chain(struct spinchain *chain)
{
	close_plotter(chain->plotter);
}


/* analyze if spins are |S| = 1
 */
void printforce_chain(struct spinchain *chain)
{
	double sum = 0.0;
	int r;
	double *spin = chain->spins;
	for (r=0; r< chain->size; r++)
	{
		sum+= (*spin)*(*spin) + *(spin+1)* *(spin+1) + *(spin+2)* *(spin+2);
		spin= spin+3;
	}
	sum=sum/chain->size;
	printf("%E \n", sum);
}




/* Create a 3dim spin chain with randomdata:
 * TODO: remove alpha before forschleife, work on A
 *       work on printfs
 */
struct spinchain *create_spinchain(int size, double timestep, double J, double Delta, int qnumber, int *id, double zcouplingmax)
{
	double a, A, alpha, singlerandomnumber, cosqr;
	int r;
	struct spinchain *chain = malloc(sizeof(struct spinchain));
	char name[30];
	int *allrandomnumber;
	
	chain->size = size;
	if(alloc_spinchain(chain) != 0) return NULL;

	allrandomnumber = calloc(5+ (3* size), sizeof(int));
	if(allrandomnumber == NULL){
		fputs("allokieren von allrandomnumber\n", stderr);
		free(chain->spins);
		free(chain->spinssecond);
		free(chain->randomzcoupling);
		free(chain);
		return NULL;
	}
	sprintf(name, "data/random%03d.dat", *id);
	if(intsofsize(name,3, 5+(3* size), allrandomnumber)<0){
		free(chain->spins);
		free(chain->spinssecond);
		free(chain->randomzcoupling);
		free(chain);
		free(allrandomnumber);
		return NULL;
	}

	//singlerandomnumber = 0.001 * (double)*(allrandomnumber );
	//double alpha=singlerandomnumber*2*3.141592653;	

	//singlerandomnumber = 0.001 * (double)*(allrandomnumber + 1);
	singlerandomnumber = 0.25;
	A=singlerandomnumber;	


	for (r=0; r< size; r++)
	{
		cosqr = cos(2 * pi * (r % qnumber) / qnumber);

		singlerandomnumber = 0.001 * (double)*(allrandomnumber + 5 + r);
		a = singlerandomnumber;
		singlerandomnumber = 0.001 * (double)*(allrandomnumber + 5 + r + size);
		alpha = 3.141592653 * 2 * singlerandomnumber;
		*(chain->spins+(r*3)  )=beginningx(cosqr, a, alpha, A, r);
		*(chain->spins+(r*3)+1)=beginningy(cosqr, a, alpha, A, r);
		*(chain->spins+(r*3)+2)=beginningz(cosqr, a, alpha, A, r);
		singlerandomnumber = 0.001 * (double)*(allrandomnumber + 5 + r + (2* size));
		*(chain->randomzcoupling+r)=((2*singlerandomnumber)-1)*zcouplingmax;
	}
	chain->timestep = timestep;
	chain->J = J;
	chain->Delta = Delta;
	chain->time = 0.0;
	chain->q = 2 * pi / qnumber;
	chain->qnumber = qnumber;
	chain->id = *id;
	chain->max_randomzcoupling = zcouplingmax;

	printf("hier kommt size    %d \n", size);
	printf("hier kommt qnumber %d \n", qnumber);
	printf("hier kommt J       %f \n", J);
	printf("hier kommt A       %f \n", A);
	printf("hier kommt alpha   %f \n", alpha);
	printf("hier kommt id      %d \n", *id);
	printf("hier kommt Delta   %f \n", chain->Delta);
	
	free(allrandomnumber);
	return chain;
}


int alloc_spinchain(struct spinchain *chain)
{
	chain->spins = calloc(3* chain->size,64);
	if(chain->spins == NULL){
		fputs("allokieren von spins", stderr);
		free(chain);
		return -1;
	}
	chain->spinssecond = calloc(3* chain->size,64);
	if(chain->spinssecond == NULL){
		fputs("allokieren von spinssecond", stderr);
		free(chain->spins);
		free(chain);
		return -2;
	}
	chain->randomzcoupling = calloc(chain->size,64);
	if(chain->randomzcoupling == NULL){
		fputs("allokieren von randomzcoupling\n", stderr);
		free(chain->spins);
		free(chain->spinssecond);
		free(chain);
		return -3;
	}
	return 0;
}


void free_spinchain(struct spinchain *chain)
{
	free(chain->spins);
	free(chain->spinssecond);
	free(chain);	
}





/* Calculate timederivation 
 * x = S^x_r     y = S^y_r     z = S^z_r
 * x2= S^x_r-1   y2= S^y_r-1
 * x3= S^x_r+1   
 */
double timedex(double y, double z, double y2, double z2, double y3, double z3, double J, double Delta, double zcoupling)
{
	return J * (        z * ( y2 + y3)  -  Delta * y * (z2 + z3) + y * zcoupling);
} 


double timedey(double x, double z, double x2, double z2, double x3, double z3, double J, double Delta, double zcoupling)
{
	return J * ( Delta * x * (z2 + z3) -            z * (x2 + x3) - x * zcoupling );
}


double timedez(double x, double y, double x2, double y2, double x3, double y3, double J, double Delta)
{
	return J * (         y * (x2 + x3)   -          x * (y2 + y3) );
}



/* To start we need a Anfangsverteilung
 * We need 3 random numbers A, a, alpha and momentum of the system
 * cos alpha * sqrt(1-A^2cos^2(qr) a^2)
 */
double beginningx(double cosqr, double  a, double  alpha, double A, int r)
{
	return cos(alpha) * sqrtf(1 - powf(beginningz( cosqr, a, alpha, A, r),2.0f));
}


double beginningy(double cosqr, double  a, double  alpha, double A, int r)
{
	return sin(alpha) * sqrtf(1 - powf(beginningz( cosqr, a, alpha, A, r),2.0f));
}


double beginningz(double cosqr, double  a, double  alpha, double A, int r)
{
	return  A * cosqr * a;
}


