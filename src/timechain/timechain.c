#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct spinchain {
	float *spins;
	float *spinssecond;
	int size;
	float timestep;
	float time;
	float J; //dipolecoupling constant
	float Delta; //Inhomogenity
	int id;
	struct inhalt *plotter;
	float qanalysis;
	float q;
};

void progress_rk(struct spinchain *chain);
void progress_eul(struct spinchain *chain);
void print_chain(struct spinchain *chain);
void printmode_chain(struct spinchain *chain, float *q);
void plotmodebegin_chain(struct spinchain *chain, float *q);
void printforce_chain(struct spinchain *chain);
struct spinchain *create_spinchain(int *size, float *timestep, float *J, float *Delta, float *q, int *id);
void free_spinchain(struct spinchain *chain);
float beginningx(float q, float  a, float  alpha, float A, int r);
float beginningy(float q, float  a, float  alpha, float A, int r);
float beginningz(float q, float  a, float  alpha, float A, int r);
float timedex(float y, float z, float y2, float z2, float y3, float z3, float J, float Delta);
float timedey(float x, float z, float x2, float z2, float x3, float z3, float J, float Delta);
float timedez(float x, float y, float x2, float y2, float x3, float y3, float J, float Delta);
void plotmode_chain(struct spinchain *chain, float *qanalysis);
void plotmodecycle_chain(struct spinchain *chain);
void plotmodeend_chain(struct spinchain *chain);


/* One step further with Runge-Kutta of S_r
 * timede = Dt S_r = J^r_- - J^r_+
 * J^r_+- = +- J S_r x (s^x_r-1 , s^y_r-1, Delta s^z_r-1)
 * x = S^x_r     y = S^y_r     z = S^z_r
 * x2= S^x_r-1   y2= S^y_r-1
 * x3= S^x_r+1   
 */
void progress_rk(struct spinchain *chain)
{
	int i;
	float sigma[9];
	float erg[chain->size * 3];
	float x,y,z,xm,ym,zm,xp,yp,zp,J,Delta, h;
	float spinlaenge;

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

		sigma[0] = timedex(y,z,ym,zm,yp,zp,J,Delta);
		sigma[1] = timedey(x,z,xm,zm,xp,zp,J,Delta);
		sigma[2] = timedez(x,y,xm,ym,xp,yp,J,Delta);
		x=x + 0.5*h*sigma[0];
		y=y + 0.5*h*sigma[1];
		z=z + 0.5*h*sigma[2];
		sigma[3] = timedex(y,z,ym,zm,yp,zp,J,Delta);
		sigma[4] = timedey(x,z,xm,zm,xp,zp,J,Delta);
		sigma[5] = timedez(x,y,xm,ym,xp,yp,J,Delta);
		x=x - 1.5*h*sigma[0] + 2*h*sigma[3];
		y=y - 1.5*h*sigma[1] + 2*h*sigma[4];
		z=z - 1.5*h*sigma[2] + 2*h*sigma[5];
		sigma[6] = timedex(y,z,ym,zm,yp,zp,J,Delta);
		sigma[7] = timedey(x,z,xm,zm,xp,zp,J,Delta);
		sigma[8] = timedez(x,y,xm,ym,xp,yp,J,Delta);
		erg[i*3]=*(chain->spins+i*3) + h*(sigma[0]+4*sigma[3]+sigma[6])/6;
		erg[i*3+1]=*(chain->spins+i*3+1) + h*(sigma[1]+4*sigma[4]+sigma[7])/6;
		erg[i*3+2]=*(chain->spins+i*3+2) + h*(sigma[2]+4*sigma[5]+sigma[8])/6;
		spinlaenge=sqrtf(erg[i*3+2]*erg[i*3+2] + erg[i*3+1]*erg[i*3+1] + erg[i*3]*erg[i*3]);
		//--erg[i*3] = erg[i*3] / spinlaenge;
		//erg[i*3+1] = erg[i*3+1] / spinlaenge;
		//erg[i*3+2] = erg[i*3+2] / spinlaenge;
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
	
	sigma[0] = timedex(y,z,ym,zm,yp,zp,J,Delta);
	sigma[1] = timedey(x,z,xm,zm,xp,zp,J,Delta);
	sigma[2] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x + 0.5*h*sigma[0];
	y=y + 0.5*h*sigma[1];
	z=z + 0.5*h*sigma[2];
	sigma[3] = timedex(y,z,ym,zm,yp,zp,J,Delta);
	sigma[4] = timedey(x,z,xm,zm,xp,zp,J,Delta);
	sigma[5] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x - 1.5*h*sigma[0] + 2*h*sigma[3];
	y=y - 1.5*h*sigma[1] + 2*h*sigma[4];
	z=z - 1.5*h*sigma[2] + 2*h*sigma[5];
	sigma[6] = timedex(y,z,ym,zm,yp,zp,J,Delta);
	sigma[7] = timedey(x,z,xm,zm,xp,zp,J,Delta);
	sigma[8] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	erg[0]=*(chain->spins) + h*(sigma[0]+4*sigma[3]+sigma[6])/6;
	erg[1]=*(chain->spins+1) + h*(sigma[1]+4*sigma[4]+sigma[7])/6;
	erg[2]=*(chain->spins+2) + h*(sigma[2]+4*sigma[5]+sigma[8])/6;
	spinlaenge=sqrtf(erg[0]*erg[0] + erg[1]*erg[1] + erg[2]*erg[2]);
	//erg[0] = erg[0] / spinlaenge;
	//erg[1] = erg[1] / spinlaenge;
	//erg[2] = erg[2] / spinlaenge;
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

	sigma[0] = timedex(y,z,ym,zm,yp,zp,J,Delta);
	sigma[1] = timedey(x,z,xm,zm,xp,zp,J,Delta);
	sigma[2] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x + 0.5*h*sigma[0];
	y=y + 0.5*h*sigma[1];
	z=z + 0.5*h*sigma[2];
	sigma[3] = timedex(y,z,ym,zm,yp,zp,J,Delta);
	sigma[4] = timedey(x,z,xm,zm,xp,zp,J,Delta);
	sigma[5] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	x=x - 1.5*h*sigma[0] + 2*h*sigma[3];
	y=y - 1.5*h*sigma[1] + 2*h*sigma[4];
	z=z - 1.5*h*sigma[2] + 2*h*sigma[5];
	sigma[6] = timedex(y,z,ym,zm,yp,zp,J,Delta);
	sigma[7] = timedey(x,z,xm,zm,xp,zp,J,Delta);
	sigma[8] = timedez(x,y,xm,ym,xp,yp,J,Delta);
	erg[chain->size*3-3]=*(chain->spins+chain->size*3-3) + h*(sigma[0]+4*sigma[3]+sigma[6])/6;
	erg[chain->size*3-2]=*(chain->spins+chain->size*3-2) + h*(sigma[1]+4*sigma[4]+sigma[7])/6;
	erg[chain->size*3-1]=*(chain->spins+chain->size*3-1) + h*(sigma[2]+4*sigma[5]+sigma[8])/6;
	spinlaenge=sqrtf(erg[chain->size*3-3]*erg[chain->size*3-3] + erg[chain->size*3-2]*erg[chain->size*3-2] + erg[chain->size*3-1]*erg[chain->size*3-1]);
	//erg[chain->size*3-3] = erg[chain->size*3-3] / spinlaenge;
	//erg[chain->size*3-2] = erg[chain->size*3-2] / spinlaenge;
	//erg[chain->size*3-1] = erg[chain->size*3-1] / spinlaenge;
	*(chain->spinssecond+chain->size*3-3) = erg[chain->size*3-3] / spinlaenge;
	*(chain->spinssecond+chain->size*3-2) = erg[chain->size*3-2] / spinlaenge;
	*(chain->spinssecond+chain->size*3-1) = erg[chain->size*3-1] / spinlaenge;
	
	//float *spin = chain->spins;
	//now lets copy our list  when using erg
	//for (i=0; i<chain->size*3; i++)
	//{
	//	*(spin++) = erg[i];
	//}

	float *spinhelp = chain->spins;
	chain->spins = chain->spinssecond;
	chain->spinssecond = spinhelp;
}



/* One step further with Euler of S_r
 * timede = Dt S_r = J^r_- - J^r_+
 * J^r_+- = +- J S_r x (s^x_r-1 , s^y_r-1, Delta s^z_r-1)
 * x = S^x_r     y = S^y_r     z = S^z_r
 * x2= S^x_r-1   y2= S^y_r-1
 * x3= S^x_r+1   
 */
void progress_eul(struct spinchain *chain)
{
	int i;
	float erg[chain->size * 3];
	float x,y,z,x2,y2,z2,x3,y3,z3,J,Delta, h;

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

		erg[i*3]=*(chain->spins+i*3) + h * timedex(y,z,y2,z2,y3,z3,J,Delta);
		erg[i*3+1]=*(chain->spins+i*3+1) + h*timedey(x,z,x2,z2,x3,z3,J,Delta);
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
	
	erg[0]=*(chain->spins) + h*timedex(y,z,y2,z2,y3,z3,J,Delta);
	erg[1]=*(chain->spins+1) + h*timedey(x,z,x2,z2,x3,z3,J,Delta);
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

	erg[chain->size*3-3]=*(chain->spins+i*3) + h*timedex(y,z,y2,z2,y3,z3,J,Delta);
	erg[chain->size*3-2]=*(chain->spins+i*3+1) + h*timedey(x,z,x2,z2,x3,z3,J,Delta);
	erg[chain->size*3-1]=*(chain->spins+i*3+2) + h*timedez(x,y,x2,y2,x3,y3,J,Delta);
	
	float *spin = chain->spins;
	//now lets copy our list
	for (i=0; i<chain->size*3; i++)
	{
		*(spin++) = erg[i];
	}
}



void print_chain(struct spinchain *chain)
{
	int r;
	float *spin = chain->spins;
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
	sprintf(headerstring,"# id:                 %d\n# size:               %d\n# time:               %E\n# coupling constant:  %E\n# inhomogenity Delta: %E\n", chain->id, chain->size, chain->time,chain->J, chain->Delta);
	print_main(data, name, headerstring);
	close_inhalt(data);
}


/* analyse q_z-mode (fourier analysis)
 */
void printmode_chain(struct spinchain *chain, float *q)
{
	float sum = 0.0;
	int r;
	float *spin = chain->spins + 2;
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
void plotmode_chain(struct spinchain *chain, float *qanalysis)
{
	struct inhalt *data;
	char name[255];
	sprintf(name, "output/heismode%.3f_%03d_%.3f_%.2f.dat", *qanalysis, chain->id, chain->q, chain->Delta);
	//printf(name);
	//printf("\n");
	char headerstring[255];
	sprintf(headerstring,"# id:                 %d\n# size:               %d\n# coupling constant:  %E\n# inhomogenity Delta: %E\n", chain->id, chain->size, chain->J, chain->Delta);
	data = (struct inhalt *)init_plotter(name, headerstring);
	chain->plotter = data;
	chain->qanalysis = *qanalysis;
	//print_main(data, name, string);
	//close_inhalt(data);
}


void plotmodecycle_chain(struct spinchain *chain)
{
	float qmode = 0.0;
	int r;
	float *spin = chain->spins + 2;
	for (r=0; r< chain->size; r++)
	{
		qmode= qmode + cosf((chain->qanalysis)*r)* (*spin);
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
	float sum = 0.0;
	int r;
	float *spin = chain->spins;
	for (r=0; r< chain->size; r++)
	{
		sum+= (*spin)*(*spin) + *(spin+1)* *(spin+1) + *(spin+2)* *(spin+2);
		spin= spin+3;
	}
	sum=sum/chain->size;
	printf("%E \n", sum);
}


/* Create a 3dim spin chain with randomdata:
 * TODO: remove alhpa before forschleife, work on A
 *       work on printfs
 */
struct spinchain *create_spinchain(int *size, float *timestep, float *J, float *Delta, float *q, int *id)
{
	float a, alpha;
	struct spinchain *chain = malloc(sizeof(struct spinchain));
	chain->size = *size;

	int *allrandomnumber;
	allrandomnumber = calloc(5+ (2* *size), sizeof(int));
	char name[20];
	sprintf(name, "data/random%03d.dat", *id);
	intsofsize(name,3, 5+(2* *size), allrandomnumber);
	float singlerandomnumber;

	//singlerandomnumber = 0.001 * (float)*(allrandomnumber );
	//float alpha=singlerandomnumber*2*3.141592653;	

	//singlerandomnumber = 0.001 * (float)*(allrandomnumber + 1);
	singlerandomnumber = 0.25;
	float A=singlerandomnumber;	

	int r;
	chain->spins = calloc(3* *size,64);
	chain->spinssecond = calloc(3* *size,64);
	for (r=0; r< *size; r++)
	{
		singlerandomnumber = 0.001 * (float)*(allrandomnumber + 5 + r);
		a = singlerandomnumber;
		singlerandomnumber = 0.001 * (float)*(allrandomnumber + 5 + r + *size);
		alpha = 3.141592653 * 2 * singlerandomnumber;
		*(chain->spins+(r*3)  )=beginningx(*q, a, alpha, A, r);
		*(chain->spins+(r*3)+1)=beginningy(*q, a, alpha, A, r);
		*(chain->spins+(r*3)+2)=beginningz(*q, a, alpha, A, r);
	}
	chain->timestep = *timestep;
	chain->J = *J;
	chain->Delta = *Delta;
	chain->time = 0.0;
	chain->q = *q;
	chain->id = *id;

	//printf("hier kommt size  %d \n", *size);
	//printf("hier kommt q     %f \n", *q);
	//printf("hier kommt J     %f \n", *J);
	//printf("hier kommt A     %f \n", A);
	//printf("hier kommt alpha %f \n", alpha);
	//printf("hier kommt id    %d \n", *id);
	//printf("hier kommt Delta %f \n", chain->Delta);
	
	free(allrandomnumber);
	return chain;
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
//float timedex(float x, float y, float z, float x2, float y2, float z2, float x3, float y3, float z3, float J, float Delta)
float timedex(float y, float z, float y2, float z2, float y3, float z3, float J, float Delta)
{
	//return J * (Delta * y * (z2 - z3) +         z * ( y2 - y3));
	return J * ( z * ( y2 + y3)  -  Delta * y * (z2 + z3) );
} 


float timedey(float x, float z, float x2, float z2, float x3, float z3, float J, float Delta)
{
	//return J * (        z * (x2 - x3) + Delta * x * (y2 - y3));
	return J * (  Delta * x * (z2 + z3) -        z * (x2 + x3) );
}


float timedez(float x, float y, float x2, float y2, float x3, float y3, float J, float Delta)
{
	//return J * (        x * (y2 - y3) +         y * (x2 - x3));
	return J * (  y * (x2 + x3)   -   x * (y2 + y3) );
}


/* To start we need a Anfangsverteilung
 * We need 3 random numbers A, a, alpha and momentum of the system
 * cos alpha * sqrt(1-A^2cos^2(qr) a^2)
 */
float beginningx(float q, float  a, float  alpha, float A, int r)
{
	return cos(alpha) * sqrtf(1 - powf(beginningz( q, a, alpha, A, r),2.0f));
}


float beginningy(float q, float  a, float  alpha, float A, int r)
{
	return sin(alpha) * sqrtf(1 - powf(beginningz( q, a, alpha, A, r),2.0f));
}


float beginningz(float q, float  a, float  alpha, float A, int r)
{
	return  A * cos(q*r) * a;
}