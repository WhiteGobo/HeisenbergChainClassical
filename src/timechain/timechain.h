struct spinchain {
	float *spins;
	float *spinssecond;
	float *randomzcoupling;
	int size;
	float timestep;
	float time;
	float J; //dipolecoupling constant
	float Delta; //Inhomogenity
	int id;
	struct inhalt *plotter;
	float qanalysis;
	float q;
	int qnumber;
	float max_randomzcoupling;
};

void progress_rk(struct spinchain *chain);
void progress_eul(struct spinchain *chain);
void print_chain(struct spinchain *chain);
void printmode_chain(struct spinchain *chain, float *q);
void plotmodebegin_chain(struct spinchain *chain, float q);
void printforce_chain(struct spinchain *chain);
struct spinchain *create_spinchain(int size, float timestep, float J, float Delta, float q, int *id, float zcouplingmax);
struct spinchain *create_spinchain2(int size, float timestep, float J, float Delta, int qnumber, int *id, float zcouplingmax);
void free_spinchain(struct spinchain *chain);
float beginningx(float q, float  a, float  alpha, float A, int r);
float beginningy(float q, float  a, float  alpha, float A, int r);
float beginningz(float q, float  a, float  alpha, float A, int r);
float beginningx2(float cosqr, float  a, float  alpha, float A, int r);
float beginningy2(float cosqr, float  a, float  alpha, float A, int r);
float beginningz2(float cosqr, float  a, float  alpha, float A, int r);
float timedex(float y, float z, float y2, float z2, float y3, float z3, float J, float Delta, float zcoupling);
float timedey(float x, float z, float x2, float z2, float x3, float z3, float J, float Delta, float zcoupling);
float timedez(float x, float y, float x2, float y2, float x3, float y3, float J, float Delta);
void plotmode_chain(struct spinchain *chain, float qanalysis);
void plotmodecycle_chain(struct spinchain *chain);
void plotmodeend_chain(struct spinchain *chain);
