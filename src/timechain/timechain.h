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
