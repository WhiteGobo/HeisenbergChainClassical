#include<stdio.h>
#include<string.h>
#include<stdlib.h>

struct inhalt {
	FILE *save;
	double *numbers;
	int columns;
	int rows;
};

int print_main(struct inhalt *data, char *name, char *string);
void print_header(struct inhalt *data, char *string);
void print_number(struct inhalt *data);
void open(struct inhalt *data, char *name);
void close(struct inhalt *data);
int sizeofinhalt(void);
struct inhalt *init_inhalt(void);
void close_inhalt(struct inhalt *data);
void ins_numbers(struct inhalt *data, double *numbers, int columns, int rows);
void copy_numbers(struct inhalt *data, double *numbers, int columns, int rows);
struct inhalt *init_plotter(char *name, char *headerstring);
void cycle_plotter(struct inhalt *data, double *xnumber, double *ynumber);
void close_plotter(struct inhalt *data);
