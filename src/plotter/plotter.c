#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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


/* This here creates the struct inhalt
 */
struct inhalt *init_inhalt(void)
{
	return (struct inhalt *)malloc(sizeofinhalt());
}



void close_inhalt(struct inhalt *data)
{
	free(data);
}


/* use this to insert the rows and columns
 */
void ins_numbers(struct inhalt *data, double *numbers, int columns, int rows)
{
	data->numbers = numbers;
	data->columns = columns;
	data->rows = rows;
}


/* use this to copy numbers to a new location
 * TODO: write this method make clear that within malloc one free the 
 * 	space at end of the program
 */
void copy_numbers(struct inhalt *data, double *numbers, int columns, int rows)
{

}

/* Call this Method to print a file. It should work automatically
 * This method handles Errors
 * Returns Errors if exists
 */
int print_main(struct inhalt *data, char *name, char *string)
{
	open(data, name);
	print_header(data, string);
	print_number(data);
	close(data);
	return 0;
}


/* This Method is called to print a header into a file
 * optimized for heisenbergchain need to implement spinchain 
 */
void print_header(struct inhalt *data, char *string)
{
	fprintf(data->save, string);
}


/* This Method ist here to write the data into the file
 * TODO: divide this method for rows and matrices (This is for Matrix)
 */
void print_number(struct inhalt *data)
{
	int i, j;
	double *tip = data->numbers;
	for (i=0; i<data->rows; i++)
	{
		for (j=0; j<data->columns; j++)
		{
			fprintf(data->save, "%E\t", *tip++);
		}
		fprintf(data->save, "\n");
	}
}


/* Open a file for the datastream
 */
void open(struct inhalt *data, char *name)
{
	FILE *save = fopen(name, "w");
	data->save = save;	
}


/*close the file opened by open()
 */
void close(struct inhalt *data)
{
	fclose(data->save);
}


/* Returns size of struct inhalt
 * TODO: remove this function
 */
int sizeofinhalt(void)
{
	return sizeof(struct inhalt);
}


/* begins plotting a file for continuous plotting 
 * insert header
 */
struct inhalt *init_plotter(char *name, char *headerstring)
{
	struct inhalt *data = (struct inhalt *)malloc(sizeofinhalt());

	open(data, name);
	print_header(data, headerstring);
	return data;
}


void cycle_plotter(struct inhalt *data, double *xnumber, double *ynumber)
{
	fprintf(data->save, "%f\t", *xnumber);
	fprintf(data->save, "%f\n", *ynumber);
}


void close_plotter(struct inhalt *data)
{
	close(data);
}
