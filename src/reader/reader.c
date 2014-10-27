#include<stdio.h>
#include<string.h>
#include<stdlib.h>



/* Read all the random number out of a random file
 * free pointer afterwards
 * ERROR: 0 Success; -1 File not big enough; -2 no File exists
 */
int intsofsize(char *file, int size, int menge, int *return_list)
{
	char string[79];
	int i, return_value, r; r=0;
	char zeichen;
	char number[size];
	FILE *save = fopen(file, "r");
	if(save != NULL){
		while(r < menge){
			for(i=0;i<size;i++){
				zeichen = (char)fgetc(save);
				if (zeichen == EOF){
					sprintf(string,"File %s hasnt enough numbers(needed: %d)\n", file, menge);
					fputs(string, stderr);
					fclose(save);
					return -1;
				}
				number[i]=zeichen;
			}
			sscanf(number,"%d", &return_value);
			*(return_list + r++) = return_value;
		}
	}
	else{
		sprintf(string, "File %s couldnt be opened.\n", file);
		fputs(string, stderr);
		return -2;
	}
	fclose(save);
	return 0;
}

