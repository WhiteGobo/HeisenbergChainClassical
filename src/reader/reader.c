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
		zeichen = (char)fgetc(save);
		while( zeichen != EOF && r < menge){
			if (zeichen != EOF){
				sprintf(string,"File %s hasnt enough numbers(needed: %d)", file, size);
				puts(string);
				fclose(save);
				return -1;
			}
			for(i=0;i<size;i++){
				number[i]=zeichen;
				zeichen = (char)fgetc(save);
			}
			sscanf(number,"%d", &return_value);
			*(return_list + r++) = return_value;
		}
	}
	else{
		sprintf(string, "No file %s found", file);
		puts(string);
		return -2;
	}
	fclose(save);
	return 0;
}

