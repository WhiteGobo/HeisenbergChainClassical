#include<stdio.h>
#include<string.h>
#include<stdlib.h>



/* Read all the random number out of a random file
 * free pointer afterwards
 */
void *intsofsize(char *file, int size, int menge, int *return_list)
{
	int i, return_value, r; r=0;
	char zeichen;
	char number[size];
	FILE *save = fopen(file, "r");
	if(save != NULL){
		zeichen = (char)fgetc(save);
		while( zeichen != EOF && r < menge){
			for(i=0;i<size;i++){
				number[i]=zeichen;
				zeichen = (char)fgetc(save);
			}
			sscanf(number,"%d", &return_value);
			*(return_list + r++) = return_value;
		}
	}
	else printf("ERROR: no file found, open %s", file);
	fclose(save);
}

