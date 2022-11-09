#include<stdio.h>
#include<string.h>
#include<math.h>

float ran1(long* idum);
float gasdev(long* idum);

int main() {
	FILE* fp;
	long idum;
	float n;
	int arr[4] = { 1000,100,10000,100000 };
	char* uniform[4] = { "uniform1000.txt","uniform100.txt", "uniform10000.txt","uniform100000.txt" };
	char* gauss[4] = { "gauss1000.txt","gauss100.txt", "gauss10000.txt","gauss100000.txt" };
	for (int i = 0; i < 4; i++) {
		fp = fopen(uniform[i], "w");
		for (int j = 0; j < arr[i]; j++) {
			n = ran1(&idum) * 5 - 3;
			char buf[10];
			sprintf(buf, "%f", n);
			fputs(buf, fp);
			fputs("\n", fp);
		}
	}
	for (int i = 0; i < 4; i++) {
		fp = fopen(gauss[i], "w");
		for (int j = 0; j < arr[i]; j++) {
			n = gasdev(&idum) * 1.5 - 0.5;
			char buf[10];
			sprintf(buf, "%f", n);
			fputs(buf, fp);
			fputs("\n", fp);
		}
	}
	fclose(fp);
}
