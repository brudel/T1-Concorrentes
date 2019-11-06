#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define swap(i, j) aux = i, i = j, j = aux;

typedef int estogram[100];

typedef struct data {
	double med, dp, median;
} data;

int A;

int max(int* vector, int n) {
	int m = vector[0];
	for(int i=1;i<n;i++) {
		if(vector[i] > m) {
			m = vector[i];
		}
	}
	return m;
}

int min(int* vector, int n) {
	int m = vector[0];
	for(int i=1;i<n;i++) {
		if(vector[i] < m) {
			m = vector[i];
		}
	}
	return m;
}

int emax(estogram est) {
	int i;
	for(i = 99; est[i] != 0; --i);
	return i;
}

int emin(estogram est) {
	int i;
	for(i = 0; est[i] != 0; ++i);
	return i;
}

void analyzecid(estogram cid, estogram reg, data* d, int n) {
	int sum = 0, i, sumquad = 0, num = 0;

	for (i = 0; num < (n + 1) / 2; ++i) {
		num += cid[i];
		sum += cid[i] * i;
		sumquad += cid[i] * i * i;
		reg[i] += cid[i];
	}

	d->median = i - 1;
	if (n % 2 == 0 && num == n / 2) {
		while (cid[++i] == 0);
		d->median = (i + d->median) / 2;
	}

	for (; i < 100; ++i) {
		sum += cid[i] * i;
		sumquad += cid[i] * i * i;
		reg[i] += cid[i];
	}

	d->med = (float) sum / n;
	d->dp = sqrt((sumquad  -  (float) sum * sum / n)    /    (n - 1));
}

void sum_est(estogram A, estogram B) {
	for (int i = 0; i < 100; ++i)
		A[i] += B[i];
}

int main(void)
{
	/// Read param from file
	int lines[4], i = 0, j;
	double tempoExec;
	double nclocksMedia, nclocksSortingC, nclocksSortingR, nclocksSortingP;

	// FILE *fp = fopen("entrada.txt", "r");
	// if(fp == NULL)
	// {
	//	 printf("Unable to open file!");
	//	 exit(1);
	// }

	while(scanf(" %d", &lines[i]) != EOF)
		i++;

	int R = lines[0];
	int C = lines[1];
	A = lines[2];
	int seed = lines[3];

	/// Generates the arrays
	srand(seed);
	estogram *matriz = calloc(1, R * C * sizeof(estogram)),
	*regioes = calloc(1, R * sizeof(estogram)),
	*Brasil = calloc(1, sizeof(estogram));

	for (i = 0; i < R*C; i++)
		for (j = 0; j < A; j++)
			++matriz[i][rand() % 100];


	data *data_cidades = (data *) malloc(R * C * sizeof(data));
	data *data_regiao = (data *) malloc(R * sizeof(data));
	data data_brasil;

	int *maior_cidades = (int *) malloc(R * C * sizeof(int));
	int *maior_regiao = (int *) malloc(R * sizeof(int));
	int maior_brasil = -1;

	int *menor_cidades = (int *) malloc(R * C * sizeof(int));
	int *menor_regiao = (int *) malloc(R * sizeof(int));
	int menor_brasil = 100;

	for(int i=0; i<R; i++) {
		for(int j=0; j<C; j++) {
			analyzecid(matriz[i*C + j], regioes[i], data_cidades + i*C + j, A);
		}
		//menor_cidades[i*C + j] = analyze(regioes[i]);
		//sum_est(*Brasil, regioes[i]);
	}
	//mediana_brasil = analyze(*Brasil);

	/// Printing results

	// Metrics for cities
	int k;
	for (i = 0; i < R; i++)
	{
		for(j = 0; j < C; j++)
		{
			k = i*C+j;
			printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], data_cidades[k].median, data_cidades[k].med, data_cidades[k].dp);
		}
		printf("\n");
	}

	// Metrics for regions
	for (i = 0; i < R; i++)
	{
		printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], data_regiao[i].median, data_regiao[i].med, data_regiao[i].dp);
	}

	printf("\n");

	// Metrics for the country
	printf("Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, data_brasil.median, data_brasil.med, data_brasil.dp);

	printf("\n");

	//printf("Melhor regiao: Regiao %d\n", melhor_regiao);
	//printf("Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade_reg, melhor_cidade);

	printf("Tempo de resposta sem considerar E/S, em segundos: %.3lfs\n", tempoExec);
	double tMedia, tC, tR, tP;
	tMedia = nclocksMedia/CLOCKS_PER_SEC;
	tC = nclocksSortingC/CLOCKS_PER_SEC;
	tR = nclocksSortingR/CLOCKS_PER_SEC;
	tP = nclocksSortingP/CLOCKS_PER_SEC;
	printf("tMedia=%.3lf\n tCidades=%.3lf\n tRegioes=%.3lf\n tPais=%.3lf\n", tMedia, tC, tR, tP);

	return(0);

}
