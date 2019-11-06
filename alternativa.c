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

void analyzepart(estogram cid, estogram reg, data* d, int n) {
	long int sum = 0, i, sumquad = 0, num = 0;

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

	d->med = (double) sum / n;
	d->dp = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
}

void justanalyze(estogram est, data* d, int n) {
	long int sum = 0, i, sumquad = 0, num = 0;

	for (i = 0; num < (n + 1) / 2; ++i) {
		num += est[i];
		sum += est[i] * i;
		sumquad += est[i] * i * i;
	}

	d->median = --i;
	if (n % 2 == 0 && num == n / 2) {
		while (est[++i] == 0);
		d->median = (i + d->median) / 2;
	}

	for (; i < 100; ++i) {
		sum += est[i] * i;
		sumquad += est[i] * i * i;
	}

	d->med = (double) sum / n;
	d->dp = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
}

int emax(estogram est) {
	int i;
	for(i = 99; est[i] == 0; --i);
	return i;
}

int emin(estogram est) {
	int i;
	for(i = 0; est[i] == 0; ++i);
	return i;
}

void sum_est(estogram A, estogram B) {
	for (int i = 0; i < 100; ++i)
		A[i] += B[i];
}

int main(void)
{
	int lines[4], i = 0, j;
	double tempoExec = clock();

	while(scanf(" %d", &lines[i]) != EOF)
		i++;

	int R = lines[0];
	int C = lines[1];
	int A = lines[2];
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
			analyzepart(matriz[i*C + j], regioes[i], data_cidades + i*C + j, A);
			maior_cidades[i*C + j] = emax(matriz[i*C + j]);
			menor_cidades[i*C + j] = emin(matriz[i*C + j]);
		}
		analyzepart(regioes[i], *Brasil, data_regiao + i, C * A);
		maior_regiao[i] = emax(regioes[i]);
		menor_regiao[i] = emin(regioes[i]);
	}
	justanalyze(*Brasil, &data_brasil, R * C * A);
	maior_brasil = emax(*Brasil);
	menor_brasil = emin(*Brasil);

	tempoExec = (clock() - tempoExec) / CLOCKS_PER_SEC;

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

	return(0);

}
