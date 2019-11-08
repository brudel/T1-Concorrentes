// COMO COMPILAR: gcc alternativapar.c -o main -lm -fopenmp
// COMO EXECUTAR: ./main < entrada.txt

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define swap(i, j) aux = i, i = j, j = aux;
#define MAX_VAL 101
#define NTHREADS 4

typedef int estogram[MAX_VAL];

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
	d->median = i-1;
	for (; i < MAX_VAL; ++i) {
		sum += cid[i] * i;
		sumquad += cid[i] * i * i;
		reg[i] += cid[i];
	}

	if (n % 2 == 0 && num == n / 2) {
		int j = d->median+1;
		while (cid[j] == 0){
			j++;
		}
		d->median = (j + d->median) / 2;
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

	d->median = i-1;
	for (; i < MAX_VAL; ++i) {
		sum += est[i] * i;
		sumquad += est[i] * i * i;
	}

	if (n % 2 == 0 && num == n / 2) {
		int j = d->median+1;
		while (est[j] == 0){
			j++;
		}
		d->median = (j + d->median) / 2;
	}

	d->med = (double) sum / n;
	d->dp = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
}

int emax(estogram est) {
	int i;
	for(i = MAX_VAL - 1; est[i] == 0; --i);
	return i;
}

int emin(estogram est) {
	int i;
	for(i = 0; est[i] == 0; ++i);
	return i;
}

void sum_est(estogram A, estogram B) {
	for (int i = 0; i < MAX_VAL; ++i)
		A[i] += B[i];
}

int main(void)
{
	int lines[4], i = 0, j;
	double tempoExec;

	while(scanf(" %d", &lines[i]) != EOF)
		i++;

	int R = lines[0];
	int C = lines[1];
	int A = lines[2];
	int seed = lines[3];

	/// Generates the arrays
	srand(seed);
	int *matrizNotas = (int *)malloc(R*C*A*sizeof(int));
	estogram *matriz = calloc(1, R * C * sizeof(estogram)),
	*regioes = calloc(1, R * sizeof(estogram)),
	*Brasil = calloc(1, sizeof(estogram));

	// Gera notas
	for(i=0; i<R; i++){
		for(j=0; j<C; j++){
			for(int k=0; k<A; k++){
				matrizNotas[i*C*A + j*A + k] = rand()%MAX_VAL;
			}
		}
	}

	tempoExec = omp_get_wtime();

	// Controi histograma de notas de cada cidade
	for (i = 0; i < R*C; i++){
		for (j = 0; j < A; j++){
			matriz[i][matrizNotas[i*A + j]]++;
		}
	}

	data *data_cidades = (data *) malloc(R * C * sizeof(data));
	data *data_regiao = (data *) malloc(R * sizeof(data));
	data data_brasil;

	int *maior_cidades = (int *) malloc(R * C * sizeof(int));
	int *maior_regiao = (int *) malloc(R * sizeof(int));
	int maior_brasil;

	int *menor_cidades = (int *) malloc(R * C * sizeof(int));
	int *menor_regiao = (int *) malloc(R * sizeof(int));
	int menor_brasil ;
	
	int melhor_cidade, melhor_cidade_reg, melhor_regiao;
	double maiorMediaCidade = -1, maiorMediaRegiao = -1;

	for(i=0; i<R; i++) {
		for(j=0; j<C; j++) {
			analyzepart(matriz[i*C + j], regioes[i], &data_cidades[i*C + j], A);
			maior_cidades[i*C + j] = emax(matriz[i*C + j]);
			menor_cidades[i*C + j] = emin(matriz[i*C + j]);
			
			if(data_cidades[i*C + j].med > maiorMediaCidade){
				maiorMediaCidade = data_cidades[i*C + j].med;
				melhor_cidade = j;
				melhor_cidade_reg = i;
			}
		}
		analyzepart(regioes[i], *Brasil, data_regiao + i, C * A);
		maior_regiao[i] = emax(regioes[i]);
		menor_regiao[i] = emin(regioes[i]);

		if(data_regiao[i].med > maiorMediaRegiao){
			maiorMediaRegiao = data_regiao[i].med;
			melhor_regiao = i;
		}
	}

	// double tt = omp_get_wtick();
	int tn = omp_get_thread_num();
	printf("thread %d i=%d\n", tn, i);

	justanalyze(*Brasil, &data_brasil, R * C * A);
	maior_brasil = emax(*Brasil);
	menor_brasil = emin(*Brasil);

    tempoExec = omp_get_wtime() - tempoExec;

	/// Printing results

	// Metrics for cities
	int k;
	// for (i = 0; i < R; i++)
	// {
	// 	for(j = 0; j < C; j++)
	// 	{
	// 		k = i*C+j;
	// 		printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], data_cidades[k].median, data_cidades[k].med, data_cidades[k].dp);
	// 	}
	// 	printf("\n");
	// }

	// // Metrics for regions
	// for (i = 0; i < R; i++)
	// {
	// 	printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], data_regiao[i].median, data_regiao[i].med, data_regiao[i].dp);
	// }

	// printf("\n");

	// Metrics for the country
	printf("Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, data_brasil.median, data_brasil.med, data_brasil.dp);

	printf("\n");

	printf("Melhor regiao: Regiao %d\n", melhor_regiao);
	printf("Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade_reg, melhor_cidade);

	printf("\nTempo de resposta sem considerar E/S, em segundos: %.3lfs\n", tempoExec);

	return(0);

}
