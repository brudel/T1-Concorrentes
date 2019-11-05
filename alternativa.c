#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include <time.h>

#define swap(i, j) aux = i, i = j, j = aux;

typedef int estogram[100];

int A;

float mediana(estogram est) {
	int sum = -1, i;

	for (i = 0; sum < A / 2; ++i)
		sum += est[i];

	return i;
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

	/// Generates the matrix
	srand(seed);
	estogram *matriz = calloc(1, R * C * sizeof(estogram)),
	*regioes = calloc(1, R * sizeof(estogram)),
	*Brasil = calloc(1, sizeof(estogram));
	for (i = 0; i < R*C; i++)
		for (j = 0; j < A; j++)
			++matriz[i][rand() % 100];



	// Print array in matrix's representation	
	// for (i = 0; i < R*C; i++)
	// {
	//	 for (j = 0; j < A; j++)
	//		 printf("%d ", matriz[i*A+j]);
	//	 printf("\n");
	// }	
   
	/// Calculate the metrics
	/// Mean and deviation
	float *media_cidades = (float *) calloc(1, R*C*sizeof(float));
	float *media_regiao = (float *) calloc(1, R*sizeof(float));
	float media_brasil;

	float *dp_cidades = (float *) calloc(1, R*C*sizeof(float));
	float *dp_regiao = (float *) calloc(1, R*sizeof(float));
	float dp_brasil;

	/// Median, max and min
	float *mediana_cidades = (float *) calloc(1, R*C*sizeof(float));
	float *mediana_regiao = (float *) calloc(1, R*sizeof(float));
	float mediana_brasil;

	int *maior_cidades = (int *) calloc(1, R*C*sizeof(int));
	int *maior_regiao = (int *) calloc(1, R*sizeof(int));
	int maior_brasil = -1;

	int *menor_cidades = (int *) calloc(1, R*C*sizeof(int));
	int *menor_regiao = (int *) calloc(1, R*sizeof(int));
	int menor_brasil = INT_MAX;

	// Calcula media, desvio padrão, max e min das cidades, regioes e país, propagando as somas parciais
	//  - Ex: a soma total de cada cidade é usada como soma parcial para a media da região
	// encontra também a melhor cidade e a melhor região
	float maior_cid = -1, maior_reg = -1;
	int melhor_cidade, melhor_cidade_reg, melhor_regiao;

	float s1, s2, s3, s4, s5 = 0, s6 = 0;
	float media;
	int elem, maiorC = -1, maiorR = -1, menorC = INT_MAX, menorR = INT_MAX;
	tempoExec = clock();
	nclocksMedia = clock();
	/*for(int i=0; i<R; i++){
		s3 = 0;
		s4 = 0;
		maiorR = -1;
		menorR = INT_MAX;
		for(int j=0; j<C; j++){
			s1 = 0;
			s2 = 0;
			maiorC = -1;
			menorC = INT_MAX;
			for(int k=0; k<A; k++){
				// calcula soma e soma do quadrado das notas da cidade
				elem = matriz[i*C*A + j*A + k];
				s1 += elem;
				s2 += elem * elem;

				// compara com maior e menor nota da cidade
				if(maiorC < elem){
					maiorC = elem;
				}

				if(menorC > elem){
					menorC = elem;
				}
			}
			// calcula media e dp
			media = s1/A;
			media_cidades[i*C + j] = media;
			dp_cidades[i*C + j] = sqrt((s2  -  (float) s1 * s1 / A)	/	(A - 1));

			// compara com melhor cidade
			if(media > maior_cid){
				maior_cid = media;
				melhor_cidade = j;
				melhor_cidade_reg = i;
			}

			// salva maior nota da cidade
			maior_cidades[i*C + j] = maiorC;
			menor_cidades[i*C + j] = menorC;

			// compara com maior e menor nota da regiao
			if(maiorR < maiorC){
				maiorR = maiorC;
			}

			if(menorR > menorC){
				menorR = menorC;
			}

			// propaga somas para o calculo da regiao
			s3 += s1;
			s4 += s2;
		}

		// calcula media e dp da regiao com as somas das cidades
		media = s3/(C*A);
		media_regiao[i] = media;
		dp_regiao[i] = sqrt((s4  -  (float) s3 * s3 / (C*A))	/	(C*A - 1));

		// compara com melhor regiao
		if(media > maior_reg){
			maior_reg = media;
			melhor_regiao = i;
		}

		// salva maior e menor nota da região
		maior_regiao[i] = maiorR;
		menor_regiao[i] = menorR;

		// compara com maior e menor nota do pais
		if(maior_brasil < maiorR){
			maior_brasil = maiorR;
		}
		
		if(menor_brasil > menorR){
			menor_brasil = menorR;
		}

		// propaga somas para o calculo do país
		s5 += s3;
		s6 += s4;
	}
	// calcula media e dp do país
	media_brasil = s5/(R*C*A);
	dp_brasil = sqrt((s6  -  (float) s5 * s5 / (R*C*A))	/	(R*C*A - 1));
	nclocksMedia = clock() - nclocksMedia; */

	// Mediana da cidade
	nclocksSortingC = clock();
	for(int i=0; i<R; i++)
		for(int j=0; j<C; j++){
			mediana_cidades[i*C + j] = mediana(matriz[i*C + j]);
			sum_est(regioes[i], matriz[i*C + j]);
		}
	
	nclocksSortingC = clock() - nclocksSortingC;

	// printf("\n");
	// for (i = 0; i < R*C; i++)
	// {
	//	 for (j = 0; j < A; j++)
	//		 printf("%d ", matriz[i*A+j]);
	//	 printf("\n");
	// }

	// Mediana da regiao
	nclocksSortingR = clock();
	//ordena_linhas(matriz, R, A*C);
	// encontra_maiores(matriz, maior_regiao, R, A*C);
	// encontra_menores(matriz, menor_regiao, R, A*C);
	//calcula_mediana(matriz, mediana_regiao, R, A*C);


	for(int i=0; i<R; i++)	{
		mediana_regiao[i] = mediana(regioes[i]);
		sum_est(*Brasil, regioes[i]);
	}
	nclocksSortingR = clock() - nclocksSortingR;

	// printf("\n");
	// for (i = 0; i < R*C; i++)
	// {
	//	 for (j = 0; j < A; j++)
	//		 printf("%d ", matriz[i*A+j]);
	//	 printf("\n");
	// }

	// Mediana do Brasil
	nclocksSortingP = clock();
	//ordena_linhas(matriz, 1, A*C*R);
	// encontra_maiores(matriz, &maior_brasil, 1, A*C*R);
	// encontra_menores(matriz, &menor_brasil, 1, A*C*R);
	//calcula_mediana(matriz, &mediana_brasil, 1, A*C*R);
	mediana_brasil = mediana(*Brasil);
	nclocksSortingP = clock() - nclocksSortingP;

	tempoExec = (clock() - tempoExec)/CLOCKS_PER_SEC;
	/*
	// printf("\n");
	// for (i = 0; i < R*C; i++)
	// {
	//	 for (j = 0; j < A; j++)
	//		 printf("%d ", matriz[i*A+j]);
	//	 printf("\n");
	// }

	/// Printing results

	// Metrics for cities
	int k;
	for (i = 0; i < R; i++)
	{
		for(j = 0; j < C; j++)
		{
			k = i*C+j;
			printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], mediana_cidades[k], media_cidades[k], dp_cidades[k]);
		}
		printf("\n");
	}

	// Metrics for regions
	for (i = 0; i < R; i++)
	{
		printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], mediana_regiao[i], media_regiao[i], dp_regiao[i]);
	}

	printf("\n");*/

	// Metrics for the country
	printf("Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, mediana_brasil, media_brasil, dp_brasil);

	printf("\n");

	printf("Melhor regiao: Regiao %d\n", melhor_regiao);
	printf("Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade_reg, melhor_cidade);

	printf("Tempo de resposta sem considerar E/S, em segundos: %.3lfs\n", tempoExec);
	double tMedia, tC, tR, tP;
	tMedia = nclocksMedia/CLOCKS_PER_SEC;
	tC = nclocksSortingC/CLOCKS_PER_SEC;
	tR = nclocksSortingR/CLOCKS_PER_SEC;
	tP = nclocksSortingP/CLOCKS_PER_SEC;
	printf("tMedia=%.3lf\n tCidades=%.3lf\n tRegioes=%.3lf\n tPais=%.3lf\n", tMedia, tC, tR, tP);

	return(0);

}
