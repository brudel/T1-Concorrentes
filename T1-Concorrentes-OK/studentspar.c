// COMO COMPILAR: gcc studentspar.c -o par -lm -fopenmp
// COMO EXECUTAR: ./par < input.in

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define MAX_VAL 101
#define NTHREADS 4

typedef int estogram[MAX_VAL];

typedef struct data {
	double med, dp, median;
} data;

void calcula_metricas(estogram cid, estogram reg, data* d, int n) {
	long int sum = 0, i, sumquad = 0, num = 0;

    // Percorre histograma somando número de ocorrências de cada nota
    // cálcula as somas da média e do dp multiplicando as ocorrências pelas notas
    // constrói também o histograma da região correspondente
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

    // Trata o caso de número par de notas
	if (n % 2 == 0 && num == n / 2) {
		int j = d->median+1;
		while (cid[j] == 0){
			j++;
		}
		d->median = (j + d->median) / 2;
	}

    // Cálcula média e dp
	d->med = (double) sum / n;
	d->dp = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
}

void calcula_metricas_pais(estogram est, data* d, int n) {
	long int sum = 0, i, sumquad = 0, num = 0;

    // Percorre histograma somando número de ocorrências de cada nota
    // cálcula as somas da média e do dp multiplicando as ocorrências pelas notas
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

    // Trata o caso de número par de notas
	if (n % 2 == 0 && num == n / 2) {
		int j = d->median+1;
		while (est[j] == 0){
			j++;
		}
		d->median = (j + d->median) / 2;
	}

    // Calcula média e dp
	d->med = (double) sum / n;
	d->dp = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
}

// Encontra nota máxima
int emax(estogram est) {
	int i;

	for(i = MAX_VAL - 1; est[i] == 0; --i);

	return i;
}

// Encontra nota mínima
int emin(estogram est) {
	int i;

	for(i = 0; est[i] == 0; ++i);

	return i;
}

int main(void)
{
    // Lẽ entrada
	int lines[4], i = 0, j;

	while(scanf(" %d", &lines[i]) != EOF)
		i++;

	int R = lines[0];
	int C = lines[1];
	int A = lines[2];
	int seed = lines[3];

	srand(seed);

	// Gera notas
	double tempoExec;

	int *matrizNotas = (int *)malloc(R*C*A*sizeof(int));
	estogram *cidades = calloc(1, R * C * sizeof(estogram)),
	*regioes = calloc(1, R * sizeof(estogram)),
	*Brasil = calloc(1, sizeof(estogram));

	for(i=0; i<R; i++){
		for(j=0; j<C; j++){
			for(int k=0; k<A; k++){
				matrizNotas[i*C*A + j*A + k] = rand()%MAX_VAL;
			}
		}
	}

	tempoExec = omp_get_wtime();

	// Constroi histograma de notas de cada cidade
	if(R*C >= NTHREADS){
		#pragma omp parallel num_threads(NTHREADS)
		{
			#pragma omp for private(i, j)
			for (i = 0; i < R*C; i++){
				for (j = 0; j < A; j++){
					cidades[i][matrizNotas[i*A + j]]++;
				}
			}
		}
	} else {
		estogram *cidades_loc = (estogram *)calloc(NTHREADS, sizeof(estogram));	
		for (i = 0; i < R*C; i++){
			#pragma omp parallel num_threads(NTHREADS) firstprivate(i) shared(cidades_loc)
			{				
				int my_rank = omp_get_thread_num();

				#pragma omp for private(j)
				for (j = 0; j < A; j++){
					cidades_loc[my_rank][matrizNotas[i*A + j]]++;
				}

			}

			for(int k=0; k<NTHREADS; k++){
				for(int l=0; l<MAX_VAL; l++){
					cidades[i][l] += cidades_loc[k][l];
					cidades_loc[k][l] = 0;
				}
			}
		}
	}

	// Struct de metricas (media, dp e mediana)
	data *dados_cidades = (data *) malloc(R * C * sizeof(data));
	data *dados_regiao = (data *) malloc(R * sizeof(data));
	data dados_brasil;

	// Vetores com maiores e menores notas
	int *maior_cidades = (int *) malloc(R * C * sizeof(int));
	int *menor_cidades = (int *) malloc(R * C * sizeof(int));
	int *maior_regiao = (int *) malloc(R * sizeof(int));
	int *menor_regiao = (int *) malloc(R * sizeof(int));
	int maior_brasil;
	int menor_brasil;
	
	// Cidade/Regiao com melhor media
	int melhor_cidade, melhor_cidade_reg, melhor_regiao;
	double maiorMediaCidade = -1, maiorMediaRegiao = -1;

	// Calcula metricas de cada Região/Cidade e país
	if(R >= NTHREADS){
		#pragma omp parallel num_threads(NTHREADS) shared(cidades, regioes, Brasil, dados_cidades, dados_regiao, maiorMediaCidade, maiorMediaRegiao, melhor_cidade, melhor_cidade_reg, melhor_regiao)
		{       
			double maiorMediaCidadeLoc = -1, maiorMediaRegiaoLoc = 1;
			int melhor_cidade_loc = 0, melhor_cidade_reg_loc = 0, melhor_regiao_loc = 0;
			
			// Paraleliza o calculo das regiões
			#pragma omp for private(i, j)
			for(i=0; i<R; i++) {
				for(j=0; j<C; j++) {
					// Calcula media, dp e mediana da cidade, montando o histograma da regiao também
					calcula_metricas(cidades[i*C + j], regioes[i], &dados_cidades[i*C + j], A);

					// Encontra maior/menor nota
					maior_cidades[i*C + j] = emax(cidades[i*C + j]);
					menor_cidades[i*C + j] = emin(cidades[i*C + j]);
					
					// Atualiza melhor cidade e melhor média locais
					if(dados_cidades[i*C + j].med > maiorMediaCidadeLoc){
						maiorMediaCidadeLoc = dados_cidades[i*C + j].med;
						melhor_cidade_loc = j;
						melhor_cidade_reg_loc = i;
					}
				}

				// Calcula media, dp e mediana da região, montando o histograma do país também
				calcula_metricas(regioes[i], *Brasil, dados_regiao + i, C * A);

				// Encontra maior/menor nota
				maior_regiao[i] = emax(regioes[i]);
				menor_regiao[i] = emin(regioes[i]);

				// Atualiza melhor região e melhor média locais
				if(dados_regiao[i].med > maiorMediaRegiaoLoc){
					maiorMediaRegiaoLoc = dados_regiao[i].med;
					melhor_regiao_loc = i;
				}
			}

			// Atualiza melhor região e melhor média globais
			#pragma omp critical (CHECK_REGIAO)
			{
				if(maiorMediaRegiaoLoc > maiorMediaRegiao){
					maiorMediaRegiao = maiorMediaRegiaoLoc;
					melhor_regiao = melhor_regiao_loc;
				}
			}

			// Atualiza melhor cidade e melhor média globais
			#pragma omp critical (CHECK_CIDADE)
			{
				if(maiorMediaCidadeLoc > maiorMediaCidade){
					maiorMediaCidade = maiorMediaCidadeLoc;
					melhor_cidade = melhor_cidade_loc;
					melhor_cidade_reg = melhor_cidade_reg_loc;
				}
			}
		}
	} else {
		// Histogramas locais de uma região para cada thread
		estogram *regioes_loc = (estogram *)calloc(NTHREADS, sizeof(estogram));

		for(i=0; i<R; i++) {
			// Paraleliza o cálculo das cidades
			#pragma omp parallel num_threads(NTHREADS) shared(cidades, regioes, Brasil, dados_cidades, dados_regiao, maiorMediaCidade, maiorMediaRegiao, melhor_cidade, melhor_cidade_reg, melhor_regiao)
			{
				double maiorMediaCidadeLoc = -1;
				int melhor_cidade_loc = 0, melhor_cidade_reg_loc = 0;

				int my_rank = omp_get_thread_num();

				#pragma omp for firstprivate(i) private(j)
				for(j=0; j<C; j++) {
					// Calcula media, dp e mediana da cidade, montando o histograma local da regiao também
					calcula_metricas(cidades[i*C + j], regioes_loc[my_rank], &dados_cidades[i*C + j], A);

					// Encontra maior/menor nota
					maior_cidades[i*C + j] = emax(cidades[i*C + j]);
					menor_cidades[i*C + j] = emin(cidades[i*C + j]);
					
					// Atualiza melhor cidade e melhor média locais
					if(dados_cidades[i*C + j].med > maiorMediaCidadeLoc){
						maiorMediaCidadeLoc = dados_cidades[i*C + j].med;
						melhor_cidade_loc = j;
						melhor_cidade_reg_loc = i;
					}
				}

				// Atualiza melhor cidade e melhor média globais
				#pragma omp critical (CHECK_CIDADE)
				{
					if(maiorMediaCidadeLoc > maiorMediaCidade){
						maiorMediaCidade = maiorMediaCidadeLoc;
						melhor_cidade = melhor_cidade_loc;
						melhor_cidade_reg = melhor_cidade_reg_loc;
					}
				}
			}

			// Monta o histograma global da região em questão
			for(int k=0; k<NTHREADS; k++){
				for(int l=0; l<MAX_VAL; l++){
					regioes[i][l] += regioes_loc[k][l];
					regioes_loc[k][l] = 0;
				}
			}

			// Calcula media, dp e mediana da região, montando o histograma do país também
			calcula_metricas(regioes[i], *Brasil, dados_regiao + i, C * A);
			maior_regiao[i] = emax(regioes[i]);
			menor_regiao[i] = emin(regioes[i]);

			// Atualiza melhor região e melhor média
			if(dados_regiao[i].med > maiorMediaRegiao){
				maiorMediaRegiao = dados_regiao[i].med;
				melhor_regiao = i;
			}
		}
	}	

	// Calcula media, dp e mediana do país
	calcula_metricas_pais(*Brasil, &dados_brasil, R * C * A);

	// Encontra maior/menor nota
	maior_brasil = emax(*Brasil);
	menor_brasil = emin(*Brasil);

    tempoExec = omp_get_wtime() - tempoExec;

	/// Imprime resultados

	// Métricas das Cidades
	int k;
	for (i = 0; i < R; i++)
	{
		for(j = 0; j < C; j++)
		{
			k = i*C+j;
			printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], dados_cidades[k].median, dados_cidades[k].med, dados_cidades[k].dp);
		}
		printf("\n");
	}

	// Métricas das regiões
	for (i = 0; i < R; i++)
	{
		printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], dados_regiao[i].median, dados_regiao[i].med, dados_regiao[i].dp);
	}

	printf("\n");

	// Métricas do país
	printf("Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, dados_brasil.median, dados_brasil.med, dados_brasil.dp);

	printf("\n");

	printf("Melhor regiao: Regiao %d\n", melhor_regiao);
	printf("Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade_reg, melhor_cidade);

	printf("\nTempo de resposta sem considerar E/S, em segundos: %.3lfs\n", tempoExec);

	return(0);
}