// COMO COMPILAR: gcc studentsseq.c -o seq -lm -fopenmp
// COMO EXECUTAR: ./seq < input.in

// INTEGRANTES:
// Marcelo Kiochi Hatanaka (10295645)
// Paulo Renato Campos Barbosa (9779475)
// Rafael Farias Roque (10295412)
// Rodrigo Mendes Andrade (10262721)
// Bruno Del Monde (10262818)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi/mpi.h>

#define MAX_VAL 101

enum {MEDIA, DESVIO_PADRAO, MEDIANA};
enum {CIDADE, REGIAO};

typedef int estogram[MAX_VAL];

typedef double data[3]; //# Implementar

#define oi printf("oi\n"); fflush(stdout);

void calcula_metricas(estogram cid, estogram reg, double* d, int n) {
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
	d[MEDIANA] = i-1;
	for (; i < MAX_VAL; ++i) {
		sum += cid[i] * i;
		sumquad += cid[i] * i * i;
		reg[i] += cid[i];
	}

	// Trata o caso de número par de notas
	if (n % 2 == 0 && num == n / 2) {
		i = d[MEDIANA]+1;
		while (cid[i] == 0){
			i++;
		}
		d[MEDIANA] = (i + d[MEDIANA]) / 2;
	}

	// Cálcula média e dp
	d[MEDIA] = (double) sum / n;
	d[DESVIO_PADRAO] = sqrt((sumquad  -  (double) sum * sum / n)	/	(n - 1));
}

void calcula_metricas_pais(estogram est, double* d, int n) {
	long int sum = 0, i, sumquad = 0, num = 0;

	// Percorre histograma somando número de ocorrências de cada nota
	// Cálcula as somas da média e do dp multiplicando as ocorrências pelas notas
	for (i = 0; num < (n + 1) / 2; ++i) {
		num += est[i];
		sum += est[i] * i;
		sumquad += est[i] * i * i;
	}
	d[MEDIANA] = i-1;
	for (; i < MAX_VAL; ++i) {
		sum += est[i] * i;
		sumquad += est[i] * i * i;
	}

	// Trata o caso de número par de notas
	if (n % 2 == 0 && num == n / 2) {
		i = d[MEDIANA]+1;
		while (est[i] == 0){
			i++;
		}
		d[MEDIANA] = (i + d[MEDIANA]) / 2;
	}

	// Calcula média e dp
	d[MEDIA] = (double) sum / n;
	d[DESVIO_PADRAO] = sqrt((sumquad  -  (double) sum * sum / n)	/	(n - 1));
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

int main(int argc, char **argv)
{
	// Informações do MPI
	MPI_Init(&argc, &argv);
	int npes, myrank, nRegProc;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// Variáveis do root
	int allR;

	// Variáveis de execução
	int R, C, A;
	int *matrizNotas;
	estogram *cidades, *regioes, *Brasil;
    estogram auxBrasil;

	// Resultados
	double *dados_cidades, *dados_regiao, *dados_brasil;
	int *maior_cidades, *menor_cidades, *maior_regiao, *menor_regiao;
	int menor_brasil, maior_brasil;    

	int melhor_cidade[2], melhor_regiao;
	double maiorMediaCidade = -1, maiorMediaRegiao = -1;
    // Auxiliares
    int melhores_cidades[npes * 2], melhores_regioes[npes];
    double maioresMediasCidade[npes], maioresMediasRegiao[npes];

	// Outros
	int i, j, k;
	double tempoExec;

	FILE *arq;
	char name[9] = "arq0.txt";
	name[3] += myrank;

	arq = fopen(name, "w+");

	printf("meu rank = %d\n", myrank);

	if(myrank == 0) {
		int seed, sizes[npes][3]/*# Tirar redundância #*/, resto, regAtual;

		scanf(" %d %d %d %d", &allR, &C, &A, &seed);
		srand(seed);

		matrizNotas = malloc(allR * C * A * sizeof(int));

		// Gera matriz de notas
		for(i=0; i<allR; i++){
			for(j=0; j<C; j++){
				for(k=0; k<A; k++){
					matrizNotas[i*C*A + j*A + k] = rand()%MAX_VAL;
					fprintf(arq, "%d ", matrizNotas[i*C*A + j*A + k]);
				}
				fprintf(arq, "\n");
			}
			fprintf(arq, "\n");
		}

		// Inicia contagem do tempo
		tempoExec = MPI_Wtime();

		// Metricas (media, dp e mediana)
		dados_cidades = malloc(allR * C * 3 * sizeof(double));
		dados_regiao = malloc(allR * 3 * sizeof(double));
		dados_brasil = malloc(3 * sizeof(double));

		// Vetores com maiores e menores notas
		maior_cidades = malloc(allR * C * sizeof(int));
		menor_cidades = malloc(allR * C * sizeof(int));
		maior_regiao = malloc(allR * sizeof(int));
		menor_regiao = malloc(allR * sizeof(int));

		nRegProc = allR/npes;
		resto = allR%npes;
		R = nRegProc + (resto > 0 ? 1 : 0);
		MPI_Request request[npes];

		regAtual = R;

		// Distribui trabalho
		for(i=1; i<npes; i++){
			sizes[i][1] = C;
			sizes[i][2] = A;
			if(i < resto){
				sizes[i][0] = nRegProc+1; 
			} else {
				sizes[i][0] = nRegProc;
			}
			
			MPI_Isend(sizes[i], 3, MPI_INT, i, 0, MPI_COMM_WORLD, &request[i]); //# Pode ser transformado em broadcast e cálculo local
			MPI_Isend(&matrizNotas[regAtual*C*A], sizes[i][0]*C*A, MPI_INT, i, 1, MPI_COMM_WORLD, &request[i]); //# Tenho quase certeza que pode ser feito com scatter

			regAtual += sizes[i][0];
		}

	} else {
		int sizes[3];
		MPI_Recv(sizes, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);

		R = sizes[0], C = sizes[1], A = sizes[2];
		matrizNotas = malloc(R*C*A*sizeof(int));

		// Metricas (media, dp e mediana)
		dados_cidades = malloc(R * C * 3 * sizeof(double));
		dados_regiao = malloc(R * 3 * sizeof(double));
		dados_brasil = malloc(3 * sizeof(double));

		// Vetores com maiores e menores notas
		maior_cidades = malloc(R * C * sizeof(int));
		menor_cidades = malloc(R * C * sizeof(int));
		maior_regiao = malloc(R * sizeof(int));
		menor_regiao = malloc(R * sizeof(int));

		MPI_Recv(matrizNotas, R*C*A, MPI_INT, 0, 1, MPI_COMM_WORLD, NULL);

		fprintf(arq, "Rank = %d\n", myrank);
		for(i=0; i<R; i++)
			for(j=0; j<C; j++)
				for(k=0; k<A; k++)
					fprintf(arq, "%d ", matrizNotas[i*C*A + j*A + k]);
	}

	cidades = calloc(1, R * C * sizeof(estogram)),
	regioes = calloc(1, R * sizeof(estogram)),
	Brasil = calloc(1, sizeof(estogram));

	for(i=0; i<R*C; i++){
		for(j=0; j<A; j++){
			cidades[i][matrizNotas[i*A + j]]++;
		}
	}

	// Calcula metricas de cada Região/Cidade e país
	for(i=0; i<R; i++) {
		for(j=0; j<C; j++) {
			// Calcula media, dp e mediana da cidade, montando o histograma da regiao também
			calcula_metricas(cidades[i*C + j], regioes[i], &dados_cidades[i*C*3 + j*3], A);

			// Encontra maior/menor nota
			maior_cidades[i*C + j] = emax(cidades[i*C + j]);
			menor_cidades[i*C + j] = emin(cidades[i*C + j]);
			
			// Atualiza melhor cidade e melhor média
			if(dados_cidades[i*C*3 + j*3 + MEDIA] > maiorMediaCidade){
				maiorMediaCidade = dados_cidades[i*C*3 + j*3 + MEDIA];
				melhor_cidade[CIDADE] = j;
				melhor_cidade[REGIAO] = i;
			}
		}

		// Calcula media, dp e mediana da região, montando o histograma do país também
		calcula_metricas(regioes[i], *Brasil, dados_regiao + i * 3, C * A);

		// Encontra maior/menor nota
		maior_regiao[i] = emax(regioes[i]);
		menor_regiao[i] = emin(regioes[i]);

		// Atualiza melhor região e melhor média
		if(dados_regiao[i*3 + MEDIA] > maiorMediaRegiao){
			maiorMediaRegiao = dados_regiao[i*3 + MEDIA];
			melhor_regiao = i;
		}
	}

	if(myrank == 0) {
		int nRegs;
		int regAtual = R;
		int resto = allR%npes;
		MPI_Request requests[npes][6];
		
		for(i=1; i<npes; i++){
			if(i < resto){
				nRegs = nRegProc+1; 
			} else {
				nRegs = nRegProc;
			}
			
			MPI_Irecv(&dados_cidades[regAtual*C*3], nRegs*C*3, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &requests[i][0]);
			MPI_Irecv(&dados_regiao[regAtual*3], nRegs*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &requests[i][1]);

			MPI_Irecv(&menor_cidades[regAtual*C], nRegs*C, MPI_INT, i, 2, MPI_COMM_WORLD, &requests[i][2]);
			MPI_Irecv(&maior_cidades[regAtual*C], nRegs*C, MPI_INT, i, 3, MPI_COMM_WORLD, &requests[i][3]);
			MPI_Irecv(&menor_regiao[regAtual], nRegs, MPI_INT,i, 4, MPI_COMM_WORLD, &requests[i][4]);
			MPI_Irecv(&maior_regiao[regAtual], nRegs, MPI_INT, i, 5, MPI_COMM_WORLD, &requests[i][5]);

			regAtual += nRegs;
		}

        MPI_Status status[npes][6];
		for(i=1; i<npes; i++){
			for(j=0; j<6; j++){
				MPI_Wait(&(requests[i][j]), &(status[i][j]));
			}
		}

	} else {
		int auxBrasil[MAX_VAL];

		MPI_Request requests[6];
		MPI_Isend(dados_cidades, R * C * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests[0]);
		MPI_Isend(dados_regiao, R * 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &requests[1]);

		MPI_Isend(menor_cidades, R*C, MPI_INT, 0, 2, MPI_COMM_WORLD, &requests[2]);
		MPI_Isend(maior_cidades, R*C, MPI_INT, 0, 3, MPI_COMM_WORLD, &requests[3]);
		MPI_Isend(menor_regiao, R, MPI_INT, 0, 4, MPI_COMM_WORLD, &requests[4]);
		MPI_Isend(maior_regiao, R, MPI_INT, 0, 5, MPI_COMM_WORLD, &requests[5]);
	}

    MPI_Reduce(*Brasil, auxBrasil, MAX_VAL, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Gather(melhor_cidade, 2, MPI_INT, melhores_cidades, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&melhor_regiao, 1, MPI_INT, melhores_regioes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&maiorMediaCidade, 1, MPI_DOUBLE, maioresMediasCidade, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&maiorMediaRegiao, 1, MPI_DOUBLE, maioresMediasRegiao, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (i = 0; i < R; i++)
	{
		for(j = 0; j < C; j++)
		{
			k = i*C+j;
			fprintf(arq, "Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], dados_cidades[k*3 + MEDIANA], dados_cidades[k*3 + MEDIA], dados_cidades[k*3 + DESVIO_PADRAO]);
		}
		fprintf(arq, "\n");
	}

	// Métricas das regiões
	for (i = 0; i < R; i++)
	{
		fprintf(arq, "Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], dados_regiao[i*3 + MEDIANA], dados_regiao[i*3 + MEDIA], dados_regiao[i*3 + DESVIO_PADRAO]);
	}

	fprintf(arq, "\n");

	if(myrank == 0) {
        int nRegs;
		int regAtual = R;
		int resto = allR%npes;
        for(int i=0; i<MAX_VAL; i++){
            (*Brasil)[i] = auxBrasil[i];
        }

		calcula_metricas_pais(*Brasil, dados_brasil, allR*C*A);

		//Encontra maior/menor nota
		maior_brasil = emax(*Brasil);
		menor_brasil = emin(*Brasil);

        regAtual = 0;
        for(int i=0; i<npes; i++){
            
            if(maioresMediasCidade[i] > maiorMediaCidade){
                maiorMediaCidade = maioresMediasCidade[i];
                melhor_cidade[CIDADE] = melhores_cidades[i*2 + CIDADE];
                melhor_cidade[REGIAO] = regAtual + melhores_cidades[i*2 + REGIAO];
            }

            if(maioresMediasRegiao[i] > maiorMediaRegiao){
                maiorMediaRegiao = maioresMediasRegiao[i];
                melhor_regiao = regAtual + melhores_regioes[i];
            }

            regAtual += nRegProc + (i < resto ? 1 : 0);
        }

		// Termina contagem do tempo
		tempoExec = MPI_Wtime() - tempoExec;
        
		for (i = 0; i < allR; i++)
		{
			for(j = 0; j < C; j++)
			{
				k = i*C+j;
				printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], dados_cidades[k*3 + MEDIANA], dados_cidades[k*3 + MEDIA], dados_cidades[k*3 + DESVIO_PADRAO]);
			}
			printf("\n");
		}

		// Métricas das regiões
		for (i = 0; i < allR; i++)
		{
			printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], dados_regiao[i*3 + MEDIANA], dados_regiao[i*3 + MEDIA], dados_regiao[i*3 + DESVIO_PADRAO]);
		}

		printf("\n");

		//Métricas do país
		printf("Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, dados_brasil[MEDIANA], dados_brasil[MEDIA], dados_brasil[DESVIO_PADRAO]);

		printf("\n");

		printf("Melhor regiao: Regiao %d\n", melhor_regiao);
		printf("Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade[REGIAO], melhor_cidade[CIDADE]);

		printf("\nTempo de resposta sem considerar E/S, em segundos: %.3lfs\n", tempoExec);
	}	

	free(dados_cidades);
	free(dados_regiao);
	free(dados_brasil);
	free(maior_cidades);
	free(menor_cidades);
	free(maior_regiao);
	free(menor_regiao);
	free(cidades);
	free(regioes);
	free(Brasil);
	free(matrizNotas);
	fclose(arq);
	MPI_Finalize();

	return(0);
}