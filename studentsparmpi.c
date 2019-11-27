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

#define MAX_VAL 11

#define MEDIA 0
#define DESVIO_PADRAO 1
#define MEDIANA 2

#define CIDADE 0
#define REGIAO 1

typedef int estogram[MAX_VAL];

typedef struct data {
	double med, dp, median;
} data;

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
		int j = d[MEDIANA]+1;
		while (cid[j] == 0){
			j++;
		}
		d[MEDIANA] = (j + d[MEDIANA]) / 2;
	}

    // Cálcula média e dp
	d[MEDIA] = (double) sum / n;
	d[DESVIO_PADRAO] = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
}

void calcula_metricas_pais(estogram est, double* d, int n) {
	long int sum = 0, i, sumquad = 0, num = 0;

    // Percorre histograma somando número de ocorrências de cada nota
    // cálcula as somas da média e do dp multiplicando as ocorrências pelas notas
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
		int j = d[MEDIANA]+1;
		while (est[j] == 0){
			j++;
		}
		d[MEDIANA] = (j + d[MEDIANA]) / 2;
	}

    // Calcula média e dp
	d[MEDIA] = (double) sum / n;
	d[DESVIO_PADRAO] = sqrt((sumquad  -  (double) sum * sum / n)    /    (n - 1));
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
    MPI_Init(&argc, &argv);

	// Lê entrada
	int lines[4], i = 0, j;
    int npes, myrank, nRegProc;   
	double tempoExec;


    MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    FILE *arq;
    char name[9] = "arq0.txt";
    name[3] += myrank;

    arq = fopen(name, "w+");

    if(myrank == 0){
        while(scanf(" %d", &lines[i]) != EOF)
            i++;

        int R = lines[0];
        int C = lines[1];
        int A = lines[2];
        int seed = lines[3];
        int sizes[npes][3]; //# Tirar redundância 

        srand(seed);

        int *matrizNotas = (int*) malloc(R * C * A * sizeof(int));
        estogram *cidades = calloc(R * C, sizeof(estogram)),
        *regioes = calloc(R, sizeof(estogram)),
        *Brasil = calloc(sizeof(estogram));

        // Gera matriz de notas
        for(i=0; i<R; i++){
            for(j=0; j<C; j++){
                for(int k=0; k<A; k++){
                    matrizNotas[i*C*A + j*A + k] = rand()%MAX_VAL;
                    printf("%d ", matrizNotas[i*C*A + j*A + k]);
                }
                printf("\n");
            }
            printf("\n");
        }

        // tempoExec = MPI_Wtime

        nRegProc = R/npes;
        int resto = R%npes;
        MPI_Request request[npes];

        int regAtual = nRegProc;
        if(resto > 0)
            regAtual = nRegProc + 1;

        // Distribui trabalho
        for(int i=1; i<npes; i++){ //# Pode ser quebrado em dois fors antes e dps do tamanho do resto
            sizes[i][1] = C;
            sizes[i][2] = A;
            if(i < resto){
                sizes[i][0] = nRegProc+1; 
            } else {
                sizes[i][0] = nRegProc;
            }
            
            MPI_Isend(sizes[i], 3, MPI_INT, i, 0, MPI_COMM_WORLD, &request[i]);
            MPI_Isend(&matrizNotas[regAtual*C*A], sizes[i][0]*C*A, MPI_INT, i, 1, MPI_COMM_WORLD, &request[i]);

            regAtual += sizes[i][0];
        }

        if(resto > 0)
            nRegProc++;

        for(int i=0; i<nRegProc*C; i++){
            for(int j=0; j<A; j++){
                cidades[i][matrizNotas[i*A + j]]++;
            }
        }

        // Struct de metricas (media, dp e mediana)
        double *dados_cidades = (double *) malloc(R * C * 3 * sizeof(double));
        double *dados_regiao = (double *) malloc(R * 3 * sizeof(double));
        double dados_brasil[3];

        // Vetores com maiores e menores notas
        int *maior_cidades = (int *) malloc(R * C * sizeof(int));
        int *menor_cidades = (int *) malloc(R * C * sizeof(int));
        int *maior_regiao = (int *) malloc(R * sizeof(int));
        int *menor_regiao = (int *) malloc(R * sizeof(int));
        int maior_brasil, menor_brasil;
        
        // Cidade/Regiao com melhor media        
        int melhor_cidade[2], melhor_regiao;
        double maiorMediaCidade = -1, maiorMediaRegiao = -1;

        // Calcula metricas de cada Região/Cidade e país
        for(i=0; i<nRegProc; i++) {
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
            calcula_metricas(regioes[i], *Brasil, dados_regiao + i, C * A);

            // Encontra maior/menor nota
            maior_regiao[i] = emax(regioes[i]);
            menor_regiao[i] = emin(regioes[i]);

            // Atualiza melhor região e melhor média
            if(dados_regiao[i*3 + MEDIA] > maiorMediaRegiao){
                maiorMediaRegiao = dados_regiao[i*3 + MEDIA];
                melhor_regiao = i;
            }
        }

        int auxBrasil[npes*MAX_VAL];

        int nRegs;
        regAtual = nRegProc;
        if(resto > 0)
            regAtual = nRegProc + 1;
        MPI_Request requests[npes][6];
        
        for(int i=1; i<npes; i++){
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

        MPI_Gather(*Brasil, MAX_VAL, MPI_INT, auxBrasil, MAX_VAL, MPI_INT, 0, MPI_COMM_WORLD);

        int melhores_cidades[npes * 2], melhores_regioes[npes];
        double maioresMediasCidade[npes], maioresMediasRegiao[npes];

        MPI_Gather(melhor_cidade, 2, MPI_INT, melhores_cidades, 2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&melhor_regiao, 1, MPI_INT, melhores_regioes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&maiorMediaCidade, 1, MPI_DOUBLE, maioresMediasCidade, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&maiorMediaRegiao, 1, MPI_DOUBLE, maioresMediasRegiao, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        printf("asjksdkfsçf");
        fflush(stdout);
        
        // regAtual = 0;
        // for(int i=0; i<npes; i++){
        //     if(maioresMediasCidade[i] > maiorMediaCidade){
        //         maiorMediaCidade = maioresMediasCidade[i];
        //         melhor_cidade[CIDADE] = melhores_cidades[i + CIDADE];
        //         melhor_cidade[REGIAO] = melhores_cidades[i + REGIAO] + regAtual;
        //     }

        //     if(maioresMediasRegiao[i] > maiorMediaRegiao){
        //         maiorMediaRegiao = maioresMediasRegiao[i];
        //         melhor_regiao = melhores_regioes[i] + regAtual;
        //     }

        //     if(i < resto){
        //         nRegs = nRegProc+1; 
        //     } else {
        //         nRegs = nRegProc;
        //     }
        //     regAtual += nRegs;
        // }

        for(int i=1; i<npes; i++){
            for(int j=0; j<MAX_VAL; j++){
                (*Brasil)[j] += auxBrasil[i*MAX_VAL + j];
            }
        }

        // printf("Hist Brasil\n");
        // for(int i=0; i<MAX_VAL; i++){
        //     printf("%d ", (*Brasil)[i]);
        // }
        // printf("\n");

        calcula_metricas_pais(*Brasil, dados_brasil, R*C*A);

        // Encontra maior/menor nota
        maior_brasil = emax(*Brasil);
        menor_brasil = emin(*Brasil);

        MPI_Status status;

        for(int i=1; i<npes; i++){
            for(int j=0; j<6; j++){
                MPI_Wait(&requests[i][j], &status);
            }
        }
        
        int k;
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

        // Métricas do país
        fprintf(arq, "Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, dados_brasil[MEDIANA], dados_brasil[MEDIA], dados_brasil[DESVIO_PADRAO]);

        fprintf(arq, "\n");

        fprintf(arq, "Melhor regiao: Regiao %d\n", melhor_regiao);
        fprintf(arq, "Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade[REGIAO], melhor_cidade[CIDADE]);

        // printf("Melhores cidades\n");
        // for(int i=0; i<npes; i++){
        //     printf("processo %d\n", i);
        //     printf("reg %d cid %d\n", melhores_cidades[i*2 + REGIAO], melhores_cidades[i*2 + CIDADE]);
        //     printf("melhor reg %d\n", melhores_regioes[i]);
        //     printf("maior media Cidade %.3lf\n", maioresMediasCidade[i]);
        //     printf("maior media Regiao %.3lf\n", maioresMediasRegiao[i]);
        //     printf("\n");
        // }

    } else {
        int sizes[3];
        MPI_Recv(sizes, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);

        int R = sizes[0], C = sizes[1], A = sizes[2];
        int *matrizNotas = (int *)malloc(R*C*A*sizeof(int));
        estogram *cidades = calloc(1, R * C * sizeof(estogram)),
        *regioes = calloc(1, R * sizeof(estogram)),
        *Brasil = calloc(1, sizeof(estogram));

        MPI_Recv(matrizNotas, R*C*A, MPI_INT, 0, 1, MPI_COMM_WORLD, NULL);

        fprintf(arq, "Rank = %d\n", myrank);
        for(i=0; i<R; i++){
            for(j=0; j<C; j++){
                for(int k=0; k<A; k++){
                    fprintf(arq, "%d ", matrizNotas[i*C*A + j*A + k]);
                }
                fprintf(arq, "\n");
            }
            fprintf(arq, "\n");
        }

        for(int i=0; i<R*C; i++){
            for(int j=0; j<A; j++){
                cidades[i][matrizNotas[i*A + j]]++;
            }
        }

        // for(int i=0; i<R*C; i++){
        //     for(int j=0; j<MAX_VAL; j++){
        //         fprintf(arq, "%d ", cidades[i][j]);
        //     }
        //     fprintf(arq, "\n");
        // }

        // Struct de metricas (media, dp e mediana)
        double *dados_cidades = (double *) malloc(R * C * 3 * sizeof(double));
        double *dados_regiao = (double *) malloc(R * 3 * sizeof(double));
        double *dados_brasil;

        // Vetores com maiores e menores notas
        int *maior_cidades = (int *) malloc(R * C * sizeof(int));
        int *menor_cidades = (int *) malloc(R * C * sizeof(int));
        int *maior_regiao = (int *) malloc(R * sizeof(int));
        int *menor_regiao = (int *) malloc(R * sizeof(int));
        int menor_brasil, maior_brasil;
        
        // Cidade/Regiao com melhor media
        int melhor_cidade[2], melhor_regiao;
        double maiorMediaCidade = -1, maiorMediaRegiao = -1;

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
            calcula_metricas(regioes[i], *Brasil, dados_regiao + i, C * A);

            // Encontra maior/menor nota
            maior_regiao[i] = emax(regioes[i]);
            menor_regiao[i] = emin(regioes[i]);

            // Atualiza melhor região e melhor média
            if(dados_regiao[i*3 + MEDIA] > maiorMediaRegiao){
                maiorMediaRegiao = dados_regiao[i*3 + MEDIA];
                melhor_regiao = i;
            }
        }

        MPI_Request requests[6];
        MPI_Isend(dados_cidades, R * C * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(dados_regiao, R * 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &requests[1]);

        MPI_Isend(menor_cidades, R*C, MPI_INT, 0, 2, MPI_COMM_WORLD, &requests[2]);
        MPI_Isend(maior_cidades, R*C, MPI_INT, 0, 3, MPI_COMM_WORLD, &requests[3]);
        MPI_Isend(menor_regiao, R, MPI_INT, 0, 4, MPI_COMM_WORLD, &requests[4]);
        MPI_Isend(maior_regiao, R, MPI_INT, 0, 5, MPI_COMM_WORLD, &requests[5]);

        int auxBrasil[npes*MAX_VAL];

        MPI_Gather(*Brasil, MAX_VAL, MPI_INT, auxBrasil, MAX_VAL, MPI_INT, 0, MPI_COMM_WORLD);

        int melhores_cidades[npes * 2], melhores_regioes[npes];
        double maioresMediasCidade[npes], maioresMediasRegiao[npes];

        MPI_Gather(melhor_cidade, 2, MPI_INT, melhores_cidades, 2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&melhor_regiao, 1, MPI_INT, melhores_regioes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&maiorMediaCidade, 1, MPI_DOUBLE, maioresMediasCidade, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&maiorMediaRegiao, 1, MPI_DOUBLE, maioresMediasRegiao, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        int k;
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

        // Métricas do país
        //fprintf(arq, "Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, dados_brasil, dados_brasil.med, dados_brasil.dp);

        //fprintf(arq, "\n");

        fprintf(arq, "Melhor regiao: Regiao %d\n", melhor_regiao);
        fprintf(arq, "Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade[REGIAO], melhor_cidade[CIDADE]);
    }

	// // Calcula media, dp e mediana do país
	// calcula_metricas_pais(*Brasil, &dados_brasil, R * C * A);

	// // Encontra maior/menor nota
	// maior_brasil = emax(*Brasil);
	// menor_brasil = emin(*Brasil);

    // tempoExec = omp_get_wtime() - tempoExec;

	// Imprime resultados

	// Métricas das Cidades
	// int k;
	// for (i = 0; i < R; i++)
	// {
	// 	for(j = 0; j < C; j++)
	// 	{
	// 		k = i*C+j;
	// 		fprintf(arq, "Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], dados_cidades[k].median, dados_cidades[k].med, dados_cidades[k].dp);
	// 	}
	// 	fprintf(arq, "\n");
	// }

	// // Métricas das regiões
	// for (i = 0; i < R; i++)
	// {
	// 	fprintf(arq, "Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", i, menor_regiao[i], maior_regiao[i], dados_regiao[i].median, dados_regiao[i].med, dados_regiao[i].dp);
	// }

	// fprintf(arq, "\n");

	// // Métricas do país
	// //fprintf(arq, "Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f e DP: %.2f\n", menor_brasil, maior_brasil, dados_brasil.median, dados_brasil.med, dados_brasil.dp);

	// //fprintf(arq, "\n");

	// fprintf(arq, "Melhor regiao: Regiao %d\n", melhor_regiao);
	// fprintf(arq, "Melhor cidade: Regiao %d, Cidade %d\n", melhor_cidade_reg, melhor_cidade);

	//printf("\nTempo de resposta sem considerar E/S, em segundos: %.3lfs\n", tempoExec);

    fclose(arq);
    MPI_Finalize();

	return(0);
}
