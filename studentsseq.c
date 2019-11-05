#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include <time.h>

#define swap(i, j) aux = i, i = j, j = aux;

int max(int* vector, int n){
    int m = vector[0];
    for(int i=1;i<n;i++){
        if(vector[i] > m){
            m = vector[i];
        }
    }
    return m;
}

int min(int* vector, int n){
    int m = vector[0];
    for(int i=1;i<n;i++){
        if(vector[i] < m){
            m = vector[i];
        }
    }
    return m;
}

int quickoddmedian(int* vet, int high, int target) {
    int pivo = vet[0], aux, i = 0, rnd = rand() % high;
    swap(vet[0], vet[rnd]);

    for (int j = 1; j < high; j++)
        if (vet[j] < pivo) {
            ++i;
            swap(vet[i], vet[j]);
        }

    swap(vet[0], vet[i]);
    if (target == i)
        return vet[i];

    if (target < i)
        return quickoddmedian(vet, i, target);
    else
        return quickoddmedian(vet + i + 1, high - i - 1, target - i - 1);
}

int quickpairmedian(int* vet, int high, int target) {
    int pivo = vet[0], aux, i = 0, rnd = rand() % high;
    swap(vet[0], vet[rnd]);

    for (int j = 1; j < high; j++)
        if (vet[j] < pivo) {
            ++i;
            swap(vet[i], vet[j]);
        }

    swap(vet[0], vet[i]);
    if (target == i - 1)
        return (float) (vet[i] + min(vet + i + 1, high - i - 1)) / 2;

    if (target == i)
        return (float) (vet[i] + max(vet, i + 1)) / 2;

    if (target < i)
        return quickpairmedian(vet, i, target);
    else
        return quickpairmedian(vet + i + 1, high - i - 1, target - i - 1);
}

void counting_sort(int *vector, int n) {
	int i, min = INT_MAX, max = INT_MIN;
	int *Aux, *R;

	// Encontrar a chave mínima e a máxima
	for (i = 0; i < n; i++) {			// c1 * n
		if (vector[i] > max) max = vector[i];
		if (vector[i] < min) min = vector[i];
	}

	// Declarar o vetor auxiliar (max-min+1)
	Aux = (int *) calloc(max-min+1, sizeof(int));   // c2 * (max-min+1)

	// Contar a partir do vetor original, quantas vezes cada chave ocorre
	for (i = 0; i < n; i++) Aux[vector[i]-min]++;  // c3 * n

	// Função cumulativa
	for (i = 1; i < max-min+1; i++) Aux[i] += Aux[i-1];   // c4 * (max-min)

	// Ordenação ocorre aqui!
	R = (int *) malloc(sizeof(int) * n);

	for (i = n-1; i >= 0; i--) R[--Aux[vector[i]-min]] = vector[i]; // c5 * n
	for (i = 0; i < n; i++) vector[i] = R[i];		// c6 * n

	free(Aux);
	free(R);
}

// Ordena lin linhas passadas (ou seja, ordena cidades, regioes e o brasil dependendo do valor de col)
void ordena_linhas(int *matriz, int lin, int col)
{
    int j;

    for (j = 0; j < lin; j++)
    {
        counting_sort(&matriz[j*col], col);
    }
}

void encontra_maiores(int *matriz, int *maiores, int lin, int col){
    for(int i=0; i<lin; i++){
        maiores[i] = matriz[col*i + col-1];
    }
}

void encontra_menores(int *matriz, int *menores, int lin, int col){
    for(int i=0; i<lin; i++){
        menores[i] = matriz[col*i];
    }
}

// Mediana de cada linha (cidade, região ou o brasil, dependendo do valor de col)
void calcula_mediana(int *matriz, float *vet, int lin, int col)
{
    int j;

    for (j = 0; j < lin; j++)
    {
        vet[j] = matriz[col*j + col/2];
        if(!(col % 2))
        {
            vet[j] += matriz[col*j + col/2-1];
            vet[j] = (float )vet[j]*0.5;
        }
    }
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
    //     printf("Unable to open file!");
    //     exit(1);
    // }

    while(scanf("%d", &lines[i]) != EOF)
        i++;

    int R = lines[0];
    int C = lines[1];
    int A = lines[2];
    int seed = lines[3];

    /// Generates the matrix
    srand(seed);
    int *matriz = (int *)malloc(R*C*A*sizeof(int));
    for (i = 0; i < R*C*A; i++)
    {
        matriz[i] = rand() % 100;
    }

    // Print array in matrix's representation    
    // for (i = 0; i < R*C; i++)
    // {
    //     for (j = 0; j < A; j++)
    //         printf("%d ", matriz[i*A+j]);
    //     printf("\n");
    // }    
   
    /// Calculate the metrics
    /// Mean and deviation
    float *media_cidades = (float *) malloc(R*C*sizeof(float));
    float *media_regiao = (float *) malloc(R*sizeof(float));
    float media_brasil;

    float *dp_cidades = (float *) malloc(R*C*sizeof(float));
    float *dp_regiao = (float *) malloc(R*sizeof(float));
    float dp_brasil;

    /// Median, max and min
    float *mediana_cidades = (float *) malloc(R*C*sizeof(float));
    float *mediana_regiao = (float *) malloc(R*sizeof(float));
    float mediana_brasil;

    int *maior_cidades = (int *) malloc(R*C*sizeof(int));
    int *maior_regiao = (int *) malloc(R*sizeof(int));
    int maior_brasil = -1;

    int *menor_cidades = (int *) malloc(R*C*sizeof(int));
    int *menor_regiao = (int *) malloc(R*sizeof(int));
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
    for(int i=0; i<R; i++){
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
            dp_cidades[i*C + j] = sqrt((s2  -  (float) s1 * s1 / A)    /    (A - 1));

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
        dp_regiao[i] = sqrt((s4  -  (float) s3 * s3 / (C*A))    /    (C*A - 1));

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
    dp_brasil = sqrt((s6  -  (float) s5 * s5 / (R*C*A))    /    (R*C*A - 1));
    nclocksMedia = clock() - nclocksMedia;

    // Mediana da cidade
    nclocksSortingC = clock();
    // ordena_linhas(matriz, R*C, A);
    // encontra_maiores(matriz, maior_cidades, R*C, A);
    // encontra_menores(matriz, menor_cidades, R*C, A);
    // calcula_mediana(matriz, mediana_cidades, R*C, A);
    if(!(A%2)){
        for(int i=0; i<R; i++){
            for(int j=0; j<C; j++){
                mediana_cidades[i*C + j] = quickpairmedian(&matriz[i*C*A + j*A], A, A/2);
            }
        }
    } else {
        for(int i=0; i<R; i++){
            for(int j=0; j<C; j++){
                mediana_cidades[i*C + j] = quickoddmedian(&matriz[i*C*A + j*A], A, A/2);
            }
        }
    }
    
    nclocksSortingC = clock() - nclocksSortingC;

    // printf("\n");
    // for (i = 0; i < R*C; i++)
    // {
    //     for (j = 0; j < A; j++)
    //         printf("%d ", matriz[i*A+j]);
    //     printf("\n");
    // }

    // Mediana da regiao
    nclocksSortingR = clock();
    ordena_linhas(matriz, R, A*C);
    // encontra_maiores(matriz, maior_regiao, R, A*C);
    // encontra_menores(matriz, menor_regiao, R, A*C);
    calcula_mediana(matriz, mediana_regiao, R, A*C);
    nclocksSortingR = clock() - nclocksSortingR;

    // printf("\n");
    // for (i = 0; i < R*C; i++)
    // {
    //     for (j = 0; j < A; j++)
    //         printf("%d ", matriz[i*A+j]);
    //     printf("\n");
    // }

    // Mediana do Brasil
    nclocksSortingP = clock();
    ordena_linhas(matriz, 1, A*C*R);
    // encontra_maiores(matriz, &maior_brasil, 1, A*C*R);
    // encontra_menores(matriz, &menor_brasil, 1, A*C*R);
    calcula_mediana(matriz, &mediana_brasil, 1, A*C*R);
    nclocksSortingP = clock() - nclocksSortingP;

    tempoExec = (clock() - tempoExec)/CLOCKS_PER_SEC;

    // printf("\n");
    // for (i = 0; i < R*C; i++)
    // {
    //     for (j = 0; j < A; j++)
    //         printf("%d ", matriz[i*A+j]);
    //     printf("\n");
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

    printf("\n");

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
