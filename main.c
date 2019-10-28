#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int maior_gen(int *matriz, int R, int C, int A, int *maior_cidades, int *maior_regiao)
{
    int i, j, maior = -1, maior_brasil = -1;
    int cidades = R*C;
    int start, end;

    for (i = 0; i < cidades; i++)
    {
        start = A * i;
        end = start + A;
        for (j = start; j < end; j++)
        {
            if(matriz[j] > maior)
                maior = matriz[j];
        }
        maior_cidades[i] = maior;
        maior = -1;
    }

    for (i = 0; i < R; i++)
    {
        start = i*C;
        end = start+C;
        for (j = start; j < end; j++)
        {
            if(maior_cidades[j] > maior)
                maior = maior_cidades[j];
        }
        maior_regiao[i] = maior;
        maior = -1;
        if(maior_regiao[i] > maior_brasil)
            maior_brasil = maior_regiao[i];
    }

    return maior_brasil;
}

int menor_gen(int *matriz, int R, int C, int A, int *menor_cidades, int *menor_regiao)
{
    int i, j, menor = 101, menor_brasil = 101;
    int cidades = R*C;
    int start, end;

    for (i = 0; i < cidades; i++)
    {
        start = A * i;
        end = start + A;
        for (j = start; j < end; j++)
        {
            if(matriz[j] < menor)
                menor = matriz[j];
        }
        menor_cidades[i] = menor;
        menor = 101;
    }

    for (i = 0; i < R; i++)
    {
        start = i*C;
        end = start+C;
        for (j = start; j < end; j++)
        {
            if(menor_cidades[j] < menor)
                menor = menor_cidades[j];
        }
        menor_regiao[i] = menor;
        menor = 101;
        if(menor_regiao[i] < menor_brasil)
            menor_brasil = menor_regiao[i];
    }

    return menor_brasil;
}

int partition (int *arr, int low, int high, int C)
{
    int i, j;
    int pivot,swap;

    // pivot (Element to be placed at right position)
    pivot = arr[high];

    i = (low - 1);  // Index of smaller element

    for (j = low; j <= high-1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j] <= pivot)
        {
            i++;    // increment index of smaller element

            // swap arr[i] and arr[j]
            swap = arr[i];
            arr[i] = arr[j];
            arr[j] = swap;
        }
    }

    //swap arr[i + 1] and arr[high]
    swap = arr[(i + 1)];
    arr[(i + 1)] = arr[high];
    arr[high] = swap;

    return (i + 1);

} // fim partition

void quicksort(int *arr, int low, int high, int C)
{
    int pi;
    if (low < high)
    {
        /* pi is partitioning index, arr[pi] is now
           at right place */
        pi = partition(arr, low, high, C);

        quicksort(arr, low, pi - 1, C);  // Before pi
        quicksort(arr, pi + 1, high, C); // After pi
    }
} // fim quicksort

// Ordena tam linhas passadas (ou seja, ordena cidades, regioes e o brasil)
void ordena_linhas(int *matriz, int lin, int col, int tam)
{
    int j;
    int inferior, superior;
    for (j = 0; j < lin/tam; j++)
    {
        inferior = tam*j*col;
        superior = tam*(j*col+col)-1;
        //manda o endereco do primeiro elemento da coluna, limites inf e sup e a largura da matriz
        quicksort(&matriz[0], inferior, superior, tam*col);
    }
}

void calcula_mediana(int *matriz, float *vet, int lin, int col, int tam)
{
    int j;
    for (j = 0; j < lin/tam; j++)
    {
        vet[j] = matriz[tam*col*j + tam*col/2];
        if(!(lin%2))
        {
            vet[j]+=matriz[tam*col*j + tam*col/2-1];
            vet[j]= (float )vet[j]*0.5;
        }
    }
}

// Tem erro propagado pq eu fiz media das medias
float calcula_media(int *matriz, float *vet_cidade, float *vet_reg, int R, int C, int A)
{
    int i,j;
    int start, end;
    float soma = 0;

    for(i=0; i<R*C; i++)
    {
        for(j=0; j<A; j++)
        {
            soma+=matriz[i*A+j];
        }
        vet_cidade[i]=(float) soma/A;
        soma=0;
    }

    float soma_reg = 0;
    for (i = 0; i < R; i++)
    {
        start = i*C;
        end = start+C;
        for (j = start; j < end; j++)
        {
            soma+=vet_cidade[j];
        }
        vet_reg[i] = (float) soma/C;
        soma_reg+=vet_reg[i];
        soma = 0;
    }

    return soma_reg/R;
}

void calcula_variancia(double *matriz, double *media,double *variancia, int lin, int col)
{
    int i,j;
    double soma;
    for(i=0;i<col;i++){
        soma=0;
        for(j=0;j<lin;j++){
            soma+=pow((matriz[j*col+i]-media[i]),2);
        }
        variancia[i]=soma/(lin-1);
    }
}

void calcula_desvio_padrao(double *variancia,double *dp, int col)
{
    int i;
    for(i=0;i<col;i++){
        dp[i]=sqrt(variancia[i]);
    }
}

int main(void)
{
    /// Read param from file
    int lines[4], i = 0, j;

    FILE *fp = fopen("entrada.txt", "r");
    if(fp == NULL)
    {
        printf("Unable to open file!");
        exit(1);
    }

    while(fscanf(fp, "%d", &lines[i]) != EOF)
        i++;

    int R = lines[0];
    int C = lines[1];
    int A = lines[2];
    int seed = lines[3];

    /// Generates the matrix
    int *matriz = (int *)malloc(R*C*A*sizeof(int));
    int *matrizOrdenada = (int *)malloc(R*C*A*sizeof(int));

    srand(seed);

    for (i = 0; i < R*C*A; i++)
    {
        matriz[i] = rand() % 100;
        matrizOrdenada[i] = matriz[i];
    }

    // Print array in matrix's representation
    for (i = 0; i < R*C; i++)
    {
        for (j = 0; j < A; j++)
            printf("%d ", matriz[i*A+j]);
        printf("\n");
    }

    /// Calculate the metrics

    /// Max and Min
    int *maior_cidades = (int *) malloc(R*C*sizeof(int));
    int *maior_regiao = (int *) malloc(R*sizeof(int));
    int *menor_cidades = (int *) malloc(R*C*sizeof(int));
    int *menor_regiao = (int *) malloc(R*sizeof(int));
    int maior_brasil = maior_gen(matriz, R, C, A, maior_cidades, maior_regiao);
    int menor_brasil = menor_gen(matriz, R, C, A, menor_cidades, menor_regiao);

    /// Median

    /*
    printf("\n\n\n");
    for (i = 0; i < R*C; i++)
    {
        for (j = 0; j < A; j++)
            printf("%d ", matrizOrdenada[i*A+j]);
        printf("\n");
    }
    */
    float *mediana_cidades = (float *) malloc(R*C*sizeof(float));
    float *mediana_regioes = (float *) malloc(R*sizeof(float));
    float mediana_brasil;

    ordena_linhas(matrizOrdenada, R*C, A, 1);
    calcula_mediana(matrizOrdenada, mediana_cidades, R*C, A,1);

    ordena_linhas(matrizOrdenada, R*C, A, C);
    calcula_mediana(matrizOrdenada, mediana_regioes, R*C, A, C);

    ordena_linhas(matrizOrdenada, R*C, A, R*C);
    calcula_mediana(matrizOrdenada, &mediana_brasil, R*C, A, R*C);



    /// Mean
    float *media_cidades = (float *) malloc(R*C*sizeof(float));
    float *media_regiao = (float *) malloc(R*sizeof(float));
    float media_brasil = calcula_media(matriz,media_cidades,media_regiao,R,C,A);

    /// Deviation


    /// Best city and region
    float maior = -1;
    int cidade, regiao_cidade, regiao;
    for(i = 0; i < R; i++){
        if(media_regiao[i] > maior){
            maior = media_regiao[i];
            regiao = i;
        }
    }
    maior = -1;
    for(i = 0; i < R*C; i++){
        if(media_cidades[i] > maior){
            maior = media_cidades[i];
            cidade = i;
        }
    }
    regiao_cidade = cidade/4;
    cidade = cidade%C;


    // Metrics for cities
    int k;
    for (i = 0; i < R; i++)
    {
        for(j = 0; j < C; j++)
        {
            k = i*C+j;
            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f\n", i, j, menor_cidades[k], maior_cidades[k], mediana_cidades[k], media_cidades[k]);
        }
        printf("\n");
    }

    // Metrics for regions
    for (i = 0; i < R; i++)
    {
        printf("Reg %d: menor: %d, maior: %d, mediana: %.2f, media: %.2f\n", i, menor_regiao[i], maior_regiao[i], mediana_regioes[i], media_regiao[i]);
    }

    printf("\n");

    // Metrics for the country
    printf("Brasil: menor: %d, maior: %d, mediana: %.2f, media: %.2f\n", menor_brasil, maior_brasil, mediana_brasil, media_brasil);

    printf("\n");

    printf("Melhor regiao: Regiao %d\n", regiao);
    printf("Melhor cidade: Regiao %d, Cidade %d", regiao_cidade, cidade);
    fclose(fp);
    return(0);

}
