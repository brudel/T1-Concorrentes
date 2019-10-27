#include <stdlib.h>
#include <stdio.h>

int maior_gen(int **matriz, int R, int C, int A, int *maior_cidades, int *maior_regiao){
    int i, j, maior = -1, maior_brasil = -1;
    int cidades = R*C;
    int start, end;

    for (i = 0; i < cidades; i++)
    {
        for (j = 0; j < A; j++)
        {
            if(matriz[i][j] > maior)
                maior = matriz[i][j];
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

int menor_gen(int **matriz, int R, int C, int A, int *menor_cidades, int *menor_regiao){
    int i, j, menor = 101, menor_brasil = 101;
    int cidades = R*C;
    int start, end;

    for (i = 0; i < cidades; i++)
    {
        for (j = 0; j < A; j++)
        {
            if(matriz[i][j] < menor)
                menor = matriz[i][j];
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
/*
int partition (double *arr, int low, int high, int C)
{
    int i, j;
    double pivot,swap;

    pivot = arr[high*C];

    i = (low - 1);

    for (j = low; j <= high-1; j++)
    {
        if (arr[j*C] <= pivot)
        {
            i++;
            swap = arr[i*C];
            arr[i*C] = arr[j*C];
            arr[j*C] = swap;
        }
    }
    swap = arr[(i + 1)*C];
    arr[(i + 1)*C] = arr[high*C];
    arr[high*C] = swap;

    return (i + 1);
}

void quicksort(double *arr, int low, int high, int C)
{
    int pi;

    if (low < high)
    {
        pi = partition(arr, low, high, C);

        quicksort(arr, low, pi - 1, C);
        quicksort(arr, pi + 1, high, C);
    }
}

void ordena_colunas(double *matriz, int lin, int col)
{
    int j;

    for (j = 0; j < col; j++)
    {
        //manda o endereco do primeiro elemento da coluna, limites inf e sup e a largura da matriz
        quicksort(&matriz[j], 0, lin-1, col);
    }
}
void calcula_media(double *matriz, double *vet, int lin, int col)
{
    int i,j;
    double soma;

    for(i=0; i<col; i++)
    {
        soma=0;
        for(j=0; j<lin; j++)
        {
            soma+=matriz[j*col+i];
        }
        vet[i]=soma/lin;
    }
}


void calcula_mediana(double *matriz, double *vet, int lin, int col)
{
    int j;
    for (j = 0; j < col; j++)
    {
        vet[j] = matriz[((lin/2)*col)+j];
        if(!(lin%2))
        {
            vet[j]+=matriz[((((lin-1)/2)*col)+j)];
            vet[j]*=0.5;
        }
    }
}

void calcula_variancia(double *matriz, double *media,double *variancia, int lin, int col)
{
    int i,j;
    double soma;
    for(i=0; i<col; i++)
    {
        soma=0;
        for(j=0; j<lin; j++)
        {
            soma+=pow((matriz[j*col+i]-media[i]),2);
        }
        variancia[i]=soma/(lin-1);
    }
}

void calcula_desvio_padrao(double *variancia,double *dp, int col)
{
    int i;
    for(i=0; i<col; i++)
    {
        dp[i]=sqrt(variancia[i]);
    }
}
*/
int main(void)
{

    //double *matrizOrdenada, *mediana,*media,*variancia,*dp; //Define a matriz (forma linear), vetores de medidas estatísticas
    //matrizOrdenada=(double *)malloc(lin*col * sizeof(double)); //Aloca a matriz ORDENADA
    //media=(double *)malloc(col * sizeof(double)); //Aloca o vetor de media
    //mediana=(double *)malloc(col * sizeof(double)); //Aloca o vetor de mediana
    //variancia=(double *)malloc(col * sizeof(double)); //Aloca o vetor de variância
    //dp=(double *)malloc(col * sizeof(double)); //Aloca o vetor de desvio padrão


    FILE *fp = fopen("entrada.txt", "r");
    if(fp == NULL)
    {
        printf("Unable to open file!");
        exit(1);
    }

    int lines[4], i = 0, j;

    while(fscanf(fp, "%d", &lines[i]) != EOF)
    {
        i++;
    }

    int R = lines[0];
    int C = lines[1];
    int A = lines[2];
    int seed = lines[3];
    int cidade;
    int **matriz = (int **)malloc(R*C*sizeof(int*));
    for(i=0; i<R*C; i++)
        matriz[i] = (int *)malloc(A*sizeof(int));

    int **matrizOrdenada = (int **)malloc(R*C*sizeof(int*));
    for(i=0; i<R*C; i++)
        matrizOrdenada[i] = (int *)malloc(A*sizeof(int));

    srand(seed);

    for (i = 0; i < R*C; i++)
    {
        for (j = 0; j < A; j++)
        {
            matriz[i][j] = rand() % 100;
            matrizOrdenada[i][j] = matriz[i][j];
            printf("%d ", matriz[i][j]);
        }
        printf("\n");
    }

    int *maior_cidades = (int *) malloc(R*C*sizeof(int));
    int *maior_regiao = (int *) malloc(R*sizeof(int));
    int *menor_cidades = (int *) malloc(R*C*sizeof(int));
    int *menor_regiao = (int *) malloc(R*sizeof(int));
    int maior_brasil = maior_gen(matriz, R, C, A, maior_cidades, maior_regiao);
    int menor_brasil = menor_gen(matriz, R, C, A, menor_cidades, menor_regiao);

    for (i = 0; i < R; i++)
    {
        for(j = 0; j < C; j++)
        {
            cidade = i*C+j;
            printf("Reg %d - Cid %d: menor: %d, maior: %d\n", i, j, menor_cidades[cidade], maior_cidades[cidade]);
        }
        printf("\n");
    }

    for (i = 0; i < R; i++)
    {
        printf("Reg %d: menor: %d, maior: %d\n", i, menor_regiao[i], maior_regiao[i]);

    }

    printf("\n");

    printf("Brasil: menor: %d, maior: %d\n", menor_brasil, maior_brasil);

    fclose(fp);
    return(0);

}
