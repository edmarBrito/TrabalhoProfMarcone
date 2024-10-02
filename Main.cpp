/*
Tecnicas Heuristicas para resolucao do Problema do Caixeiro Viajante
Autor: Marcone Jamilson Freitas Souza
Local: Departamento de Computacao - Universidade Federal de Ouro Preto
Homepage: www.decom.ufop.br/prof/marcone
Data: 21-05-2007
*/

//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "Util.h"
#include "Arquivos.h"
#include "Estruturas.h"
#include "Construcao.h"
#include "Listas.h"
#include "Menus.h"
#include "Descida.h"
#include "SimulatedAnnealing.h"
#include "ILS.h"
#include "AG.h"
#include "GRASP.h"
//---------------------------------------------------------------------------


int main(int argc, char* argv[])
{
  int n,                    // numero de cidades
      *s;                   // vetor solucao corrente
  float **d,                // matriz de distancias entre as cidades
        fo,                 // funcao objetivo corrente
        melhor_fo_lit;      // melhor fo da literatura
  clock_t inicio_CPU,       // clock no inicio da aplicacao do metodo
          fim_CPU;          // clock no final da aplicacao do metodo


  char C50INFO[] = "C50INFO.TXT";
  char C50[] = "C50.TXT";

  obter_parametros_pcv(C50INFO, &n, &melhor_fo_lit);
  s = cria_vetor(n); // vetor solucao
  d = cria_matriz_float(n, n); // matriz de distancias
  le_arq_matriz(C50, n, d);

  int escolha = 0;
  do {
    escolha = menu_principal();
    switch (escolha) {
    case 1: /* Geracao de uma solucao inicial */
           switch(menu_solucao_inicial()) {
           case 1: /* Geracao gulosa de uma solucao inicial via metodo do vizinho mais proximo */
                 constroi_solucao_gulosa_vizinho_mais_proximo(n,s,d);
                 fo = calcula_fo(n, s, d);
                 printf("\nSolucao construida de forma gulosa (Vizinho Mais Proximo):\n");
                 imprime_rota(s, n);
                 printf("Funcao objetivo = %f\n",fo);
	         break;
           case 2: /* Geracao parcialmente gulosa de uma solucao inicial via metodo do vizinho mais proximo */
                 constroi_solucao_parcialmente_gulosa_vizinho_mais_proximo(n,s,d,0.05);
                 fo = calcula_fo(n, s, d);
                 printf("\nSolucao construida de forma gulosa (Vizinho Mais Proximo):\n");
                 imprime_rota(s, n);
                 printf("Funcao objetivo = %f\n",fo);
	         break;
           case 3: /* Geracao gulosa de uma solucao inicial via metodo da insercao mais barata */
                 constroi_solucao_gulosa_insercao_mais_barata(n,s,d);
                 fo = calcula_fo(n, s, d);
                 printf("\nSolucao construida de forma gulosa (Insercao Mais Barata):\n");
                 imprime_rota(s, n);
                 printf("Funcao objetivo = %f\n",fo);
	         break;
           case 4: /* Geracao parcialmente gulosa de uma solucao inicial via insercao mais barata */
                 printf("Ainda nao implementado...\n");
	         break;
           case 5: /* Geracao aleatoria de uma solucao inicial */
                 //srand(1000); // pega o numero 1000 como semente de numeros aleatorios
                 srand((unsigned) time(NULL)); // pega a hora do relogio como semente
                 constroi_solucao(n,s,d);
                 embaralha_vetor(s,n);
                 fo = calcula_fo(n, s, d);
                 printf("\nSolucao construida de forma aleatoria:\n");
                 imprime_rota(s, n);
                 printf("Funcao objetivo = %f\n",fo);
	         break;
           }
           break;
    case 2: /* Descida */
           inicio_CPU = clock();
           fo = descida(n,s,d);
           fim_CPU = clock();
           printf("Tempo execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           printf("Distancia percorrida = %f \n",fo);
           imprime_rota(s,n);
           break;

    case 3: /* Descida Randomica */
           inicio_CPU = clock();
           fo = descida_randomica(n,s,d,0.7*n*(n-1)/2);
           fim_CPU = clock();
           printf("Distancia percorrida = %f \n",fo);
           printf("Tempo execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           imprime_rota(s,n);
           break;

    case 4: /* Descida com primeiro de melhora */
           inicio_CPU = clock();
           fo = descida_primeiro_melhora(n, s, d);
           fim_CPU = clock();
           printf("Distancia percorrida = %f \n",fo);
           printf("Tempo execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           imprime_rota(s,n);
           break;

    case 5: /* Multi-Start */
           srand((unsigned) time(NULL));
           inicio_CPU = clock();
           // fo = MultiStart(n,s,d,100);
           fim_CPU = clock();
           // printf("fo = %10.2f \n",fo);
           printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           // imprime_rota(s,n);
           printf("Nao implementado\n");
           break;

    case 6: /* Simulated Annealing */
           //srand(1000);
           srand((unsigned) time(NULL));
           inicio_CPU = clock();
           float temp_inicial;
           temp_inicial = calcula_temperatura_inicial(n,s,d,1.1,0.95,500);
           fo = SimulatedAnnealing(n,s,d,0.998,10*n,temp_inicial,0.001);
           fim_CPU = clock();
           printf("fo = %10.2f \n",fo);
           printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           imprime_rota(s,n);
           break;

    case 7: /* Busca Tabu */
           inicio_CPU = clock();
           // fo = BT(n,s,d,2,100);
           fim_CPU = clock();
           printf("Tempo execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           printf("Distancia percorrida = %f \n",fo);
           imprime_rota(s,n);
           break;

    case 8: /* Iterated Local Search */
           srand(1000);
           //srand((unsigned) time(NULL));
           inicio_CPU = clock();
           fo = ILS(n,s,d,15,30);
           fim_CPU = clock();
           printf("fo = %10.2f \n",fo);
           printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           imprime_rota(s,n);
           break;

    case 9: /* GRASP */
           switch(menu_GRASP()) {
           case 1: /* Geracao Parcialmente gulosa de uma solucao inicial via metodo do vizinho mais proximo */
                 inicio_CPU = clock();
                 fo = GRASP(n,s,d,0.05,100,1);
                 fim_CPU = clock();
                 printf("\nSolucao gerada pelo Metodo GRASP:\n");
                 imprime_rota(s, n);
                 printf("\nMelhor fo encontrada = %f \n",fo);
                 printf("Tempo de CPU = %f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
                 break;
           case 2: /* Geracao parcialmente gulosa de uma solucao inicial via metodo da insercao mais barata */
	         break;
           }
           break;

    case 10: /* VND */
           printf("Nao implementado\n");
           break;

    case 11: /* VNS */
           srand(1000);
           //srand((unsigned) time(NULL));
           inicio_CPU = clock();
           // fo = VNS(n,s,d,20,10);
           fim_CPU = clock();
           // printf("fo = %10.2f \n",fo);
           printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           printf("Nao implementado\n");
           // imprime_rota(s,n);
           break;

    case 12: /* Algoritmos Geneticos */
           switch(menu_AG()) {
           case 1: /* Algoritmos Geneticos usando operador OX */
                 //srand(1000);
                 srand((unsigned) time(NULL));
                 inicio_CPU = clock();
                 fo = AG(n,s,d,100,0.03,0.85,0.01,1);
                 fim_CPU = clock();
                 printf("Solucao por AG usando operador OX:\n");
                 printf("fo = %10.2f \n",fo);
                 printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
                 imprime_rota(s,n);
                 break;
           case 2: /* Algoritmos Geneticos usando operador ERX */
                 //srand(1000);
                 srand((unsigned) time(NULL));
                 inicio_CPU = clock();
                 fo = AG(n,s,d,100,0.03,0.85,0.01,2);
                 fim_CPU = clock();
                 printf("Solucao por AG usando operador ERX:\n");
                 printf("fo = %10.2f \n",fo);
                 printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
                 imprime_rota(s,n);
                 break;
           }
           break;

    case 13: /* Algoritmos Memeticos */
           srand((unsigned) time(NULL));
           inicio_CPU = clock();
           // fo = Memeticos(n,s,d,100,0.03,50,0.85,0.01,1);
           fim_CPU = clock();
           printf("Solucao por Memeticos usando operador OX:\n");
           // printf("fo = %10.2f \n",fo);
           printf("Tempo de execucao = %10.2f segundos\n",(float)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC);
           // imprime_rota(s,n);
           printf("Nao implementado\n");
           break;

    case 14: /* Colonia de Formigas */
           printf("\n Colonia de Formigas ainda nao implementado ... \n");
           break;

    default:
          break;
    }
  } while (escolha != 0);




  libera_vetor(s);
  libera_matriz_float(d, n);
  return 0;
}
//---------------------------------------------------------------------------
