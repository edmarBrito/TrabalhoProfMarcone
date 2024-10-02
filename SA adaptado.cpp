//SA Adaptado

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "Util.h"
#include "Construcao.h"
#include "Arquivos.h"
#include "Descida.h"
#include "SimulatedAnnealing.h"



float SimulatedAnnealing(int n,
                         int *s,
                         float **d,
                         float alpha,
                         int SAmax,
                         float temp_inicial,
                         float temp_final)
{
  float f_star, fo, temp, delta, fo_viz;
  int *s_star, iter, i, j, aux;
  clock_t inicio_CPU, fim_CPU;

  s_star = cria_vetor(n);
  atualiza_vetor(s_star, s, n);
  f_star = fo = calcula_fo(n,s_star,d);
  temp = temp_inicial;

  char SAsaida[] = "SA.txt";

  limpa_arquivo(SAsaida);
  inicio_CPU = fim_CPU = clock();
  imprime_fo(SAsaida, (fim_CPU - inicio_CPU)/CLOCKS_PER_SEC,f_star,0);
  printf("f_star = %8.2f \t Temperatura = %8.2f \n", f_star, temp);
  while (temp > temp_final){
    iter = 0;
    while (iter < SAmax){
      iter++;
      // Escolhe um vizinho qualquer
      i = rand() % (n);
      do{
        j = rand() % (n);
      }while (j == i);

      // Calcula o custo das arestas envolvidas ANTES do movimento
      float delta1 = calcula_delta(n,s,d,i,j);

      // Faz o movimento
      aux = s[i];
      s[i] = s[j];
      s[j] = aux;

      // Calcula o custo das arestas envolvidas DEPOIS do movimento
      float delta2 = calcula_delta(n,s,d,i,j);

      fo_viz = fo - delta1 + delta2;
      // Calcula a varia��o de energia
      delta = fo_viz - fo;
      if (delta < 0){
        fo = fo_viz;
        if (fo < f_star){
          f_star = fo;
          atualiza_vetor(s_star, s, n);
          f_star = descida_primeiro_melhora(n, s_star, d); //Refinamento da Solução com First Improvement
          fim_CPU = clock();
          imprime_fo(SAsaida, (fim_CPU - inicio_CPU)/CLOCKS_PER_SEC,f_star,0);
          printf("f_star = %8.2f \t Temperatura = %8.2f \n", f_star, temp);
        }
      }
      else{
        double x;
        x = randomico(0,1);
        if (x < exp(-delta/temp)){
          fo = fo_viz;
        }
        else{
          /* Desfaz o movimento caso o vizinho nao seja de melhora
             ou n�o passe no teste de aceita��o */
          aux = s[i];
          s[i] = s[j];
          s[j] = aux;
        }
      }
    } // final de SAmax iteracoes
    temp = temp * alpha;
  } // temperatura de congelamento do sistema
  atualiza_vetor(s, s_star, n);
  fim_CPU = clock();
  imprime_fo(SAsaida, (fim_CPU - inicio_CPU)/CLOCKS_PER_SEC,f_star,0);
  return f_star;
}


float calcula_temperatura_inicial(int n,
                                  int *s,
                                  float **d,
                                  float beta,  // taxa de reaquecimento
                                  float gamma, // taxa de resfriamento
                                  int SAmax)
{
  int i, j, aux, iterT, aceitos;
  float temperatura, delta;
  bool continua;

  temperatura = 10;  // chute inicial para a temperatura inicial
  continua = true;
  while (continua){
    aceitos = 0;
    iterT = 0;
    while (iterT < SAmax){
      iterT++;
      // Gere um vizinho qualquer
      i = rand() % (n);
      do{
        j = rand() % (n);
      } while (j == i);
      float delta1 = calcula_delta(n, s, d, i, j);
      // Faz o movimento
      aux = s[j];
      s[j] = s[i];
      s[i] = aux;
      float delta2 = calcula_delta(n, s, d, i, j);
      delta = - delta1 + delta2;
      if (delta < 0){
         aceitos++;
      }
      else{
        float x = randomico(0,1);
        if (x < exp(-delta/temperatura)){
           aceitos++;
        }
      }
      // Desfaz o movimento
      aux = s[j];
      s[j] = s[i];
      s[i] = aux;
    }
    if (aceitos < gamma * SAmax)
      temperatura = beta * temperatura;
    else
      continua = false;
  }
  return temperatura;
}