#include <stdio.h>
#include <malloc.h>
#include <math.h>

typedef float **tabela;
typedef float *vetorA;
typedef float *vetorB;

float Lagrange(int n, tabela pontos, float ponto) {
	float *L, P = 0, num = 1, den = 1;
	int i, j;
	
	L = (float*) malloc(n * sizeof(float));
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i != j) {
				num = num * (ponto - pontos[0][j]);
				den = den * (pontos[0][i] - pontos[0][j]);
			}
		}
		L[i] = num / den;
		num = 1;
		den = 1;
		P = P + (pontos[1][i] * L[i]);
	}
	return P;
}



float Newton(int n, tabela pontos, float ponto) {
	float **dif, num, den, *anterior, produto, P;
	int i, j, lim, soma = 1;
	
	lim = n - 1;
	
	dif = (float**) malloc((n - 1) * sizeof(float*));
	anterior = (float*) malloc(n * sizeof(float));
	for(i = 0; i < n; i++) {
		dif[i] = (float*) malloc(n * sizeof(float));
		anterior[i] = pontos[1][i];
	}
	
	for(i = 0; i < (n - 1); i++) {
		for(j = 0; j < lim; j++) {
			dif[i][j] = (anterior[j+1] - anterior[j]) / (pontos[0][j + soma] - pontos[0][j]);
			anterior[j] = dif[i][j];
		}
		soma++;
		lim--;
	}
	
	P = pontos[1][0];
	lim = 1;
	produto = 1;
	
	for(i = 0; i < (n - 1); i++) {
		for(j = 0; j < lim; j++) {
			produto = produto * (ponto - pontos[0][j]);
		}
		lim++;
		P = P + produto * dif[i][0];
		produto = 1;
	}
	
	return P;
}



int Fatorial(int n) {
	if(n == 1)
		return 1;
	else {
		return n * Fatorial(n - 1);
	}
}

float NewtonGregory(int n, tabela pontos, float ponto) {
	float **dif;
	float P;
	float h;
	float diferenca;
	float *anterior;
	float produto;
	int i;
	int j;
	int lim;
	
	lim = n - 1;
	
	produto = 1;
	
	//Checa se números são espacados igualmente
	if(n >= 2) {
		h = pontos[0][1] - pontos[0][0];
		
		for(i = 0; i < (n - 1); i++) {
			diferenca = pontos[0][i+1] - pontos[0][i];
			if(diferenca != h) {
				printf("\n Os valores nao sao espacados igualmente! \n");
				return (0);
			}
		}
	}
	else {
		printf("\n Valor unico! \n");
		return (0);
	}
	
	//Constroi tabela das diferencas
	dif = (float**) malloc((n - 1) * sizeof(float*));
	for(i = 0; i < (n - 1); i++) {
		dif[i] = (float*) malloc(n * sizeof(float));
	}
	
	//Define o y como anterior
	anterior = (float*) malloc(n * sizeof(float));
	for(i = 0; i < n; i++) {
		anterior[i] = pontos[1][i];
	}
	
	//Contabiliza tabela de diferencas
	for(i = 0; i < (n - 1); i++) {
		for(j = 0; j < lim; j++) {
			dif[i][j] = anterior[j+1] - anterior[j];
			//printf("(%f - %f = %f)", anterior[j+1], anterior[j], dif[i][j]);
			anterior[j] = dif[i][j];
			//printf("%f ", dif[i][j]);
		}
		printf("\n");
		lim--;
	}
	
	P = pontos[1][0];
	lim = 1;
	
	for(i = 0; i < (n - 1); i++) {
		for(j = 0; j < lim; j++) {
			produto = produto * (ponto - pontos[0][j]);
			printf("%f     ", produto);
		}
		lim++;
		P = P + produto * (dif[i][0] / (Fatorial(i + 1) * pow(h, (i + 1))));
		produto = 1;
	}
	
	return (P);
}

float CoefDeterminacao(int n, tabela pontos, float *vetorY) {
	int i;
	float e, quadrado, y, yquadrado, calculo;
	
	e = 0;
	y = 0;
	yquadrado = 0;
	
	for(i = 0; i < n; i++) {
		quadrado = pontos[1][i] - vetorY[i];
		e = e + (quadrado * quadrado);
		
		y = y + (pontos[1][i]) * (pontos[1][i]);
		
		yquadrado = yquadrado + pontos[1][i];
	}
	
	yquadrado = yquadrado * yquadrado;
	
	calculo = 1 - ((n * e) / (n * y - yquadrado));
	
	return (calculo);
}


int main() {
	tabela pontos;
	int i;
	float *vetorY;
	
	pontos = (float**) malloc(2 * sizeof(float*));
	vetorY = (float*) malloc(5 * sizeof(float));
	for(i = 0; i < 2; i++) {
		pontos[i] = (float*) malloc(4 * sizeof(float));
	}
	
	pontos[0][0] = 0;
	pontos[0][1] = 0.2;
	pontos[0][2] = 0.5;
	pontos[0][3] = 0.8;
	pontos[1][0] = 0;
	pontos[1][1] = 0.1987;
	pontos[1][2] = 0.4794;
	pontos[1][3] = 0.7174;
	
	vetorY[0] = 1.9302;
	vetorY[1] = 1.7093;
	vetorY[2] = 1.4884;
	vetorY[3] = 1.0466;
	vetorY[4] = 0.8257;
	
	//printf("%f", Lagrange(4, pontos, 2));
	//printf("%f", NewtonGregory(4, pontos, 5));
	//printf("%f", CoefDeterminacao(5, pontos, vetorY));
	printf("%f", Newton(4, pontos, 0.35));
	
	return 0;
}
