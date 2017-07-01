#include <stdio.h>
#include <math.h>
#include <malloc.h>

float Determinante(int n, float **m){
	int sinal = 1;
	int i, j, k;
	float det = 0;
	float **A;

	if(n == 1){
		return m[0][0];
	}
	else {

		A =(float**) malloc((n-1)*sizeof(float*));
		for (i=0;i<n;i++)
		{
			A[i] = (float*) malloc((n-1)*sizeof(float));
		}

		for(k=0;k<n;k++){

			for (j=0;j<k;j++){
				for(i=1;i<n;i++){
					A[i-1][j] = m[i][j];
				}
			}
			for (j=k+1;j<n;j++){
				for(i=1;i<n;i++){
					A[i-1][j-1] = m[i][j];
				}
			}

			det += sinal*(m[0][k]*Determinante(n-1, A));
			sinal *= -1;
		}
	}

	return det;
}

void SistemaTriangularInferior(int ordem, float **matriz, float *termos, float **solucao) {
	int i, j;
	float s = 0;
	
	for(i = 0; i < ordem; i++) {
		s = termos[i];
		for(j = 0; j < i; j++) {
			s = s - (matriz[i][j] * (*solucao)[j]);
		}
		(*solucao)[i] = s/matriz[i][i];
	}
}

void SistemaTriangularSuperior(int ordem, float **matriz, float *termos, float **solucao) {
	int i, j;
	float s = 0;

	for(i = ordem - 1; i >= 0; i--) {
		s = termos[i];
		for(j = ordem - 1; j > i; j--) {
			s = s - (matriz[i][j] * (*solucao)[j]);
		}
		(*solucao)[i] = s/matriz[i][i];
	} 
}

void DecomposicaoLU(int ordem, float **matriz, float *termos, float **solucao) {
	int i, j, linha, coluna, uc, lc, uf, lf, k, n, o = 0;
	float **m, **l, **u, *y, soma;
	
	for(i = ordem; i >= 1; i--) {
		if(Determinante(i, matriz) == 0) {
			printf("\nSistema nao converge!\n");
			o = 1;
			break;
		}
	}
	
	if(o == 0) {
		u = (float**) malloc(ordem * sizeof(float*));
		for(i = 0; i < ordem; i++) {
			u[i] = (float*) malloc(ordem * sizeof(float));
			for(j = 0; j < ordem; j++) {
				u[i][j] = 0;
			}
		}
		
		l = (float**) malloc(ordem * sizeof(float*));
		for(i = 0; i < ordem; i++) {
			l[i] = (float*) malloc(ordem * sizeof(float));
			for(j = 0; j < ordem; j++) {
				l[i][j] = 0;
			}
		}
		
		soma = 0;
		
		uc = 0;
		lc = 0;
		
		linha = 0;
		coluna = 0;
		
		for(n = 0; n < ordem; n++) {
			for(j = uc; j < ordem; j++) {
				for(k = 0; k <= linha - 1; k++) {
					soma = soma + l[linha][k] * u[k][j];
				}
				u[linha][j] = matriz[linha][j] - soma;
				soma = 0;
			}
			linha = linha + 1;
			uc = uc + 1;
			
			for(i = lc; i < ordem; i++) {
				for(k = 0; k <= coluna - 1; k++) {
					soma = soma + l[i][k] * u[k][coluna];
				}
				l[i][coluna] = (matriz[i][coluna] - soma) / u[coluna][coluna];
				soma = 0;
			}
			coluna = coluna + 1;
			
			lc = lc + 1;
		}
		
		y = (float*) malloc(ordem * sizeof(float));
		
		SistemaTriangularInferior(ordem, l, termos, &y);
		
		SistemaTriangularSuperior(ordem, u, y, solucao);
	}
}

void GaussCompacto(int ordem, float **matriz, float **termos, float **solucao) {
	int i, j, o = 0, n, uc, lc, k, linha, coluna;
	float **a, **u, **l, soma;
	
	for(i = ordem; i >= 1; i--) {
		if(Determinante(i, matriz) == 0) {
			printf("\nSistema nao converge!\n");
			o = 1;
			break;
		}
	}
	
	if(o == 0) {
		a = (float**) malloc(ordem * sizeof(float*));
		for(i = 0; i < ordem; i++) {
			a[i] = (float*) malloc((ordem + 1) * sizeof(float));
			for(j = 0; j < ordem; j++) {
				a[i][j] = matriz[i][j];
			}
			a[i][ordem] = (*termos)[i];
		}
	
		u = (float**) malloc(ordem * sizeof(float*));
		for(i = 0; i < ordem; i++) {
			u[i] = (float*) malloc((ordem + 1) * sizeof(float));
			for(j = 0; j < ordem + 1; j++) {
				u[i][j] = 0;
			}
		}
		
		l = (float**) malloc(ordem * sizeof(float*));
		for(i = 0; i < ordem; i++) {
			l[i] = (float*) malloc(ordem * sizeof(float));
			for(j = 0; j < ordem; j++) {
				l[i][j] = 0;
			}
		}
		
		soma = 0;
	
		uc = 0;
		lc = 1;
	
		for(n = 0; n < ordem; n++) {
			for(j = uc; j < ordem + 1; j++) {
				for(k = 0; k <= n - 1; k++) {
					soma = soma + l[n][k] * u[k][j];	
				}
				u[n][j] = a[n][j] - soma;
				soma = 0;
			}
			uc = uc + 1;
		
			for(i = lc; i < ordem; i++) {
				for(k = 0; k <= n - 1; k++) {
					soma = soma + l[i][k] * u[k][n];
				}
				l[i][n] = (a[i][n] - soma) / u[n][n];
				soma = 0;
			}
			lc = lc + 1;
		}
	
		for(i = 0; i < ordem; i++) {
			(*termos)[i] = u[i][ordem];
		}
	
		SistemaTriangularSuperior(ordem, u, *termos, solucao);
	}
}






void GaussSeidel(int ordem, float **matriz, float *termos, float e, float *inicial, int maximo, float **solucao, int *iteracoes) {
	float soma, a, d, max, *b, *x, numerador, denominador, erro, *anterior;
	int i, j, o = 0, k;
	
	soma = 0;
	
	for(i = 0; i < ordem; i++) {
		if(matriz[i][i] == 0) {
			printf("\nSistema nao converge!\n");
			o = 1;
			break;
		}
	}
	
	if(o == 0) {
		d = Determinante(ordem, matriz);
		if(d == 0) {
			printf("\nSistema nao converge!\n");
			o = 1;
		}
		if(o == 0) {
			for(i = 0; i < ordem; i++) {
				for(j = 0; j < ordem; j++) {
					if(j != i) {
						a = (matriz[i][j] / matriz[i][j]);
						if(a < 0) {
							a = a * (-1);
						}
						soma = soma + a;
					}
				}
				if(i == 0) {
					max = soma;
				}
				else {
					if(soma > max) {
						max = soma;
					}
				}
				soma = 0;
			}
			if(max >= 1) {
				o = 2;
			}
		}
		if(o == 2) {
			b = (float*) malloc(ordem * sizeof(float));
			for(i = 0; i < ordem; i++) {
				soma = 0;
				for(j = 0; j < i - 1; j++) {
					a = matriz[i][j] / matriz[i][i];
					if(a < 0) {
						a = a * (-1);
					}
					soma = soma + a * b[j];
				}
				for(j = i + 1; j < ordem; j++) {
					a = matriz[i][j] / matriz[i][i];
					if(a < 0) {
						a = a * (-1);
					}
					soma = soma + a;
				}
				b[i] = soma;
				if(i == 0) {
					max = b[i];
				}
				else {
					if(b[i] > max) {
						max = b[i];
					}
				}				
			}
			if(max < 1) {
				o = 0;
			}
			free(b);
		}
		if(o == 0) {
			
			x = inicial;
			k = 0;
			erro = e;
			max = -1;
			anterior = (float*) malloc(ordem * sizeof(float));
			do {
				for(i = 0; i < ordem; i++) {
					x[i] = termos[i];
					for(j = 0; j < i; j++) {
						x[i] = x[i] - (matriz[i][j]) * x[j]; 
					}
					for(j = i + 1; j < ordem; j++) {
						x[i] = x[i] - (matriz[i][j]) * x[j];
					}
					x[i] = x[i] / matriz[i][i];
				}
				if(k > 0) {
					for(i = 0; i < ordem; i++) {
						numerador = x[i] - anterior[i];
						if(numerador < 0) {
							numerador = numerador * (-1);
						}
						if(numerador > max) {
							max = numerador;
						}
					}
					numerador = max;
					max = -1;
					for(i = 0; i < ordem; i++) {
						denominador = x[i];
						if(denominador < 0) {
							denominador = denominador * (-1);
						}
						if(denominador > max) {
							max = denominador;
						}
					}
					denominador = max;
					max = -1;
					erro = numerador/denominador;
				}
				
				k++;
				for(i = 0; i < ordem; i++) {
					anterior[i] = x[i];
				}
			} while((erro >= e)&&(k + 1 <= maximo));
			
			(*solucao) = x;
			(*iteracoes) = k;
		}
	}
}





int main() {
	int i;
	int ordem;
	float **matriz;
	float *termos;
	float *solucao;
	float determinante;
	float *inicial;
	int iteracoes;

	ordem = 3;
	matriz = (float**)malloc(ordem*sizeof(int));
	
	for(i = 0; i < ordem; i++) {
		matriz[i] = (float*) malloc(ordem * sizeof(float));
	}
	matriz[0][0] = 5;
	matriz[0][1] = 1;
	matriz[0][2] = 1;
	matriz[1][0] = 3;
	matriz[1][1] = 4;
	matriz[1][2] = 1;
	matriz[2][0] = 3;
	matriz[2][1] = 3;
	matriz[2][2] = 6;
	
	solucao = (float*) malloc(ordem*sizeof(float));
	
	termos = (float*) malloc(ordem*sizeof(float));
	termos[0] = 5;
	termos[1] = 6;
	termos[2] = 0;
	
	inicial = (float*) malloc(ordem * sizeof(float));
	inicial[0] = 0;
	inicial[1] = 0;
	inicial[2] = 0;
	
	//DecomposicaoLU(ordem, matriz, termos, &solucao);
	
	//GaussCompacto(ordem, matriz, &termos, &solucao);
	
	GaussSeidel(ordem, matriz, termos, 0.01, inicial, 5, &solucao, &iteracoes);
	
	for(i = 0; i < ordem; i++) {
		printf("%f ", solucao[i]);
	}
	
	printf("Iteracoes = ", iteracoes);
	
	return (0);
}
