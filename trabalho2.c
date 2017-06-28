#include <stdio.h>
#include <math.h>
#include <malloc.h>

float Determinante(int ordem, float **matriz) {
	float soma;
	float cofator;
	float **m;
	int i;
	int j;
	int k;
	int f;
	
	soma = 0;
	
	if(ordem == 1) {
		return (matriz[0][0]);
	}
	else {
		m = (float**) malloc((ordem - 1) * sizeof(float*));
		for(i = 0; i < ordem - 1; i++) {
			m[i] = (float*) malloc((ordem - 1) * sizeof(float));	
		}
		
		f = 0;
		
		for(k = 0; k < ordem; k++) {
			for(i = 0; i < ordem - 1; i++) {
				for(j = 0; j < ordem; j++) {
					if(f != j) {
						m[i][j] = matriz[i][j];
					}
				}
			}

			f = f + 1;
			cofator = pow(-1, k) * Determinante(ordem - 1, m);
			soma = soma + matriz[ordem - 1][k] * cofator;
			return (soma);
		}
	}
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
	
	/*for(i = ordem; i >= 1; i--) {
		if(i != ordem) {
			m = (float**) malloc(i * sizeof(float*));
			for(j = 0; j < i; j++) {
				m[j] = (float*) malloc(i * sizeof(float));
			}
		}
		else {
			m = matriz;
		}
		if(Determinante(i, m) == 0) {
			printf("\nO sistema nao converge!\n");
			o = 1;
			break;
		}
	}*/
	
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

void GaussCompacto(int ordem, float **matriz, float *termos, float **solucao) {
	
}

int main() {
	int i;
	int ordem;
	float **matriz;
	float *termos;
	float *solucao;
	float determinante;

	ordem = 3;
	matriz = (float**)malloc(ordem*sizeof(int));
	
	for(i = 0; i < ordem; i++) {
		matriz[i] = (float*) malloc(ordem * sizeof(float));
	}
	matriz[0][0] = 5;
	matriz[0][1] = 2;
	matriz[0][2] = 1;
	matriz[1][0] = 3;
	matriz[1][1] = 1;
	matriz[1][2] = 4;
	matriz[2][0] = 1;
	matriz[2][1] = 1;
	matriz[2][2] = 3;
	
	solucao = (float*) malloc(ordem*sizeof(float));
	
	termos = (float*) malloc(ordem*sizeof(float));
	termos[0] = 0;
	termos[1] = -7;
	termos[2] = -5;
	
	DecomposicaoLU(ordem, matriz, termos, &solucao);
	
	for(i = 0; i < ordem; i++) {
		printf("%f ", solucao[i]);
	}
	
	return (0);
}
