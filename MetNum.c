
#include <stdio.h>
#include <malloc.h>
#include <math.h>


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
	int i, j, linha, coluna, uc, lc, k, n, o = 0;
	float **l, **u, *y, soma;
	
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
		
		
		for(i = 0; i < ordem; i++) {
			printf("%f   ", (*solucao)[i]);
		}
	}
}

void GaussCompacto(int ordem, float **matriz, float **termos, float **solucao) {
	int i, j, o = 0, n, uc, lc, k;
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
		
		for(i = 0; i < ordem; i++) {
			printf("%f   ", (*solucao)[i]);
		}
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

			printf("\nx%d = {", *iteracoes);
			for(i=0;i<ordem;i++){
				printf( " %f ",x[i]);
			} 
			printf("}");
		}
	}
}

void GaussJordan(int n, float **s, float *b, float *x){
	int i, j, k;
	float m;
	//Calcular determinante Calcular determinante Calcular determinante Calcular determinante Calcular determinante
	float **matriz =(float**) malloc(n*sizeof(float*));
	for (i=0;i<n;i++)
	{
		matriz[i] = (float*) malloc((n+1)*sizeof(float));
	}

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			matriz[i][j] = s[i][j];
		}
		matriz[i][j] = b[i];
	}

	for(i=0;i<n;i++){
		for (j=0;j<i;j++){
			m = matriz[j][i]/matriz[i][i];
			for(k=i;k<=n;k++){
				matriz[j][k] = matriz[j][k] - m*matriz[i][k];
			}
		}
		for (j=i+1;j<n;j++){
			m = matriz[j][i]/matriz[i][i];
			for(k=i;k<=n;k++){
				matriz[j][k] = matriz[j][k] - m*matriz[i][k];
			}
		}
	}

	for(i=0;i<n;i++){
		x[i] = matriz[i][n]/matriz[i][i];
	}
}

void Jacobi(int n, float **s, float *b, float e, float *x0, int k, float *x, int *ite){
	int i, j;
	float max = 0, max2 = 0, aux = 0;
	char criterio = 0;
	float *anterior, *atual;

	float **matriz =(float**) malloc(n*sizeof(float*));
	for (i=0;i<n;i++)
	{
		matriz[i] = (float*) malloc((n)*sizeof(float));
	}

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			matriz[i][j] = s[i][j];
		}
	}

	for(i=0;i<n;i++){
		if(s[i][i] == 0){
			printf("\nImpossivel realizar operacao. (Divisao por zero)\n");
			return;
		} 
	}

	if(!Determinante(n,s)){
		printf("\nImpossivel realizar operacao. (Determinante igual a zero)\n");
		return;
	}
	//crit?io das linhas
	for(i=0;i<n;i++){
		for (j=0;j<i;j++){
			aux += fabs(matriz[i][j]/matriz[i][i]);
		}
		for (j=i+1;j<n;j++){
			aux += fabs(matriz[i][j]/matriz[i][i]);
		}
		if(aux > max) max = aux;
		aux = 0;
	}
	if(max < 1) criterio = 1;
	else{ //crit?io das colunas
		max = 0;
		for(j=0;j<n;j++){
			for (i=0;i<j;i++){
				aux += fabs(matriz[i][j]/matriz[j][j]);
			}
			for (i=j+1;i<n;i++){
				aux += fabs(matriz[i][j]/matriz[j][j]);
			}
			if(aux > max) max = aux;
			aux = 0;
		}
		if(max < 1) criterio = 1;
	}

	if(criterio){
		atual = (float*) malloc((n)*sizeof(float));
		anterior = x0;

		for((*ite)=0;(*ite)<k;(*ite)++){
			for(i=0;i<n;i++){
				atual[i] = b[i];
				for (j=0;j<i;j++){
					atual[i] -= matriz[i][j]*anterior[j]; 
				}
				for (j=i+1;j<n;j++){
					atual[i] -= matriz[i][j]*anterior[j]; 
				}
				atual[i] /= matriz[i][i];
			}

			max = 0;
			for(i=0;i<n;i++){
				aux = fabs(atual[i] - anterior[i]);
				if(aux > max) max = aux;
			}
			for(i=0;i<n;i++){
				aux = fabs(atual[i]);
				if(aux > max2) max2 = aux;
			}

			if(max/max2 < e) {
				for(i=0;i<n;i++){
					x[i] = atual[i];
				}
				return;
			}

			for(i=0;i<n;i++){
				anterior[i] = atual[i];
			}
		}
		for(i=0;i<n;i++){
			x[i] = atual[i];
		}
		printf("\nx%d = {", *ite);
		for(j=0;j<n;j++){
			printf( " %f ",x[j]);
		} 
		printf("}");
	}
	else printf("\nNao satisfaz criterio.\n");
}

void MatrizInversa(int n, float **m, float ** inv){
	int i, j, op;
	float **id = (float**) malloc(n*sizeof(float*));
	float *sol = (float*) malloc(n*sizeof(float));
	for(i=0;i<n;i++){
		id[i] = (float*) malloc(n*sizeof(float)); 
	}
	for(i=0;i<n;i++){
		for(j=0;j<i;j++){
			id[i][j] = 0; 
		}
		id[i][i] = 1;
		for(j=i+1;j<n;j++){
			id[i][j] = 0; 
		} 
	}
	printf("Determinar a inversa utilizando o Metodo da Decomposicao LU ou o Metodo de Gauss Compacto?\n1.Decomposicao LU\n2.Metodo de Gauss Compacto\n");
	fflush(stdin);
	scanf("%d", &op);
	if(op == 1){
		for(i=0;i<n;i++){
			printf("Coluna %d = {",i+1);
			DecomposicaoLU(n,m,id[i],&sol);
			for(j=0;j<n;j++){
				inv[j][i] = sol[j];
			} 
			printf(" }");
		}
	}
	if(op == 2){
		for(i=0;i<n;i++){
			printf("Coluna %d = {",i+1);
			GaussCompacto(n,m,&id[i],&sol);
			for(j=0;j<n;j++){
				inv[j][i] = sol[j];
			} 
			printf(" }\n");
		}
	}	
}

int main(){

	int i, j, ite, k, op = 15;
	int n;
	float *b, *x0, *x, **A, **Inv, e;

	while(op != 0) {
		switch(op) {
			case 1:
				printf("Digite a ordem da matriz:");
				fflush(stdin);
				scanf("%d", &n);
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]);
					} 
				}
				printf("\nDeterminante:%f",Determinante(n,A));
				op = 15;
				break;
			case 2:
				printf("Digite a ordem da matriz de coeficientes:");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						if(i>=j){
							printf("\nDigite o elemento a%d%d da matriz de coeficientes [NOTA: lembre que i>=j]:",i+1,j+1);
							fflush(stdin);
							scanf("%f", &A[i][j]);
						}
						else A[i][j] = 0;
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				} 
				SistemaTriangularInferior(n,A,b,&x);
				printf("\nx = {");
				for(j=0;j<n;j++){
					printf(" %f  ",x[j]);
				} 
				printf("}");
				op = 15;
				break;
			case 3:
					printf("Digite a ordem da matriz::");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						if(i<=j){
							printf("\nDigite o elemento a%d%d da matriz de coeficientes [NOTA: lembre que i<=j]:",i+1,j+1);
							fflush(stdin);
							scanf("%f", &A[i][j]);
						}
						else A[i][j] = 0;
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				} 
				SistemaTriangularSuperior(n,A,b,&x);
				printf("\nx = {");
				for(j=0;j<n;j++){
					printf("  %f  ",x[j]);
				} 
				printf("}");
				op = 15;
				break;
			case 4:
				printf("Digite a ordem da matriz de coeficientes:");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz de coeficientes:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]);
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				}
				printf("Solucao = {");
				DecomposicaoLU(n,A,b,&x);
				printf("}");
				op = 15;
				break;
				break;
			case 5:
				printf("Digite a ordem da matriz de coeficientes:");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz de coeficientes:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]);
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				} 
				printf("Solucao = {");
				GaussCompacto(n,A,&b,&x);
				printf("}");
				op = 15;
				break;
			case 6:
				printf("Digite a ordem da matriz de coeficientes:");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz de coeficientes:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]);
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				} 
				GaussJordan(n, A, b, x);
				printf("\nx = {");
				for(j=0;j<n;j++){
					printf( "%f  ",x[j]);
				} 
				printf("}");
				op = 15;
				break;
			case 7:
				printf("Digite a ordem da matriz de coeficientes:");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				x0 =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz de coeficientes:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]);
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				}
				printf("\nDigite a precisao desejada:");
				fflush(stdin);
				scanf("%f", &e);

				for(j=0;j<n;j++){
					printf("\nDigite o elemento x%d do vetor de aproximação inicial:",j+1);
					fflush(stdin);
					scanf("%f", &x0[j]);
				}

				printf("\nDigite o numero máximo de iteracoes:");
				fflush(stdin);
				scanf("%d", &k);

				Jacobi(n, A, b, e, x0, k, x, &ite);

				op = 15;
				
				break;
			case 8:
				printf("Digite a ordem da matriz de coeficientes:");
				fflush(stdin);
				scanf("%d", &n);
				x =(float*) malloc(n*sizeof(float));
				x0 =(float*) malloc(n*sizeof(float));
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				b = (float*) malloc((n)*sizeof(float));
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz de coeficientes:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]);
					} 
				}
				for(j=0;j<n;j++){
					printf("\nDigite o elemento b%d do vetor dos termos independentes:",j+1);
					fflush(stdin);
					scanf("%f", &b[j]);
				}
				printf("\nDigite a precisao desejada:");
				fflush(stdin);
				scanf("%f", &e);

				for(j=0;j<n;j++){
					printf("\nDigite o elemento x%d do vetor de aproximação inicial:",j+1);
					fflush(stdin);
					scanf("%f", &x0[j]);
				}

				printf("\nDigite o numero máximo de iteracoes:");
				fflush(stdin);
				scanf("%d", &k);

				GaussSeidel(n, A, b, e, x0, k, &x, &ite);

				op = 15;
				break;
			case 9:
				printf("Digite a ordem da matriz:");
				fflush(stdin);
				scanf("%d", &n);
				A =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					A[i] = (float*) malloc((n)*sizeof(float));
				}
				Inv =(float**) malloc(n*sizeof(float*));
				for (i=0;i<n;i++)
				{
					Inv[i] = (float*) malloc((n)*sizeof(float));
				}
				for(i=0;i<n;i++){
					for(j=0;j<n;j++){
						printf("\nDigite o elemento a%d%d da matriz:",i+1,j+1);
						fflush(stdin);
						scanf("%f", &A[i][j]); 				
					}
				}
				MatrizInversa(n,A,Inv);

				op = 15;
				break;
			default:
				printf("\n\n");
				printf("Opcoes\n");
				printf("1. Rotina Determinante\n");
				printf("2. Rotina SistemaTriangularInferior\n");
				printf("3. Rotina SistemaTriangularSuperior\n");
				printf("4. Rotina DecomposicaoLU\n");
				printf("5. Rotina GaussCompacto\n");
				printf("6. Rotina GaussJordan\n");
				printf("7. Rotina Jacobi\n");
				printf("8. Rotina GaussSeidel\n");
				printf("9. Rotina MatrizInversa\n");
				printf("0. SAIR\n");
				printf("Opcao: ");
				fflush(stdin);
				scanf("%d", &op);
				break;
		}
	}

	return 0;
}
