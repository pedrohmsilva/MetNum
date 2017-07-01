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
			printf("\nImpossível realizar operação. (Divisão por zero)\n");
			return;
		} 
	}

	if(!Determinante(n,s)){
		printf("\nImpossível realizar operação. (Determinante igual a zero)\n");
		return;
	}
	//critério das linhas
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
	else{ //critério das colunas
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
	}
	else printf("Não satisfaz critério.\n");
}

int main(){

	int i, op;
	int n = 3;
	float b[3] = {14,11,8};
	float chute[3] = {0, 0, 0};
	float x[3];
	float **s =(float**) malloc(n*sizeof(float*));
	for (i=0;i<n;i++)
	{
		s[i] = (float*) malloc((n)*sizeof(float));
	}

	s[0][0] = 10;
	s[0][1] = 2;
	s[0][2] = 1;
	s[1][0] = 1;
	s[1][1] = 5;
	s[1][2] = 1;
	s[2][0] = 2;
	s[2][1] = 3;
	s[2][2] = 10;

	Jacobi(3, s, b, 0.01, chute, 20, x, &i);
	printf("%f\n", x[0]);
	printf("%f\n", x[1]);
	printf("%f\n", x[2]);
	printf("%d\n", i);

	GaussJordan(3, s, b, x);
	printf("\n%f\n", x[0]);
	printf("%f\n", x[1]);
	printf("%f\n", x[2]);



	/*
	while(op != 0) {
		switch(op) {
			case 1:
				break;
			case 2:
				break;
			case 3:
				break;
			case 4:
				break;
			case 5:
				break;
			case 6:
				GaussJordan(n, s, b, x);
				for(i=0;i<n;i++){
					printf("%f\n", x[i]);
				}
				break;
			case 7:
				Jacobi(3, s, b, 0.01, chute, 20, x, &i);
				printf("%f\n", x[0]);
				printf("%f\n", x[1]);
				printf("%f\n", x[2]);
				break;
			case 8:
				break;
			case 9:
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
	*/

	return 0;
}