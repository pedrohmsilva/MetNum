#include <stdio.h>
#include <malloc.h>

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

int main(){

	int i, op;
	int n = 3;
	float b[3] = {4,0,-1};
	float x[3];
	float **s =(float**) malloc(n*sizeof(float*));
	for (i=0;i<n;i++)
	{
		s[i] = (float*) malloc((n)*sizeof(float));
	}

	s[0][0] = 1;
	s[0][1] = 1;
	s[0][2] = 2;
	s[1][0] = 2;
	s[1][1] = -1;
	s[1][2] = -1;
	s[2][0] = 1;
	s[2][1] = -1;
	s[2][2] = -1;
	
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
	

	return 0;
}