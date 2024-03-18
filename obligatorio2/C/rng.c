#include <stdio.h>
#include <time.h>
#include "gsl_rng.h"  // Se usará para generar números aleatorios

gsl_rng *tau;  // Se define globalmente


int rng(int seed, int fSize)
{
	// Inicializamos el generador de números aleatorios gsl_rng
	extern gsl_rng *tau;
	tau = gsl_rng_alloc(gsl_rng_taus); 
	gsl_rng_set(tau, seed);  

	// Generamos un número aleatorio entre 0 y 1:
	int x = gsl_rng_uniform(tau);
	printf("%i", x);

	// Abrimos el fichero en el que escribiremos los números aleatorios
	FILE *rngFile;
	rngFile = fopen("fichero.txt", "w");

	// Rellenamos el fichero con números aleatorios
	int i = 0;
	for(int j=0; j<fSize; j++)
	{
		gsl_rng_set(tau, seed+j);  // Genera el número aleatorio (y cambia la semilla con cada iteración)
		fprintf(rngFile, "%f/n", gsl_rng_uniform(tau));  // Introduce el número aleatorio en el fichero
	}
}


int main(void)
{
	// Input de la semilla
	int seed;
	printf("Introduce la semilla del generador de números aleatorios: ");
	scanf("%i", &seed);

	// Input del tamaño del fichero
	int fSize;
	printf("Introduce el tamaño del fichero: ");
	scanf("%i", &fSize);

	// Aquí llamamos a la función rng() para generar un número aleatorio
	rng(seed, fSize);

	return 0;
}