# Fortran95 Special Functions


-Bessel functions
-Error function
-Exponential integrals
-Incomplete Gamma


(BUILDING)


FUNCIONES IMPLEMENTADAS EN LA LIBRERÍA / NOMBRE DE MÓDULOS
--------------------------------------
   
 * Funciones Trigonométricas    	- Complex_Hyperbolic_Trigonometric_Functions
 					- Complex_Inverse_Hyperbolic_Trigonometric_Functions
 					- Complex_Inverse_Trigonometric_Functions
 									
 * Integrales Exponenciales		- Exponential_Integrals  (Necesita el módulo Incomplete_Gamma para compilar)
 
 * Función Error			- Error_Function
 
 * Función Incompleta de Gamma		- Incomplete_Gamma
 
 * Funciones de Bessel			- Bessel_Functions
 
 
UTILIZACIÓN
--------------------------------------

 * Importar módulo al código principal (main):

		use [module_name.f90]
		

 * Colocar el módulo en la ubicación correcta
 
 * Indicar la localización del módulo al compilar
 
 
  OTROS
--------------------------------------
 
 * Los nombres de los módulos son los mismos que los de las carpetas que lo contienen 
 
 * Para más información, recurrir al informe de Prácticas
 