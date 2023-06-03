# Proyecto de Memoria de Título, Septiembre 2022
ESTIMACIÓN DE DISTANCIA JENSEN-SHANNON USANDO SKETCHES APLICADA EN ANÁLISIS GENÓMICO

El análisis y la comparación de secuencias genómicas es un área de investigación que va de la mano con la implementación de nuevos métodos capaces de lograr la computación de conjuntos de datos cada vez más masivos. Por esta razón, existe le necesidad de implementar y evaluar algoritmos y estructuras de datos que soporten la comparación de la información genómica de una manera más eficiente. En este contexto, se evalúa la implementación del uso de sketches, que son estructuras de datos probabilísticas con la ventaja de ser más eficientes que las estructuras convencionales. Los sketches aceptan consultas rápidas con baja complejidad lineal y espacio sublineal, no obstante, están sujetos a un error de estimación.
Esta memoria de título presenta la implementación y evaluación de estimación de la métrica de similitud entre distribuciones llamada Distancia Jensen-Shannon, la cual está basada en la estimación de la Entropía de Shannon usando sketches que se propone en el paper "A High-Throughput Hardware Accelerator for
Network Entropy Estimation Using Sketches" de junio de 2021. La métrica de distancia es evaluada en el contexto de la comparación de similitud entre distintos conjuntos de secuencias genómicas, para lo cual se realiza una búsqueda de los parámetros óptimos mediante la experimentación y el análisis de resultados. Los resultados obtenidos son positivos, por un lado, en relación con la precisión de estimación de la métrica, y también en cuanto al uso de recursos de memoria utilizados.
## Autor
Juan Albornoz
##Requisitos y compilación
###Pre-requisitos
Compilador gcc/g++ y Python >= 3.7
El proyecto fue corrido en Ubuntu 20.04.4 LTS (GNU/Linux 5.15.90.1-microsoft-standard-WSL2 x86_64), pero puede ser compilado y ejecutado análogamente desde windows mediante, por ejemplo, el paquete MinGW.
###Compilación
g++ -o estimacionJSD est_JSD.cpp count_min_sketch.cpp PQ.cpp hll_sketch.h
###Ejecución y parámetros
./estimacionJSD genoma1.fna genoma2.fna k_lenght p_bits pq_height pq_width cu_width cu_depth
*k_length: Largo k de la de la subcadena contenida dentro del conjunto de todos los k-mers en una secuenciación genómica. En el contexto de esta memoria se utilizan secuencias de k-mer de largo 10 y 15, esto es, k=10 y k=15.
*p_bits: Parámetro de precisión del HyperLogLog sketch, y se relaciona a la dimensión del vector del sketch involucrado en el algoritmo de estimación de cardinalidad.
*pq_height: 
*pq_width:
*cu_width:
*cu_depth:

GCF_0041.fna GCF_0019.fna 10 12 4 2 13 3
