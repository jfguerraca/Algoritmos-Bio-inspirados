# Algoritmos-Bio-inspirados (En construcción)
Códigos en MATLAB

## Algoritmo Colonia artificial de abejas (Artifitial Bee Colony)

Inspirado en el comportamiento de las abejas en la búsqueda de alimento, el ABC se se basa en los roles de trabajo que toman
estas para completar su búsqueda.

Las abejas obreras explotan las soluciones conocidas, mientras que las exploradoras realizan una búsqueda de nuevas soluciones.
Su implementación depende de ajustar cuatro parámetros que afectan principalmente
la fase de exploración.

> Límite de exploración (L): Si una fuente de alimento no mejora después de L búsquedas, es abandonada, para problemas de gran escala o alta complejidad L = 100.

> Criterio de abandono (Fu,Fl): Establece una condición por medio de un rango para establecer una búsqueda más exhaustiva, a cambio de aumentar el tiempo de ejecución.

> Coeficiente de exploración (prob): Porcentaje que influye si una abeja obrera se convierte en exploradora.

![bee](https://github.com/user-attachments/assets/807749d3-9ce2-4933-b179-95eea86448de)

### Ejemplo de ABC
 > [ABC](https://github.com/jfguerraca/Algoritmos-Bio-inspirados/blob/main/ABC.m)

## Optimización de colonia de hormigas (ACO)

Inspirado en el comportamiento de las
hormigas en la búsqueda de alimento, el
ACO se caracteriza por las tareas y los mecanismos
que estas llevan a cabo con este
fin. Siguiendo los rastros de feromonas
que indican los caminos más rápidos a las
fuentes de alimento.

El mecanismo principal de este algoritmo
consiste en las feromonas.

> Tasa de evaporación (tse): Controla
la cantidad de feromona que se
evapora en cada iteración. Los valores
bajos aumentan la explotación
y los altos la exploración.

> Rastro de feromona (𝛼,𝛽): Establecen
el grado de aleatoriedad en la
que la nueva poblacion de agentes
cambia, afectando la calidad de las
soluciones por lo general 𝛼, 𝛽 ≤ 2.

> Factor de calidad (q): Porcentaje
que influye en la probabilidad de
que una hormiga elĳa un camino
nuevo.

![ant](https://github.com/user-attachments/assets/f9a1a65c-167c-4375-82d2-6db7010e61e2)

### Ejemplo de ACO
 > [ACO](https://github.com/jfguerraca/Algoritmos-Bio-inspirados/blob/main/ACO.m)

## Optimización Hormiga León (ALO)

Inspirado en el comportamiento de caza de las larvas de la hormiga león, el ALO utiliza técnicas de búsqueda aleatoria (como el vuelo de Lévy) para evitar la convergencia prematura y aumentar la exploración de soluciones

![antilion_or](https://github.com/user-attachments/assets/ae690773-b0f1-4db3-8d85-67899d9cd1c5)

### Ejemplo de ALO

 > [ALO](https://github.com/jfguerraca/Algoritmos-Bio-inspirados/blob/main/ALO.m)

## Algoritmo de murciélagos (BA)

Simulando el comportamiento natural de los murciélagos, el BA se basa principalmente en la ecolocación, controlando la frecuencia y amplitud del ultrasonido emitido por los murciélagos, afectando la intensidad de la búsqueda.

![bat](https://github.com/user-attachments/assets/02e0919d-b16b-4c4f-9071-ff9262ed532c)

 > [BA](https://github.com/jfguerraca/Algoritmos-Bio-inspirados/blob/main/BA.m)

## Optimización de la viuda negra (BWO)

El algoritmo BWO es una técnica inspirada en la forma de apareamiento de las arañas viudas negras.

![black](https://github.com/user-attachments/assets/c90f302d-e2f4-4e3c-9a9e-39ddd175373e)

 > [BWO](https://github.com/jfguerraca/Algoritmos-Bio-inspirados/blob/main/BWO.m)




