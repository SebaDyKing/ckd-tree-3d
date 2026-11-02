# cKd-tree tridimensional escalable

Código del Proyecto de Título *Diseño e implementación de un cKd-tree tridimensional escalable
mediante arquitectura de 64 bits y filtrado por cuadrantes*, Ingeniería Civil en Informática,
Universidad del Bío-Bío, 2026.

- Autores: Christian Retamal Reyes y Magdiel Fernández Delgado.
- Profesor guía: Dr. Rodrigo Torres Avilés.

El cKd-tree es un kd-tree compacto: guarda cada nodo como un código espiral respecto de su padre
y comprime los códigos con Códigos Direccionables Directamente (DACs). Su versión tridimensional
dejaba de funcionar en espacios de más de unos 2048 puntos por lado, porque los códigos no cabían
en los 32 bits del DAC. Este proyecto la lleva a 64 bits y agrega una espiral unidireccional que
numera solo el semiespacio del hijo, con lo que funciona en espacios de hasta 1.664.509 celdas por eje.

## Contenido

| Archivo | Qué contiene |
|---|---|
| `src/Dac64.hpp` | DAC de 64 bits: anchos óptimos de bloque en O(n), construcción y acceso O(1). Para valores de 32 bits produce la misma estructura que `dacs.c`, bit a bit. |
| `src/Dac64Extras.hpp` | Descompresión completa, guardar y cargar con validación del archivo. |
| `src/Spiral3dUD.hpp` | Espiral unidireccional: codifica a cada hijo solo en el semiespacio que le corresponde. |
| `src/Spiral3d.hpp` | Espiral bidireccional (cubo completo), usada como comparación. |
| `src/cKDtree3d.hpp` | cKd-tree tridimensional sobre el Dac64, con consultas de punto y de rango. |
| `test/` | Pruebas automáticas (68) y programas de medición. |
| `ref/cKDtreeUD/` | DACs de referencia en C (`dacs.c`, de Brisaboa, Ladra y Navarro), usados para comprobar la equivalencia y reconstruir la versión de Silva Pereira. |

## Uso

Requiere `g++` con C++17 y `bash`. Todo se compila y ejecuta con `Start3d`:

| Comando | Qué hace |
|---|---|
| `./Start3d` | Compila todo y ejecuta las 68 pruebas automáticas |
| `./Start3d pruebas` | Detalle de las pruebas de la espiral |
| `./Start3d dac` | Detalle de las pruebas del Dac64 contra `dacs.c` |
| `./Start3d tiempo` | Costo de codificar y decodificar según la distancia al padre |
| `./Start3d arbol` | Pruebas del árbol y comparación entre espirales e iKd-tree en espacios grandes |
| `./Start3d original` | Comparación con la versión tridimensional de Silva Pereira |
| `K3TREE=ruta ./Start3d k3` | Comparación con el k³-tree (ver abajo) |

Las mediciones repiten cada caso 5 veces con datos nuevos e informan el promedio y la desviación
estándar. Si alguna estructura responde distinto que las demás, la medición se detiene.

### Comparación con el k³-tree

La comparación usa la implementación del k³-tree de Fernando Silva-Coira, que no se incluye en este
repositorio. Se descarga y compila por separado (requiere `cmake`):

```bash
git clone https://gitlab.lbd.udc.es/lbd-open/k3-tree.git
cd k3-tree && sh compile.sh && cd ..
K3TREE="$PWD/k3-tree" ./Start3d k3
```

## Referencias

- G. Gutiérrez, R. Torres-Avilés y M. Caniupán, "cKd-tree: A Compact Kd-tree", *IEEE Access*, vol. 12, pp. 28665–28676, 2024.
- N. Silva Pereira, *Investigación sobre la viabilidad de proponer e implementar una versión multidimensional del cKd-tree*, Memoria de Título, Universidad del Bío-Bío, 2024.
- N. R. Brisaboa, S. Ladra y G. Navarro, "DACs: Bringing Direct Access to Variable-Length Codes", *Information Processing & Management*, vol. 49, no. 1, pp. 392–404, 2013.
- F. Silva-Coira, *k³-tree*, código fuente, Universidade da Coruña: https://gitlab.lbd.udc.es/lbd-open/k3-tree
