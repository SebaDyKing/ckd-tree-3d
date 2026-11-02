/**
 * Benchmark contra el k3-tree de Silva-Coira (gitlab.lbd.udc.es/lbd-open/k3-tree, sobre sdsl-lite).
 * El k3-tree guarda un conjunto de celdas y elimina los puntos repetidos, así que los dos árboles
 * se construyen con los mismos puntos ya sin repetir y la memoria se cuenta por punto distinto.
 * Se mide la memoria de las tres variantes del k3-tree y, para la básica, construcción y
 * consultas de punto y de rango, con la misma función medir<> que los otros benchmarks.
 *
 *   ./ckdtree3d_k3_bench [repeticiones]     por defecto 5
 */

#include "Bench3d.hpp"

#include <k3_tree.hpp>
#include <k3_tree_level.hpp>
#include <k3_tree_points.hpp>

typedef sdsl::k3_tree_base<>::point_type PuntoK3;

/** Lado del cubo del caso que se mide; el k3-tree lo necesita al construirse. */
static int64_t ladoDelCaso = 0;

/** Los puntos en el formato del k3-tree: tuplas (x, y, z) sin signo. */
static std::vector<PuntoK3> aTuplas(const std::vector<Point3D> &puntos) {
    std::vector<PuntoK3> tuplas;
    tuplas.reserve(puntos.size());
    for (const Point3D &p : puntos) tuplas.emplace_back(p.x, p.y, p.z);
    return tuplas;
}

/** k3-tree básico con la interfaz de los árboles de la tesis. */
struct K3TreeBasico {
    std::vector<PuntoK3> tuplas;                                     //solo para construir
    mutable sdsl::k3_tree<> arbol;                                   //get y get_region no son const
    mutable std::vector<sdsl::k3_tree_base<>::tuple_result> region;

    explicit K3TreeBasico(const std::vector<Point3D> &puntos)
        : tuplas(aTuplas(puntos)), arbol(tuplas, (uint64_t)ladoDelCaso) {
        std::vector<PuntoK3>().swap(tuplas);
    }
    bool contains(const Point3D &q) const { return arbol.get(q.x, q.y, q.z); }
    uint64_t countInRange(const Box3D &caja) const {
        region.clear();
        return arbol.get_region(caja.low.x, caja.low.y, caja.low.z, caja.high.x, caja.high.y, caja.high.z, region);
    }
    uint64_t sizeInBits() const { return 8 * sdsl::size_in_bytes(arbol); }
};

/** Bits por punto de otra variante del k3-tree; antes revisa que encuentre los puntos guardados. */
template <typename K3>
static double bitsPorPuntoK3(const std::vector<Point3D> &puntos, const Consultas &consultas, const char *nombre) {
    std::vector<PuntoK3> tuplas = aTuplas(puntos);
    K3 arbol(tuplas, (uint64_t)ladoDelCaso);
    for (const Point3D &q : consultas.presentes)
        if (!arbol.get(q.x, q.y, q.z)) abortarMedicion(nombre, "no encontró todos los puntos guardados");
    return 8.0 * (double)sdsl::size_in_bytes(arbol) / (double)puntos.size();
}

static std::vector<Point3D> sinRepetidos(std::vector<Point3D> puntos) {
    std::sort(puntos.begin(), puntos.end(), menorXYZ);
    puntos.erase(std::unique(puntos.begin(), puntos.end()), puntos.end());
    return puntos;
}

/** Resultados de un caso: la tesis, el k3-tree básico medido completo y la memoria de las variantes. */
struct Caso {
    Resultado tesis, k3;
    Estadistica bitsK3Puntos, bitsK3Niveles, distintos;
};

static void medirCaso(Caso &caso, const std::vector<Point3D> &generados, int64_t lado, std::mt19937_64 &azar) {
    ladoDelCaso = lado;
    const std::vector<Point3D> puntos = sinRepetidos(generados);
    caso.distintos.agregar((double)puntos.size());
    Consultas consultas = generarConsultas(puntos, lado, azar);
    medir<cKDtree3dUD>(puntos, consultas, caso.tesis, "UD + Dac64");
    medir<K3TreeBasico>(puntos, consultas, caso.k3, "k3-tree básico");
    caso.bitsK3Puntos.agregar(bitsPorPuntoK3<sdsl::k3_tree_points<>>(puntos, consultas, "k3-tree con límite de puntos"));
    caso.bitsK3Niveles.agregar(bitsPorPuntoK3<sdsl::k3_tree_level<>>(puntos, consultas, "k3-tree con límite de niveles"));
}

/** Tablas de resumen: memoria, consulta de punto y de rango. */
struct Resumen {
    Tabla memoria, punto, rango;
    Resumen() {
        memoria.encabezados = { "Datos", "Espacio", "Puntos", "Tesis", "k3 básico", "k3 puntos", "k3 niveles",
                                "Tesis/mejor k3" };
        punto.encabezados   = { "Datos", "Espacio", "Puntos", "Tesis (µs)", "k3 (µs)", "Tesis/k3" };
        rango.encabezados   = { "Datos", "Espacio", "Puntos", "Tesis (µs)", "k3 (µs)", "Tesis/k3" };
    }
    void agregar(const char *datos, const std::string &espacio, const std::string &puntos, const Caso &c) {
        const double tesis = c.tesis.bitsPorPunto.promedio();
        const double mejorK3 = std::min({ c.k3.bitsPorPunto.promedio(), c.bitsK3Puntos.promedio(),
                                          c.bitsK3Niveles.promedio() });
        memoria.filas.push_back({ datos, espacio, puntos, numero(tesis, 2), numero(c.k3.bitsPorPunto.promedio(), 2),
                                  numero(c.bitsK3Puntos.promedio(), 2), numero(c.bitsK3Niveles.promedio(), 2),
                                  numero(100.0 * tesis / mejorK3, 1) + "%" });
        punto.filas.push_back({ datos, espacio, puntos, numero(c.tesis.encontrarNs.promedio() * 1e-3, 2),
                                numero(c.k3.encontrarNs.promedio() * 1e-3, 2),
                                numero(100.0 * c.tesis.encontrarNs.promedio() / c.k3.encontrarNs.promedio(), 1) + "%" });
        const int cajaGrande = TAMANOS_CAJA - 1;
        rango.filas.push_back({ datos, espacio, puntos, numero(c.tesis.rangoUs[cajaGrande].promedio(), 1),
                                numero(c.k3.rangoUs[cajaGrande].promedio(), 1),
                                numero(100.0 * c.tesis.rangoUs[cajaGrande].promedio()
                                       / c.k3.rangoUs[cajaGrande].promedio(), 1) + "%" });
    }
};

static void imprimirCaso(const Caso &caso) {
    Tabla tabla = tablaDeResultados();
    tabla.filas.push_back(filaDeResultado("Tesis (UD + Dac64)", caso.tesis));
    tabla.filas.push_back(filaDeResultado("k3-tree básico", caso.k3));
    tabla.imprimir();
    Tabla tablaRango = tablaDeRango();
    tablaRango.filas.push_back(filaDeRango("Tesis (UD + Dac64)", caso.tesis));
    tablaRango.filas.push_back(filaDeRango("k3-tree básico", caso.k3));
    tablaRango.imprimir();
    printf("  Otras variantes del k3-tree: límite de puntos %s bits/punto, límite de niveles %s bits/punto\n",
           conDesviacion(caso.bitsK3Puntos, 2).c_str(), conDesviacion(caso.bitsK3Niveles, 2).c_str());
}

int main(int argc, char **argv) {
    const int repeticiones = argc > 1 ? atoi(argv[1]) : 5;

    printf("\nkd-tree 3D de la tesis (UD + Dac64) contra el k3-tree de Silva-Coira (sdsl-lite)\n");
    printf("Los dos con los mismos puntos sin repetir; bits/punto se cuenta por punto distinto.\n");
    imprimirLeyenda(repeticiones);
    Resumen resumen;
    const Distribucion distribuciones[] = { UNIFORME, GAUSSIANA };

    //PASO 1: los cubos de la comparación con Silva Pereira, mismas semillas
    const int64_t ladosChicos[] = { 64, 128, 256 };
    const double ocupaciones[] = { 0.01, 0.10 };
    for (Distribucion distribucion : distribuciones) {
        for (int64_t lado : ladosChicos) {
            for (double ocupacion : ocupaciones) {
                const uint64_t cantidad = (uint64_t)(ocupacion * (double)lado * (double)lado * (double)lado);
                Caso caso;
                for (int repeticion = 0; repeticion < repeticiones; ++repeticion) {
                    std::mt19937_64 azar((uint64_t)lado * 1000 + (uint64_t)(ocupacion * 100) * 10
                                         + 100000 * distribucion + repeticion);
                    medirCaso(caso, generarEnCubo(distribucion, lado, cantidad, azar), lado, azar);
                }
                const std::string ocupado = numero(100.0 * ocupacion, 0) + "%";
                printf("\n%s, espacio %s, %s de las celdas con un punto (%s puntos, %s distintos)\n",
                       nombreDistribucion(distribucion), espacioLegible(lado).c_str(), ocupado.c_str(),
                       cantidadLegible(cantidad).c_str(), cantidadLegible((uint64_t)caso.distintos.promedio()).c_str());
                imprimirCaso(caso);
                resumen.agregar(nombreDistribucion(distribucion), cantidadLegible((uint64_t)lado) + "³",
                                ocupado + " de celdas", caso);
            }
        }
    }

    //PASO 2: los espacios grandes, 100.000 puntos, mismas semillas que ckdtree3d_bench
    const int64_t ladosGrandes[] = { 1000, 10000, 100000, 1000000, 1500000 };
    const uint64_t cantidad = 100000;
    for (Distribucion distribucion : distribuciones) {
        for (int64_t lado : ladosGrandes) {
            Caso caso;
            for (int repeticion = 0; repeticion < repeticiones; ++repeticion) {
                std::mt19937_64 azar((uint64_t)lado * 1000 + 10 * distribucion + repeticion);
                medirCaso(caso, generarEnCubo(distribucion, lado, cantidad, azar), lado, azar);
            }
            printf("\n%s, espacio %s, %s puntos (%s distintos)\n", nombreDistribucion(distribucion),
                   espacioLegible(lado).c_str(), cantidadLegible(cantidad).c_str(),
                   cantidadLegible((uint64_t)caso.distintos.promedio()).c_str());
            imprimirCaso(caso);
            resumen.agregar(nombreDistribucion(distribucion), cantidadLegible((uint64_t)lado) + "³",
                            cantidadLegible(cantidad), caso);
        }
    }

    printf("\nResumen 1: memoria por punto distinto (bits)\n\n");
    resumen.memoria.imprimir();
    printf("  k3 puntos y k3 niveles: variantes con límite de puntos por hoja y de niveles.\n");
    printf("  Tesis/mejor k3: lo que ocupa la tesis respecto de la variante más chica (menos de 100%% = ocupa menos).\n");
    printf("\nResumen 2: tiempo en encontrar un punto\n\n");
    resumen.punto.imprimir();
    printf("  Tesis/k3: tiempo de la tesis respecto del k3-tree básico (menos de 100%% = más rápida).\n");
    printf("\nResumen 3: consulta de rango, caja del %s del espacio\n\n", NOMBRE_CAJA[TAMANOS_CAJA - 1]);
    resumen.rango.imprimir();
    printf("  Tesis/k3: tiempo de la tesis respecto del k3-tree básico (menos de 100%% = más rápida).\n\n");
    return 0;
}
