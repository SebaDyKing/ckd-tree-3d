/**
 * Benchmark contra la versión 3D de Silva Pereira (espiral global BI sobre dacs.c de 32 bits).
 * Dos árboles con los mismos puntos y el mismo recorrido:
 *   BI + dacs.c   Silva Pereira
 *   UD + Dac64    la tesis: espiral UD y DAC de 64 bits
 * Se miden consultas de punto (que están y que no) y de rango, como en el paper [2].
 * dacs.c reserva un histograma de (código máximo + 1) enteros, así que solo
 * se puede medir en cubos chicos; se informa cuánta memoria reservó.
 *
 *   ./ckdtree3d_dacs_bench [repeticiones]     por defecto 5
 */

#include "Bench3d.hpp"

extern "C" {
#include "dacs.h"
//función interna de dacs.c, no está en dacs.h
unsigned int memoryUsage(FTRep *rep);
}

/** kd-tree 3D sobre dacs.c: los mismos códigos que CompactKdTree3d<Codec>, en el DAC de 32 bits. */
template <typename Codec>
struct KdTree3dDacOriginal {
    Point3D raiz;
    FTRep *rep = nullptr;
    uint64_t cantidad = 0;
    uint64_t codigoMaximo = 0;

    explicit KdTree3dDacOriginal(std::vector<Point3D> puntos) : cantidad(puntos.size()) {
        const std::vector<uint64_t> codigos = kdtree3dEncodeTree<Codec>(puntos, raiz);
        std::vector<unsigned int> codigos32(codigos.size());
        for (size_t i = 0; i < codigos.size(); ++i) {
            if (codigos[i] > 0xFFFFFFFFULL) {
                fprintf(stderr, "\nERROR: código %llu no cabe en dacs.c\n", (unsigned long long)codigos[i]);
                exit(1);
            }
            codigos32[i] = (unsigned int)codigos[i];
            codigoMaximo = std::max(codigoMaximo, codigos[i]);
        }
        rep = createFT(codigos32.data(), (unsigned int)codigos32.size());
    }
    ~KdTree3dDacOriginal() { destroyFT(rep); }
    KdTree3dDacOriginal(const KdTree3dDacOriginal &) = delete;
    KdTree3dDacOriginal &operator=(const KdTree3dDacOriginal &) = delete;

    uint64_t leer(uint64_t indice) const { return accessFT(rep, (unsigned int)(indice + 1)); }

    bool contains(const Point3D &q) const {
        const auto nodeAt = [this](uint64_t indice, const Point3D &padre, Spiral3dAxis eje, bool izquierdo) {
            return Codec::decode(leer(indice), padre, eje, izquierdo);
        };
        return kdtree3dSearch(nodeAt, q, raiz, 0, (int64_t)cantidad - 1, 0, true);
    }
    uint64_t countInRange(const Box3D &caja) const {
        const auto nodeAt = [this](uint64_t indice, const Point3D &padre, Spiral3dAxis eje, bool izquierdo) {
            return Codec::decode(leer(indice), padre, eje, izquierdo);
        };
        return kdtree3dRangeCount(nodeAt, caja, raiz, 0, (int64_t)cantidad - 1, 0, true, Box3D::everything());
    }
    uint64_t sizeInBits() const { return 8ULL * memoryUsage(rep) + 3 * 64; }
};

typedef KdTree3dDacOriginal<Spiral3dBICodec> KdTree3dSilva;

/** Mide un árbol completo y su desglose; leer(arbol, i) saca el código i de su DAC. */
template <typename Arbol, typename Codec, typename Leer>
static void medirTodo(const std::vector<Point3D> &puntos, Consultas &consultas,
                      Resultado &resultado, const char *nombre, const Leer &leer) {
    medir<Arbol>(puntos, consultas, resultado, nombre);
    const Arbol arbol(puntos);
    medirDesglose<Codec>(puntos, [&](uint64_t i) { return leer(arbol, i); }, consultas.presentes, resultado);
}

int main(int argc, char **argv) {
    const int repeticiones = argc > 1 ? atoi(argv[1]) : 5;

    printf("\nkd-tree 3D: Silva Pereira (BI + dacs.c de 32 bits) contra la tesis (UD + Dac64)\n");
    printf("Cada espacio es un cubo de lado L con una fracción de sus celdas ocupadas.\n");
    imprimirLeyenda(repeticiones);

    //dos tablas angostas para el resumen: tiempo y memoria
    Tabla tiempo, memoria, rango;
    tiempo.encabezados  = { "Datos", "Espacio", "Ocupado", "Silva (µs)", "Tesis (µs)", "Ahorro (ms)" };
    memoria.encabezados = { "Datos", "Espacio", "Ocupado", "Silva (bits)", "Tesis (bits)", "Tesis/Silva" };
    rango.encabezados   = { "Datos", "Espacio", "Ocupado", "Silva (µs)", "Tesis (µs)", "Tesis/Silva" };

    const auto leerDacOriginal = [](const KdTree3dSilva &arbol, uint64_t i) { return arbol.leer(i); };
    const auto leerDac64 = [](const cKDtree3dUD &arbol, uint64_t i) { return dac64Access(arbol.codes(), i); };

    //dacs.c no aguanta cubos mucho más grandes: su histograma crece con el código más grande
    const int64_t lados[] = { 64, 128, 256 };
    const double ocupaciones[] = { 0.01, 0.10 };
    const Distribucion distribuciones[] = { UNIFORME, GAUSSIANA };
    for (Distribucion distribucion : distribuciones) {
        for (int64_t lado : lados) {
            for (double ocupacion : ocupaciones) {
                const uint64_t cantidad = (uint64_t)(ocupacion * (double)lado * (double)lado * (double)lado);
                Resultado silva, tesis;
                uint64_t codigoMaximo = 0;
                Consultas consultas;
                for (int repeticion = 0; repeticion < repeticiones; ++repeticion) {
                    std::mt19937_64 azar((uint64_t)lado * 1000 + (uint64_t)(ocupacion * 100) * 10
                                         + 100000 * distribucion + repeticion);
                    const std::vector<Point3D> puntos = generarEnCubo(distribucion, lado, cantidad, azar);
                    consultas = generarConsultas(puntos, lado, azar);

                    medirTodo<KdTree3dSilva, Spiral3dBICodec>(puntos, consultas, silva, "BI + dacs.c", leerDacOriginal);
                    medirTodo<cKDtree3dUD, Spiral3dUDCodec>(puntos, consultas, tesis, "UD + Dac64", leerDac64);
                    codigoMaximo = std::max(codigoMaximo, KdTree3dSilva(puntos).codigoMaximo);
                }

                const std::string ocupado = numero(100.0 * ocupacion, 0) + "%";
                printf("\n%s, espacio %s, %s de las celdas con un punto (%s puntos)\n",
                       nombreDistribucion(distribucion), espacioLegible(lado).c_str(), ocupado.c_str(),
                       cantidadLegible(cantidad).c_str());
                printf("  código más grande de Silva %s: dacs.c reserva %.0f MB solo para contarlos\n",
                       cantidadLegible(codigoMaximo).c_str(), 4.0 * (double)(codigoMaximo + 1) / (1024.0 * 1024.0));
                Tabla tabla = tablaDeResultados();
                tabla.filas.push_back(filaDeResultado("Silva (BI + dacs.c)", silva));
                tabla.filas.push_back(filaDeResultado("Tesis (UD + Dac64)", tesis));
                tabla.imprimir();
                Tabla tablaRango = tablaDeRango();
                tablaRango.filas.push_back(filaDeRango("Silva (BI + dacs.c)", silva));
                tablaRango.filas.push_back(filaDeRango("Tesis (UD + Dac64)", tesis));
                tablaRango.imprimir();
                imprimirPuntosPorCaja(consultas);
                Tabla desglose = tablaDeDesglose();
                desglose.filas.push_back(filaDeDesglose("Silva (BI + dacs.c)", silva));
                desglose.filas.push_back(filaDeDesglose("Tesis (UD + Dac64)", tesis));
                desglose.imprimir();

                const std::string espacio = cantidadLegible((uint64_t)lado) + "³";
                const double ahorroMs = (silva.encontrarNs.promedio() - tesis.encontrarNs.promedio()) * CONSULTAS * 1e-6;
                tiempo.filas.push_back({ nombreDistribucion(distribucion), espacio, ocupado,
                                         numero(silva.encontrarNs.promedio() * 1e-3, 2),
                                         numero(tesis.encontrarNs.promedio() * 1e-3, 2), numero(ahorroMs, 1) });
                memoria.filas.push_back({ nombreDistribucion(distribucion), espacio, ocupado,
                                          numero(silva.bitsPorPunto.promedio(), 2),
                                          numero(tesis.bitsPorPunto.promedio(), 2),
                                          numero(100.0 * tesis.bitsPorPunto.promedio()
                                                 / silva.bitsPorPunto.promedio(), 1) + "%" });
                const int cajaGrande = TAMANOS_CAJA - 1;
                rango.filas.push_back({ nombreDistribucion(distribucion), espacio, ocupado,
                                        numero(silva.rangoUs[cajaGrande].promedio(), 1),
                                        numero(tesis.rangoUs[cajaGrande].promedio(), 1),
                                        numero(100.0 * tesis.rangoUs[cajaGrande].promedio()
                                               / silva.rangoUs[cajaGrande].promedio(), 1) + "%" });
            }
        }
    }

    printf("\nResumen 1: tiempo en encontrar un punto\n\n");
    tiempo.imprimir();
    printf("  Ahorro: milisegundos que ahorra la tesis respecto de Silva en %s búsquedas.\n",
           cantidadLegible(CONSULTAS).c_str());
    printf("\nResumen 2: memoria por punto\n\n");
    memoria.imprimir();
    printf("  Espacio L³: cubo de lado L. Ocupado: porcentaje de celdas con un punto.\n");
    printf("  Tesis/Silva: lo que ocupa la tesis respecto de Silva (menos de 100%% = ocupa menos).\n");
    printf("  En cubos de lado mayor a ~2048 Silva no funciona: sus códigos no caben en 32 bits.\n");
    printf("\nResumen 3: consulta de rango, caja del %s del espacio\n\n", NOMBRE_CAJA[TAMANOS_CAJA - 1]);
    rango.imprimir();
    printf("  Tesis/Silva: tiempo de la tesis respecto de Silva (menos de 100%% = más rápida).\n\n");
    return 0;
}
