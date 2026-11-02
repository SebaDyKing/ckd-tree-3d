/**
 * Benchmark del kd-tree 3D compacto: espiral UD contra BI, ambas sobre el Dac64,
 * y el iKd-tree (implícito, sin comprimir: 3 x 64 bits por punto) con el mismo recorrido.
 * Cada caso es un espacio (cubo de lado L) con puntos uniformes o gaussianos;
 * se repite con datos nuevos y se informa promedio y desviación estándar.
 * Se miden consultas de punto (que están y que no) y de rango, como en el paper [2].
 *
 *   ./ckdtree3d_bench [repeticiones] [puntos]     por defecto 5 y 100000
 */

#include "Bench3d.hpp"

/** iKd-tree: kd-tree implícito sin comprimir; los mismos nodos y el mismo recorrido, sin decodificar. */
struct IKdTree3d {
    std::vector<Point3D> nodos;

    explicit IKdTree3d(std::vector<Point3D> puntos) : nodos(std::move(puntos)) {
        if (!nodos.empty()) kdtree3dBuildOrder(nodos, 0, (int64_t)nodos.size() - 1, 0);
    }
    bool contains(const Point3D &q) const {
        if (nodos.empty()) return false;
        const std::vector<Point3D> &n = nodos;
        const auto nodeAt = [&n](uint64_t indice, const Point3D &, Spiral3dAxis, bool) { return n[indice]; };
        const int64_t ultimo = (int64_t)nodos.size() - 1;
        return kdtree3dSearch(nodeAt, q, nodos[ultimo / 2], 0, ultimo, 0, true);
    }
    uint64_t countInRange(const Box3D &caja) const {
        if (nodos.empty()) return 0;
        const std::vector<Point3D> &n = nodos;
        const auto nodeAt = [&n](uint64_t indice, const Point3D &, Spiral3dAxis, bool) { return n[indice]; };
        const int64_t ultimo = (int64_t)nodos.size() - 1;
        return kdtree3dRangeCount(nodeAt, caja, nodos[ultimo / 2], 0, ultimo, 0, true, Box3D::everything());
    }
    uint64_t sizeInBits() const { return nodos.size() * 3 * 64; }
};

/** "1.000³": el espacio en una columna angosta. */
static std::string espacioCorto(int64_t lado) {
    return cantidadLegible((uint64_t)lado) + "³";
}

int main(int argc, char **argv) {
    const int repeticiones = argc > 1 ? atoi(argv[1]) : 5;
    const uint64_t cantidad = argc > 2 ? strtoull(argv[2], nullptr, 10) : 100000;

    printf("\nkd-tree 3D compacto: espiral UD contra BI (las dos sobre el Dac64) y el iKd-tree (sin comprimir)\n");
    printf("Cada espacio es un cubo de lado L con %s puntos.\n", cantidadLegible(cantidad).c_str());
    imprimirLeyenda(repeticiones);

    //el último lado pasa el tope de BI (1.321.122 por eje) y no el de UD (1.664.509)
    const int64_t lados[] = { 1000, 10000, 100000, 1000000, 1500000 };

    Tabla tiempo, memoria, rango;
    tiempo.encabezados  = { "Datos", "Espacio", "UD (µs)", "BI (µs)", "Más rápido" };
    memoria.encabezados = { "Datos", "Espacio", "UD (bits)", "BI (bits)", "UD ahorra" };
    rango.encabezados   = { "Datos", "Espacio", "UD (µs)", "BI (µs)", "iKd-tree (µs)" };

    const Distribucion distribuciones[] = { UNIFORME, GAUSSIANA };
    for (Distribucion distribucion : distribuciones) {
        for (int64_t lado : lados) {
            Resultado ud, bi, ikd;
            bool biCabe = true;
            Consultas consultas;
            for (int repeticion = 0; repeticion < repeticiones; ++repeticion) {
                std::mt19937_64 azar((uint64_t)lado * 1000 + 10 * distribucion + repeticion);
                const std::vector<Point3D> puntos = generarEnCubo(distribucion, lado, cantidad, azar);
                consultas = generarConsultas(puntos, lado, azar);

                medir<cKDtree3dUD>(puntos, consultas, ud, "UD + Dac64");
                {
                    const cKDtree3dUD arbol(puntos);
                    medirDesglose<Spiral3dUDCodec>(puntos, [&arbol](uint64_t i) { return dac64Access(arbol.codes(), i); },
                                                   consultas.presentes, ud);
                }
                try {
                    medir<cKDtree3dBI>(puntos, consultas, bi, "BI + Dac64");
                    const cKDtree3dBI arbol(puntos);
                    medirDesglose<Spiral3dBICodec>(puntos, [&arbol](uint64_t i) { return dac64Access(arbol.codes(), i); },
                                                   consultas.presentes, bi);
                } catch (const std::invalid_argument &) { biCabe = false; }
                medir<IKdTree3d>(puntos, consultas, ikd, "iKd-tree");
            }

            printf("\n%s, espacio %s (%s de las celdas con un punto)\n", nombreDistribucion(distribucion),
                   espacioLegible(lado).c_str(), densidadLegible(lado, cantidad).c_str());
            Tabla tabla = tablaDeResultados();
            tabla.filas.push_back(filaDeResultado("UD + Dac64", ud));
            if (biCabe) tabla.filas.push_back(filaDeResultado("BI + Dac64", bi));
            else        tabla.filas.push_back({ "BI + Dac64", "no cabe", "no cabe", "no cabe", "no cabe" });
            tabla.filas.push_back(filaDeResultado("iKd-tree", ikd));
            tabla.imprimir();

            Tabla tablaRango = tablaDeRango();
            tablaRango.filas.push_back(filaDeRango("UD + Dac64", ud));
            if (biCabe) tablaRango.filas.push_back(filaDeRango("BI + Dac64", bi));
            tablaRango.filas.push_back(filaDeRango("iKd-tree", ikd));
            tablaRango.imprimir();
            imprimirPuntosPorCaja(consultas);

            Tabla desglose = tablaDeDesglose();
            desglose.filas.push_back(filaDeDesglose("UD + Dac64", ud));
            if (biCabe) desglose.filas.push_back(filaDeDesglose("BI + Dac64", bi));
            desglose.imprimir();
            if (!biCabe)
                printf("  BI no se puede construir: sus códigos no caben en 64 bits pasado m = %s\n",
                       cantidadLegible(SPIRAL3D_MAX_LAYER).c_str());

            //filas del resumen; el rango usa la caja más grande
            const int cajaGrande = TAMANOS_CAJA - 1;
            rango.filas.push_back({ nombreDistribucion(distribucion), espacioCorto(lado),
                                    numero(ud.rangoUs[cajaGrande].promedio(), 1),
                                    biCabe ? numero(bi.rangoUs[cajaGrande].promedio(), 1) : "no cabe",
                                    numero(ikd.rangoUs[cajaGrande].promedio(), 1) });
            const double tiempoUD = ud.encontrarNs.promedio() * 1e-3;
            const double bitsUD = ud.bitsPorPunto.promedio();
            if (!biCabe) {
                tiempo.filas.push_back({ nombreDistribucion(distribucion), espacioCorto(lado), numero(tiempoUD, 2),
                                         "no cabe", "solo UD" });
                memoria.filas.push_back({ nombreDistribucion(distribucion), espacioCorto(lado), numero(bitsUD, 2),
                                          "no cabe", "--" });
                continue;
            }
            const double tiempoBI = bi.encontrarNs.promedio() * 1e-3;
            const double bitsBI = bi.bitsPorPunto.promedio();
            const double diferencia = 100.0 * std::fabs(tiempoUD - tiempoBI) / std::max(tiempoUD, tiempoBI);
            const std::string masRapido = diferencia < 5.0 ? "empate"
                                        : std::string(tiempoUD < tiempoBI ? "UD " : "BI ") + numero(diferencia, 0) + "%";
            tiempo.filas.push_back({ nombreDistribucion(distribucion), espacioCorto(lado), numero(tiempoUD, 2),
                                     numero(tiempoBI, 2), masRapido });
            memoria.filas.push_back({ nombreDistribucion(distribucion), espacioCorto(lado), numero(bitsUD, 2),
                                      numero(bitsBI, 2), numero(bitsBI - bitsUD, 2) + " bits" });
        }
    }

    printf("\nResumen 1: tiempo en encontrar el mismo punto (%s puntos)\n\n", cantidadLegible(cantidad).c_str());
    tiempo.imprimir();
    printf("  empate: diferencia menor al 5%%, dentro de la variación entre repeticiones.\n");
    printf("\nResumen 2: memoria por punto (el iKd-tree ocupa siempre 192 bits)\n\n");
    memoria.imprimir();
    printf("  Espacio L³: cubo de lado L. Tope por eje entre padre e hijo: UD %s, BI %s.\n",
           cantidadLegible(SPIRAL3DUD_MAX_LAYER).c_str(), cantidadLegible(SPIRAL3D_MAX_LAYER).c_str());
    printf("\nResumen 3: consulta de rango, caja del %s del espacio\n\n", NOMBRE_CAJA[TAMANOS_CAJA - 1]);
    rango.imprimir();
    printf("  Microsegundos promedio en contar los puntos de una caja.\n\n");
    return 0;
}
