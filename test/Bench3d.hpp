#ifndef BENCH3D_HPP
#define BENCH3D_HPP

/* Mediciones compartidas por los benchmarks del kd-tree 3D, con la metodología
   del paper del cKd-tree [2]: buscar puntos que están y que no están, y contar
   los puntos de cajas del 0,01 %, 0,1 % y 1 % del espacio. Incluye el desglose
   de encontrar un punto (leer el código del DAC, decode y el resto), el encode,
   promedio y desviación, y las tablas. */

#include "cKDtree3d.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

static const uint64_t CONSULTAS = 100000;          //búsquedas que están, y otras tantas que no, por repetición
static const uint64_t CONSULTAS_DESGLOSE = 20000;  //búsquedas que se graban paso a paso
static const uint64_t CAJAS = 1000;                //cajas por tamaño, como en el paper

/** Tamaños de caja de la consulta de rango: fracción del volumen del espacio. */
static const int TAMANOS_CAJA = 3;
static const double FRACCION_CAJA[TAMANOS_CAJA] = { 0.0001, 0.001, 0.01 };
static const char *const NOMBRE_CAJA[TAMANOS_CAJA] = { "0,01 %", "0,1 %", "1 %" };

/** Promedio y desviación estándar muestral. */
struct Estadistica {
    std::vector<double> muestras;
    void agregar(double valor) { muestras.push_back(valor); }
    double promedio() const {
        double suma = 0;
        for (double m : muestras) suma += m;
        return suma / (double)muestras.size();
    }
    double desviacion() const {
        if (muestras.size() < 2) return 0;
        const double media = promedio();
        double suma = 0;
        for (double m : muestras) suma += (m - media) * (m - media);
        return std::sqrt(suma / (double)(muestras.size() - 1));
    }
};

struct Resultado {
    Estadistica bitsPorPunto, construccionMs, encontrarNs, descartarNs;  //árbol completo
    Estadistica rangoUs[TAMANOS_CAJA];                                   //consulta de rango
    Estadistica pasos, leerNs, decodeNs, encodeMs;                       //desglose
};

static bool menorXYZ(const Point3D &a, const Point3D &b) {
    return kdtree3dLess(a, b, SPIRAL3D_AXIS_X);
}

/**
 * Consultas de una repetición, iguales para todas las estructuras.
 * La primera estructura medida fija cuántos puntos tiene cada caja; las demás deben coincidir.
 */
struct Consultas {
    std::vector<Point3D> presentes, ausentes;
    std::vector<Box3D> cajas[TAMANOS_CAJA];
    std::vector<uint64_t> dentro[TAMANOS_CAJA];
};

/**
 * Puntos guardados elegidos al azar, puntos del cubo que no están y cajas de cada tamaño.
 * @param lado lado del cubo donde están los puntos, coordenadas de 0 a lado - 1
 */
static Consultas generarConsultas(const std::vector<Point3D> &puntos, int64_t lado, std::mt19937_64 &azar) {
    Consultas consultas;
    consultas.presentes.resize(CONSULTAS);
    for (Point3D &q : consultas.presentes) q = puntos[azar() % puntos.size()];

    std::vector<Point3D> ordenados(puntos);
    std::sort(ordenados.begin(), ordenados.end(), menorXYZ);
    std::uniform_int_distribution<int64_t> coordenada(0, lado - 1);
    while (consultas.ausentes.size() < CONSULTAS) {
        const Point3D q(coordenada(azar), coordenada(azar), coordenada(azar));
        if (!std::binary_search(ordenados.begin(), ordenados.end(), q, menorXYZ)) consultas.ausentes.push_back(q);
    }

    //una caja cúbica con la fracción pedida del volumen, en una posición al azar del cubo
    for (int t = 0; t < TAMANOS_CAJA; ++t) {
        const int64_t ladoCaja = std::max<int64_t>(1, std::llround((double)lado * std::cbrt(FRACCION_CAJA[t])));
        std::uniform_int_distribution<int64_t> inicio(0, lado - ladoCaja);
        for (uint64_t i = 0; i < CAJAS; ++i) {
            const Point3D bajo(inicio(azar), inicio(azar), inicio(azar));
            const Point3D alto(bajo.x + ladoCaja - 1, bajo.y + ladoCaja - 1, bajo.z + ladoCaja - 1);
            consultas.cajas[t].push_back(Box3D(bajo, alto));
        }
    }
    return consultas;
}

static double nanosegundosDesde(std::chrono::steady_clock::time_point inicio) {
    return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - inicio).count();
}

[[noreturn]] static void abortarMedicion(const char *nombre, const char *motivo) {
    fprintf(stderr, "\nERROR: %s %s\n", nombre, motivo);
    exit(1);
}

/**
 * Construye y mide: memoria, construcción, buscar puntos que están y que no, y rango.
 * Aborta si alguna respuesta está mal, así ningún tiempo sale de una estructura con errores.
 */
template <typename Arbol>
static void medir(const std::vector<Point3D> &puntos, Consultas &consultas, Resultado &resultado, const char *nombre) {
    const auto inicio = std::chrono::steady_clock::now();
    const Arbol arbol(puntos);
    resultado.construccionMs.agregar(nanosegundosDesde(inicio) * 1e-6);
    resultado.bitsPorPunto.agregar((double)arbol.sizeInBits() / (double)puntos.size());

    uint64_t encontrados = 0;
    auto t = std::chrono::steady_clock::now();
    for (const Point3D &q : consultas.presentes) encontrados += arbol.contains(q);
    resultado.encontrarNs.agregar(nanosegundosDesde(t) / CONSULTAS);
    if (encontrados != CONSULTAS) abortarMedicion(nombre, "no encontró todos los puntos guardados");

    t = std::chrono::steady_clock::now();
    for (const Point3D &q : consultas.ausentes) encontrados += arbol.contains(q);
    resultado.descartarNs.agregar(nanosegundosDesde(t) / CONSULTAS);
    if (encontrados != CONSULTAS) abortarMedicion(nombre, "encontró puntos que no están");

    for (int c = 0; c < TAMANOS_CAJA; ++c) {
        std::vector<uint64_t> dentro(CAJAS);
        t = std::chrono::steady_clock::now();
        for (uint64_t i = 0; i < CAJAS; ++i) dentro[i] = arbol.countInRange(consultas.cajas[c][i]);
        resultado.rangoUs[c].agregar(nanosegundosDesde(t) * 1e-3 / CAJAS);
        if (consultas.dentro[c].empty()) consultas.dentro[c] = dentro;
        else if (consultas.dentro[c] != dentro) abortarMedicion(nombre, "contó distinto que las otras estructuras");
    }
}

/** Promedio de puntos dentro de cada caja de un tamaño, según la primera estructura medida. */
static double puntosPorCaja(const Consultas &consultas, int tamano) {
    double suma = 0;
    for (uint64_t d : consultas.dentro[tamano]) suma += (double)d;
    return suma / (double)CAJAS;
}

/** Un paso del recorrido: qué nodo se leyó y con qué datos se decodificó. */
struct Paso {
    uint64_t indice, codigo;
    Point3D padre;
    Spiral3dAxis eje;
    bool izquierdo;
};

/**
 * Desglose de encontrar un punto: se graban los pasos de CONSULTAS_DESGLOSE búsquedas
 * y después se mide por separado leer esos códigos del DAC y decodificarlos.
 * También mide el encode: generar los códigos de todos los puntos.
 * @param leer leer(indice) devuelve el código guardado en el DAC del árbol
 */
template <typename Codec, typename Leer>
static void medirDesglose(std::vector<Point3D> puntos, const Leer &leer,
                          const std::vector<Point3D> &busquedas, Resultado &resultado) {
    kdtree3dBuildOrder(puntos, 0, (int64_t)puntos.size() - 1, 0);
    const Point3D raiz = puntos[((int64_t)puntos.size() - 1) / 2];
    std::vector<uint64_t> codigos(puntos.size(), 0);
    auto t = std::chrono::steady_clock::now();
    kdtree3dGenerateCodes<Codec>(puntos, codigos, raiz, 0, (int64_t)puntos.size() - 1, 0);
    resultado.encodeMs.agregar(nanosegundosDesde(t) * 1e-6);

    std::vector<Paso> pasos;
    const auto grabar = [&](uint64_t indice, const Point3D &padre, Spiral3dAxis eje, bool izquierdo) {
        const uint64_t codigo = leer(indice);
        pasos.push_back({ indice, codigo, padre, eje, izquierdo });
        return Codec::decode(codigo, padre, eje, izquierdo);
    };
    for (uint64_t i = 0; i < CONSULTAS_DESGLOSE; ++i)
        kdtree3dSearch(grabar, busquedas[i], raiz, 0, (int64_t)puntos.size() - 1, 0, true);
    resultado.pasos.agregar((double)pasos.size() / CONSULTAS_DESGLOSE);

    uint64_t control = 0;
    t = std::chrono::steady_clock::now();
    for (const Paso &p : pasos) control += leer(p.indice);
    resultado.leerNs.agregar(nanosegundosDesde(t) / CONSULTAS_DESGLOSE);

    t = std::chrono::steady_clock::now();
    for (const Paso &p : pasos) control += (uint64_t)Codec::decode(p.codigo, p.padre, p.eje, p.izquierdo).x;
    resultado.decodeNs.agregar(nanosegundosDesde(t) / CONSULTAS_DESGLOSE);
    if (control == 42) printf(" "); //evita que el compilador descarte los ciclos medidos
}

//caracteres visibles de un texto UTF-8: los acentos ocupan dos bytes y una columna
static int anchoVisible(const std::string &texto) {
    int ancho = 0;
    for (unsigned char c : texto) if ((c & 0xC0) != 0x80) ++ancho;
    return ancho;
}

/** Tabla con bordes; cada columna toma el ancho de su texto más largo. */
struct Tabla {
    std::vector<std::string> encabezados;
    std::vector<std::vector<std::string>> filas;

    void imprimir() const {
        std::vector<int> anchos;
        for (const std::string &e : encabezados) anchos.push_back(anchoVisible(e));
        for (const std::vector<std::string> &fila : filas)
            for (size_t c = 0; c < fila.size(); ++c) anchos[c] = std::max(anchos[c], anchoVisible(fila[c]));

        const auto linea = [&anchos](const char *izquierda, const char *medio, const char *derecha) {
            printf("  %s", izquierda);
            for (size_t c = 0; c < anchos.size(); ++c) {
                for (int i = 0; i < anchos[c] + 2; ++i) printf("─");
                printf("%s", c + 1 < anchos.size() ? medio : derecha);
            }
            printf("\n");
        };
        //la primera columna va a la izquierda y los números a la derecha
        const auto fila = [&anchos](const std::vector<std::string> &celdas) {
            printf("  │");
            for (size_t c = 0; c < celdas.size(); ++c) {
                const int relleno = anchos[c] - anchoVisible(celdas[c]);
                if (c == 0) printf(" %s%*s │", celdas[c].c_str(), relleno, "");
                else        printf(" %*s%s │", relleno, "", celdas[c].c_str());
            }
            printf("\n");
        };

        linea("┌", "┬", "┐");
        fila(encabezados);
        linea("├", "┼", "┤");
        for (const std::vector<std::string> &f : filas) fila(f);
        linea("└", "┴", "┘");
    }
};

/** "14.97 ± 0.05" con los decimales pedidos. */
static std::string conDesviacion(const Estadistica &e, int decimales, double escala = 1.0) {
    char texto[64];
    snprintf(texto, sizeof texto, "%.*f ± %.*f", decimales, e.promedio() * escala, decimales, e.desviacion() * escala);
    return texto;
}

static std::string numero(double valor, int decimales) {
    char texto[64];
    snprintf(texto, sizeof texto, "%.*f", decimales, valor);
    return texto;
}

/** Puntos con separador de miles: 1.000.000. */
static std::string cantidadLegible(uint64_t cantidad) {
    std::string texto = std::to_string(cantidad);
    for (int i = (int)texto.size() - 3; i > 0; i -= 3) texto.insert((size_t)i, ".");
    return texto;
}

enum Distribucion { UNIFORME, GAUSSIANA };

static const char *nombreDistribucion(Distribucion distribucion) {
    return distribucion == UNIFORME ? "uniforme" : "gaussiana";
}

/**
 * Puntos dentro de un cubo de lado L, con coordenadas de 0 a L - 1.
 * Uniforme: cualquier celda con la misma probabilidad.
 * Gaussiana: concentrados en el centro, desviación L / 8 por eje, recortados al cubo.
 */
static std::vector<Point3D> generarEnCubo(Distribucion distribucion, int64_t lado, uint64_t cantidad,
                                          std::mt19937_64 &azar) {
    std::vector<Point3D> puntos(cantidad);
    if (distribucion == UNIFORME) {
        std::uniform_int_distribution<int64_t> coordenada(0, lado - 1);
        for (Point3D &p : puntos) p = Point3D(coordenada(azar), coordenada(azar), coordenada(azar));
        return puntos;
    }
    std::normal_distribution<double> normal(lado / 2.0, lado / 8.0);
    const auto coordenada = [&]() {
        return std::min<int64_t>(lado - 1, std::max<int64_t>(0, std::llround(normal(azar))));
    };
    for (Point3D &p : puntos) p = Point3D(coordenada(), coordenada(), coordenada());
    return puntos;
}

/** "1.000 × 1.000 × 1.000" */
static std::string espacioLegible(int64_t lado) {
    const std::string l = cantidadLegible((uint64_t)lado);
    return l + " × " + l + " × " + l;
}

/** Porcentaje de celdas del cubo que tienen un punto (antes de quitar repetidos). */
static std::string densidadLegible(int64_t lado, uint64_t cantidad) {
    const double densidad = 100.0 * (double)cantidad / ((double)lado * (double)lado * (double)lado);
    char texto[32];
    if (densidad >= 0.01) snprintf(texto, sizeof texto, "%.2f%%", densidad);
    else                  snprintf(texto, sizeof texto, "%.0e%%", densidad);
    return texto;
}

static void imprimirLeyenda(int repeticiones) {
    printf("\nQué significa cada columna:\n");
    printf("  Bits/punto      memoria que ocupa cada punto guardado (menos es mejor)\n");
    printf("  Construir (ms)  milisegundos en armar el árbol con todos los puntos\n");
    printf("  Está (µs)       microsegundos promedio en encontrar un punto guardado\n");
    printf("  No está (µs)    microsegundos promedio en confirmar que un punto no está\n");
    printf("  Rango (µs)      microsegundos promedio en contar los puntos de una caja cúbica\n");
    printf("                  que ocupa el 0,01 %%, el 0,1 %% o el 1 %% del espacio (%s cajas)\n",
           cantidadLegible(CAJAS).c_str());
    printf("  Desglose de Está, por búsqueda:\n");
    printf("    Pasos         nodos que baja el árbol hasta dar con el punto\n");
    printf("    Leer (µs)     sacar el código de cada nodo del DAC\n");
    printf("    Decode (µs)   pasar cada código a coordenadas con la espiral\n");
    printf("    Resto (µs)    comparar coordenadas y bajar (Encontrar - Leer - Decode)\n");
    printf("  Encode (ms)     generar los códigos de todos los puntos al construir\n");
    printf("  Cada valor es el promedio de %d repeticiones con datos nuevos; ± es cuánto\n", repeticiones);
    printf("  varía entre ellas. Cada repetición busca %s puntos que están y %s que no.\n",
           cantidadLegible(CONSULTAS).c_str(), cantidadLegible(CONSULTAS).c_str());
}

static Tabla tablaDeResultados() {
    Tabla tabla;
    tabla.encabezados = { "Estructura", "Bits/punto", "Construir (ms)", "Está (µs)", "No está (µs)" };
    return tabla;
}

static std::vector<std::string> filaDeResultado(const char *nombre, const Resultado &r) {
    return { nombre, conDesviacion(r.bitsPorPunto, 2), conDesviacion(r.construccionMs, 1),
             conDesviacion(r.encontrarNs, 2, 1e-3), conDesviacion(r.descartarNs, 2, 1e-3) };
}

static Tabla tablaDeRango() {
    Tabla tabla;
    tabla.encabezados = { "Estructura" };
    for (int c = 0; c < TAMANOS_CAJA; ++c) tabla.encabezados.push_back(std::string("Rango ") + NOMBRE_CAJA[c] + " (µs)");
    return tabla;
}

static std::vector<std::string> filaDeRango(const char *nombre, const Resultado &r) {
    std::vector<std::string> fila = { nombre };
    for (int c = 0; c < TAMANOS_CAJA; ++c) fila.push_back(conDesviacion(r.rangoUs[c], 1));
    return fila;
}

/** "  Puntos por caja: 0,01 % ≈ 1, 0,1 % ≈ 10, 1 % ≈ 100" */
static void imprimirPuntosPorCaja(const Consultas &consultas) {
    printf("  Puntos por caja:");
    for (int c = 0; c < TAMANOS_CAJA; ++c)
        printf("%s %s ≈ %s", c ? "," : "", NOMBRE_CAJA[c], numero(puntosPorCaja(consultas, c), 0).c_str());
    printf("\n");
}

static Tabla tablaDeDesglose() {
    Tabla tabla;
    tabla.encabezados = { "Estructura", "Pasos", "Leer (µs)", "Decode (µs)", "Resto (µs)", "Encode (ms)" };
    return tabla;
}

static std::vector<std::string> filaDeDesglose(const char *nombre, const Resultado &r) {
    const double resto = r.encontrarNs.promedio() - r.leerNs.promedio() - r.decodeNs.promedio();
    return { nombre, numero(r.pasos.promedio(), 1), numero(r.leerNs.promedio() * 1e-3, 2),
             numero(r.decodeNs.promedio() * 1e-3, 2), numero(std::max(resto, 0.0) * 1e-3, 2),
             numero(r.encodeMs.promedio(), 1) };
}

#endif
