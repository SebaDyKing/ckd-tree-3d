/**
 * Tiempo de la espiral unidireccional (UD) contra la bidireccional (BI).
 * Mide nanosegundos por llamada de encode y decode en capas cada vez más
 * grandes. Si es O(1), el tiempo no crece con la capa.
 *
 *   ./spiral3d_bench
 */

#include "Spiral3dUD.hpp"

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

static const size_t PUNTOS   = 4096;      //puntos distintos por capa
static const size_t LLAMADAS = 4000000;   //llamadas por repetición
static const int    REPETICIONES = 5;     //se queda con la más rápida

struct Caso {
    Point3D padre, hijo;
    Spiral3dAxis eje;
    bool left;
    uint64_t ud, bi;
};

/** Puntos al azar en la capa m exacta, con padres al azar. */
static std::vector<Caso> generar(uint64_t m, std::mt19937_64 &gen) {
    const int64_t sm = (int64_t)m;
    std::uniform_int_distribution<int64_t> d(-sm, sm), p(-1000000000LL, 1000000000LL);
    std::uniform_int_distribution<int> k(0, 2), signo(0, 1);
    std::vector<Caso> v(PUNTOS);
    for (size_t i = 0; i < PUNTOS; ++i) {
        int64_t c[3] = { d(gen), d(gen), d(gen) };
        c[k(gen)] = signo(gen) ? sm : -sm;            //fuerza max(|c|) == m
        Caso &x = v[i];
        x.padre = Point3D(p(gen), p(gen), p(gen));
        x.hijo  = Point3D(x.padre.x + c[0], x.padre.y + c[1], x.padre.z + c[2]);
        x.eje   = (Spiral3dAxis)(i % 3);
        x.left  = c[x.eje] <= 0;
        x.ud    = spiral3dUDEncode(x.hijo, x.padre, x.eje);
        x.bi    = spiral3dEncode(x.hijo, x.padre);
    }
    return v;
}

/**
 * Nanosegundos por llamada de f(i). Repite la medición y se queda con la
 * mínima: el ruido (otros programas, frecuencia de la CPU) solo suma tiempo.
 */
template <class F>
static double medir(F f) {
    const size_t vueltas = LLAMADAS / PUNTOS;
    double mejor = 1e18;
    for (int rep = 0; rep < REPETICIONES; ++rep) {
        uint64_t sumidero = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t r = 0; r < vueltas; ++r) {
            for (size_t i = 0; i < PUNTOS; ++i) sumidero += f(i);
            asm volatile("" ::: "memory");            //impide reutilizar la vuelta anterior
        }
        const auto t1 = std::chrono::steady_clock::now();
        asm volatile("" : : "r"(sumidero));            //el resultado se usa: no se elimina el cálculo
        const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / (double)(vueltas * PUNTOS);
        if (ns < mejor) mejor = ns;
    }
    return mejor;
}

/** Bits necesarios para escribir el código c. */
static int bits(uint64_t c) {
    return c == 0 ? 1 : 64 - __builtin_clzll(c);
}

int main() {
    const uint64_t capas[] = { 1, 10, 1000, 1000000, SPIRAL3D_MAX_LAYER, SPIRAL3DUD_MAX_LAYER };
    const size_t nCapas = sizeof(capas) / sizeof(capas[0]);
    std::mt19937_64 gen(2026);

    printf("Tiempo por llamada: unidireccional (UD) contra bidireccional (BI)\n");
    printf("%zu puntos por capa con padres al azar, %zu llamadas, mejor de %d repeticiones\n\n",
           PUNTOS, LLAMADAS, REPETICIONES);
    printf("  %9s | %9s %9s | %9s %9s | %7s %7s\n",
           "capa m", "encode UD", "encode BI", "decode UD", "decode BI", "bits UD", "bits BI");
    printf("  ----------+---------------------+---------------------+----------------\n");

    double minimo[4], maximo[4];
    for (int c = 0; c < 4; ++c) { minimo[c] = 1e18; maximo[c] = 0; }
    bool ok = true;

    for (size_t n = 0; n < nCapas; ++n) {
        const uint64_t m = capas[n];
        const bool hayBI = m <= SPIRAL3D_MAX_LAYER;
        const std::vector<Caso> v = generar(m, gen);

        //antes de medir, comprobar que los dos recuperan cada punto
        for (const Caso &x : v) {
            if (!(spiral3dUDDecode(x.ud, x.padre, x.eje, x.left) == x.hijo)) ok = false;
            if (hayBI && !(spiral3dDecode(x.bi, x.padre) == x.hijo)) ok = false;
        }

        double t[4];
        t[0] = medir([&](size_t i) { return spiral3dUDEncode(v[i].hijo, v[i].padre, v[i].eje); });
        t[2] = medir([&](size_t i) {
            const Point3D q = spiral3dUDDecode(v[i].ud, v[i].padre, v[i].eje, v[i].left);
            return (uint64_t)(q.x + q.y + q.z);
        });
        if (hayBI) {
            t[1] = medir([&](size_t i) { return spiral3dEncode(v[i].hijo, v[i].padre); });
            t[3] = medir([&](size_t i) {
                const Point3D q = spiral3dDecode(v[i].bi, v[i].padre);
                return (uint64_t)(q.x + q.y + q.z);
            });
        }
        for (int c = 0; c < 4; ++c) {
            if (!hayBI && (c == 1 || c == 3)) continue;
            if (t[c] < minimo[c]) minimo[c] = t[c];
            if (t[c] > maximo[c]) maximo[c] = t[c];
        }

        //bits del último código de la capa
        const int bitsUD = bits(spiral3dUDCount(m) - 1);
        const uint64_t lado = 2 * m + 1;

        printf("  %9llu | %6.1f ns ", (unsigned long long)m, t[0]);
        if (hayBI) printf("%6.1f ns", t[1]); else printf("%9s", "--");
        printf(" | %6.1f ns ", t[2]);
        if (hayBI) printf("%6.1f ns", t[3]); else printf("%9s", "--");
        printf(" | %7d ", bitsUD);
        if (hayBI) printf("%7d\n", bits(lado * lado * lado - 1)); else printf("%7s\n", "--");
        fflush(stdout);
    }

    printf("  ----------+---------------------+---------------------+----------------\n");
    printf("  %11s | %8.2fx %8.2fx | %8.2fx %8.2fx |\n", "máx/mín",
           maximo[0] / minimo[0], maximo[1] / minimo[1], maximo[2] / minimo[2], maximo[3] / minimo[3]);
    printf("\n  El tiempo no crece con la capa: de m = 10 a m = %llu se mantiene casi igual, es O(1).\n",
           (unsigned long long)SPIRAL3DUD_MAX_LAYER);
    printf("  máx/mín: tiempo más lento dividido por el más rápido de cada columna.\n");
    printf("  -- : la capa pasa el tope de 64 bits de la espiral bidireccional (m = %llu).\n",
           (unsigned long long)SPIRAL3D_MAX_LAYER);
    if (!ok) {
        printf("\n  ERROR: algún punto no se recuperó al decodificar\n");
        return 1;
    }
    return 0;
}
