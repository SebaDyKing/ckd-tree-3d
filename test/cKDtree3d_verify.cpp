/* Pruebas del kd-tree 3D compacto (src/cKDtree3d.hpp) con la espiral UD y la BI.

     ./ckdtree3d_verify

   Cada respuesta de contains se compara contra un conjunto ordenado con los
   mismos puntos: todos los puntos guardados se encuentran y los demás no.
   La consulta de rango se compara contra contar uno por uno.
   También se prueba guardar y cargar, los bordes de cantidad y de tamaño del
   cubo, y que en cubos grandes aparecen códigos que solo caben en 64 bits. */

#include "cKDtree3d.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include <algorithm>

static int errores = 0;

static void reportar(const char *nombre, bool ok, const char *detalle = "") {
    //printf alinea por bytes, así que los acentos hay que contarlos aparte
    int ancho = 0;
    for (const unsigned char *p = (const unsigned char *)nombre; *p; ++p)
        if ((*p & 0xC0) != 0x80) ++ancho;

    printf("  [%-5s] %s", ok ? "OK" : "FALLA", nombre);
    for (int i = ancho; i < 56; ++i) putchar(' ');
    printf(" %s\n", detalle);
    if (!ok) ++errores;
}

static bool menorXYZ(const Point3D &a, const Point3D &b) {
    return kdtree3dLess(a, b, SPIRAL3D_AXIS_X);
}

//conjunto de referencia: los puntos ordenados, sin repetidos
struct Conjunto {
    std::vector<Point3D> puntos;
    explicit Conjunto(std::vector<Point3D> p) : puntos(std::move(p)) {
        std::sort(puntos.begin(), puntos.end(), menorXYZ);
        puntos.erase(std::unique(puntos.begin(), puntos.end()), puntos.end());
    }
    bool contiene(const Point3D &q) const { return std::binary_search(puntos.begin(), puntos.end(), q, menorXYZ); }
};

static std::vector<Point3D> puntosAzar(std::mt19937_64 &azar, uint64_t cantidad, int64_t lado) {
    std::uniform_int_distribution<int64_t> coordenada(0, lado - 1);
    std::vector<Point3D> puntos(cantidad);
    for (Point3D &p : puntos) p = Point3D(coordenada(azar), coordenada(azar), coordenada(azar));
    return puntos;
}

//los puntos guardados más consultas al azar en el cubo y un poco afuera
template <typename Arbol>
static bool respondeIgual(const Arbol &arbol, const Conjunto &conjunto, std::mt19937_64 &azar,
                          int64_t lado, uint64_t consultas, uint64_t &encontrados) {
    for (const Point3D &p : conjunto.puntos) if (!arbol.contains(p)) return false;
    std::uniform_int_distribution<int64_t> coordenada(-1, lado);
    encontrados = 0;
    for (uint64_t i = 0; i < consultas; ++i) {
        const Point3D q(coordenada(azar), coordenada(azar), coordenada(azar));
        const bool esperado = conjunto.contiene(q);
        if (arbol.contains(q) != esperado) return false;
        if (esperado) ++encontrados;
    }
    return true;
}

template <typename Arbol>
static void probarArbol(const char *nombreEspiral, std::mt19937_64 &azar) {
    char nombre[96], detalle[160];
    printf("\nÁrbol con espiral %s\n", nombreEspiral);

    //ejemplo chico: los 16 puntos de testUD.cpp con una z fija por punto
    {
        const int plano[16][2] = { {0,1},{7,3},{2,4},{7,6},{3,1},{2,7},{7,4},{8,9},
                                   {1,3},{4,2},{1,9},{3,7},{6,5},{8,1},{6,8},{9,9} };
        std::vector<Point3D> puntos;
        for (int i = 0; i < 16; ++i) puntos.push_back(Point3D(plano[i][0], plano[i][1], (plano[i][0] * plano[i][1]) % 5));
        const Arbol arbol(puntos);
        const Conjunto conjunto(puntos);
        bool ok = true;
        for (int64_t x = -1; x <= 10 && ok; ++x)
            for (int64_t y = -1; y <= 10 && ok; ++y)
                for (int64_t z = -1; z <= 5 && ok; ++z)
                    if (arbol.contains(Point3D(x, y, z)) != conjunto.contiene(Point3D(x, y, z))) ok = false;
        reportar("16 puntos de testUD.cpp en 3D, todo el cubo", ok, "12 x 12 x 7 consultas");
    }

    //cubos chicos (muchos empates y repetidos) hasta cubos de un millón por lado
    const int64_t lados[] = { 8, 1000, 1000000 };
    for (int64_t lado : lados) {
        const std::vector<Point3D> puntos = puntosAzar(azar, 20000, lado);
        const Arbol arbol(puntos);
        const Conjunto conjunto(puntos);
        uint64_t encontrados = 0;
        const bool ok = respondeIgual(arbol, conjunto, azar, lado, 100000, encontrados);
        snprintf(nombre, sizeof nombre, "20000 puntos en un cubo de lado %lld", (long long)lado);
        snprintf(detalle, sizeof detalle, "%zu distintos, %llu de 100000 consultas están",
                 conjunto.puntos.size(), (unsigned long long)encontrados);
        reportar(nombre, ok, detalle);
    }

    //rango: el conteo del árbol contra contar uno por uno, repetidos incluidos
    for (int64_t lado : lados) {
        const std::vector<Point3D> puntos = puntosAzar(azar, 20000, lado);
        const Arbol arbol(puntos);
        std::uniform_int_distribution<int64_t> coordenada(-1, lado);
        bool ok = true;
        uint64_t dentroTotal = 0;
        for (int caja = 0; caja < 2000 && ok; ++caja) {
            int64_t a[3], b[3];
            for (int eje = 0; eje < 3; ++eje) {
                a[eje] = coordenada(azar);
                b[eje] = coordenada(azar);
                if (a[eje] > b[eje]) std::swap(a[eje], b[eje]);
            }
            const Box3D box(Point3D(a[0], a[1], a[2]), Point3D(b[0], b[1], b[2]));
            uint64_t esperado = 0;
            for (const Point3D &p : puntos) esperado += box.contains(p);
            if (arbol.countInRange(box) != esperado) ok = false;
            dentroTotal += esperado;
        }
        snprintf(nombre, sizeof nombre, "rango en un cubo de lado %lld", (long long)lado);
        snprintf(detalle, sizeof detalle, "2000 cajas, %llu puntos dentro en total", (unsigned long long)dentroTotal);
        reportar(nombre, ok, detalle);
    }
    {
        std::vector<Point3D> puntos = puntosAzar(azar, 5000, 50);
        puntos.push_back(puntos[0]); //el mismo punto dos veces
        const Arbol arbol(puntos);
        const Arbol vacio(std::vector<Point3D>{});
        uint64_t repeticiones = 0;
        for (const Point3D &p : puntos) repeticiones += p == puntos[0];
        const bool ok = vacio.countInRange(Box3D::everything()) == 0
                     && arbol.countInRange(Box3D::everything()) == puntos.size()
                     && arbol.countInRange(Box3D(puntos[0], puntos[0])) == repeticiones
                     && arbol.countInRange(Box3D(Point3D(60, 60, 60), Point3D(90, 90, 90))) == 0;
        reportar("rango: vacío, todo el espacio, una celda y afuera", ok);
    }

    //bordes de cantidad
    {
        const Arbol vacio(std::vector<Point3D>{});
        const Arbol uno(std::vector<Point3D>{ Point3D(5, -3, 2) });
        const Arbol dos(std::vector<Point3D>{ Point3D(5, -3, 2), Point3D(-7, 4, 9) });
        const bool ok = !vacio.contains(Point3D(0, 0, 0))
                     && uno.contains(Point3D(5, -3, 2)) && !uno.contains(Point3D(5, -3, 3))
                     && dos.contains(Point3D(5, -3, 2)) && dos.contains(Point3D(-7, 4, 9))
                     && !dos.contains(Point3D(-7, 4, 8));
        reportar("árbol vacío, de un punto y de dos", ok);
    }

    //coordenadas negativas y lejanas: solo importa la distancia padre-hijo
    {
        std::vector<Point3D> puntos = puntosAzar(azar, 5000, 1000000);
        for (Point3D &p : puntos) p = Point3D(p.x - 4000000000000LL, p.y + 7000000000000LL, -p.z);
        const Arbol arbol(puntos);
        const Conjunto conjunto(puntos);
        bool ok = true;
        for (const Point3D &p : conjunto.puntos)
            if (!arbol.contains(p) || arbol.contains(Point3D(p.x, p.y, p.z + 1)) != conjunto.contiene(Point3D(p.x, p.y, p.z + 1)))
                ok = false;
        reportar("coordenadas negativas y lejos del origen", ok, "x - 4e12, y + 7e12, -z");
    }

    //guardar y cargar con tmpfile, que se borra solo
    {
        const std::vector<Point3D> puntos = puntosAzar(azar, 20000, 1000000);
        const Arbol arbol(puntos);
        Arbol cargado;
        FILE *archivo = tmpfile();
        bool ok = arbol.write(archivo);
        rewind(archivo);
        ok = ok && cargado.read(archivo) && cargado.size() == arbol.size();
        fclose(archivo);
        std::uniform_int_distribution<int64_t> coordenada(0, 999999);
        for (const Point3D &p : puntos) if (!cargado.contains(p)) ok = false;
        for (int i = 0; i < 100000 && ok; ++i) {
            const Point3D q(coordenada(azar), coordenada(azar), coordenada(azar));
            if (cargado.contains(q) != arbol.contains(q)) ok = false;
        }
        reportar("guardar y cargar responde igual", ok, "20000 puntos, 100000 consultas");
    }

    //en un cubo de un millón por lado hay códigos que no caben en 32 bits
    {
        const std::vector<Point3D> puntos = puntosAzar(azar, 20000, 1000000);
        const Arbol arbol(puntos);
        uint64_t sobre32Bits = 0, maximo = 0;
        for (uint64_t codigo : dac64Decompress(arbol.codes())) {
            if (codigo > 0xFFFFFFFFULL) ++sobre32Bits;
            maximo = std::max(maximo, codigo);
        }
        snprintf(detalle, sizeof detalle, "%llu códigos >= 2^32, máximo %llu",
                 (unsigned long long)sobre32Bits, (unsigned long long)maximo);
        reportar("códigos de 64 bits en un cubo de lado 10^6", sobre32Bits > 0, detalle);
    }

    //el cubo envolvente puede llegar justo a la capa máxima, pero no pasarla
    {
        const int64_t tope = (int64_t)Arbol::Spiral::maxDelta;
        const std::vector<Point3D> enElTope = { Point3D(0, 0, 0), Point3D(tope, 0, tope), Point3D(0, tope, 3) };
        bool ok = true;
        try {
            const Arbol arbol(enElTope);
            for (const Point3D &p : enElTope) if (!arbol.contains(p)) ok = false;
        } catch (const std::invalid_argument &) { ok = false; }

        bool rechazado = false;
        try { const Arbol arbol(std::vector<Point3D>{ Point3D(0, 0, 0), Point3D(tope + 1, 0, 0) }); }
        catch (const std::invalid_argument &) { rechazado = true; }
        snprintf(detalle, sizeof detalle, "capa máxima %llu", (unsigned long long)tope);
        reportar("cubo en la capa máxima sí, uno más no", ok && rechazado, detalle);
    }
}

int main() {
    printf("\nVerificación del kd-tree 3D compacto sobre el DAC de 64 bits\n");
    std::mt19937_64 azar(20261026);
    probarArbol<cKDtree3dUD>("unidireccional (UD)", azar);
    probarArbol<cKDtree3dBI>("bidireccional (BI)", azar);

    if (errores) printf("\nVERIFICACIÓN FALLIDA (%d prueba%s con error)\n\n",
                        errores, errores == 1 ? "" : "s");
    else         printf("\nTODAS LAS PRUEBAS PASARON\n\n");
    return errores ? 1 : 0;
}
