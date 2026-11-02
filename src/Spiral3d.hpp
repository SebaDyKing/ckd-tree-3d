#ifndef SPIRAL3D_HPP
#define SPIRAL3D_HPP

#include <cstdint>
#include <cmath>
#include <cstdlib>

//Espiral 3D en 64 bits. La capa m ocupa [(2m-1)^3, (2m+1)^3-1], tope m = 1321122.

struct Point3D {
    int64_t x, y, z;
    Point3D() : x(0), y(0), z(0) {}
    Point3D(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}
    bool operator==(const Point3D &o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

//mayor impar s con s^3 < 2^64, y la capa que le corresponde
static const uint64_t SPIRAL3D_MAX_ODD_SIDE = 2642245ULL;
static const uint64_t SPIRAL3D_MAX_LAYER    = (SPIRAL3D_MAX_ODD_SIDE - 1) / 2;

//índice dentro de una cara, recorrida en filas alternadas
inline uint64_t spiral3dFaceIndex(uint64_t u, uint64_t v, uint64_t side) {
    return (u & 1ULL) ? u * side + v
                      : (u + 1) * side - v - 1;
}

//inverso: del índice a la posición dentro de la cara
inline void spiral3dFaceCoords(uint64_t index, uint64_t side,
                               uint64_t &u, uint64_t &v) {
    u = index / side;
    const uint64_t r = index % side;
    v = (u & 1ULL) ? r : side - 1 - r;
}

inline uint64_t spiral3dEncode(const Point3D &current, const Point3D &parent) {
    //calcular diferencias
    const int64_t dx = current.x - parent.x;
    const int64_t dy = current.y - parent.y;
    const int64_t dz = current.z - parent.z;
    if (dx == 0 && dy == 0 && dz == 0) return 0;

    //encontrar la capa
    int64_t m = dx < 0 ? -dx : dx;
    const int64_t ay = dy < 0 ? -dy : dy;
    const int64_t az = dz < 0 ? -dz : dz;
    if (ay > m) m = ay;
    if (az > m) m = az;

    const uint64_t a    = 2 * (uint64_t)m + 1;  //ancho de la cara
    const uint64_t s    = 2 * (uint64_t)m - 1;  //ancho de la capa anterior
    const uint64_t face = a * (a - 1);          //lo que se suma por cara lateral
    uint64_t base = s * s * s;                  //primer código de la capa

    //el orden importa: una arista se la queda la primera cara que calza
    if (dz ==  m) return base + spiral3dFaceIndex((uint64_t)(m - dx), (uint64_t)(dy + m), a);
    base += face;
    if (dx == -m) return base + spiral3dFaceIndex((uint64_t)(m - dz), (uint64_t)(dy + m), a);
    base += face;
    if (dz == -m) return base + spiral3dFaceIndex((uint64_t)(dx + m), (uint64_t)(dy + m), a);
    base += face;
    if (dx ==  m) return base + spiral3dFaceIndex((uint64_t)(dz + m), (uint64_t)(dy + m), a);
    base += face;

    //las dos caras restantes van sin las orillas, ya tomadas por las laterales
    const int64_t  im = m - 1;
    const uint64_t b  = 2 * (uint64_t)im + 1;   //ancho interior, es a-2
    if (dy == -m) return base + spiral3dFaceIndex((uint64_t)(im - dx), (uint64_t)(dz + im), b);
    base += b * b;

    //si nada calzó arriba entonces dy == m, porque m es la diferencia mayor
    return base + spiral3dFaceIndex((uint64_t)(dx + im), (uint64_t)(dz + im), b);
}

//capa desde el código: s = 2m-1 es el mayor impar cuyo cubo no pasa del código
//sin bucles, cbrt yerra a lo más 1 y los dos if lo corrigen
inline uint64_t spiral3dLayer(uint64_t code) {
    uint64_t t = (uint64_t)std::cbrt((double)code);
    if (t > 0 && t * t * t > code) --t;
    if (t < SPIRAL3D_MAX_ODD_SIDE && (t + 1) * (t + 1) * (t + 1) <= code) ++t;
    const uint64_t s = (t & 1ULL) ? t : t - 1;
    return (s + 1) / 2;
}

inline Point3D spiral3dDecode(uint64_t code, const Point3D &parent) {
    if (code == 0) return parent;

    const uint64_t m  = spiral3dLayer(code);
    const uint64_t a  = 2 * m + 1;
    const uint64_t s  = 2 * m - 1;
    const uint64_t d  = code - s * s * s;   //desplazamiento dentro de la capa
    const uint64_t a2 = a * a;
    const uint64_t b  = a - 2;              //ancho interior
    const int64_t  sm = (int64_t)m;
    const int64_t  im = sm - 1;

    uint64_t u, v;
    int64_t dx, dy, dz;

    //cada cara emite a^2, a^2-a, a^2-a, a^2-2a, b^2 y b^2 códigos, que suman
    //6a^2-12a+8, las celdas de la capa. Los cortes tienen holgura de una fila:
    //caen sobre aristas compartidas y ahí ambas caras dan el mismo punto
    if (d < a2) {                                   //cara 1, z+
        spiral3dFaceCoords(d, a, u, v);
        dx = sm - (int64_t)u;   dy = (int64_t)v - sm;   dz =  sm;
    } else if (d < 2 * a2 - a) {                    //cara 2, x-
        spiral3dFaceCoords(d - (a2 - a), a, u, v);
        dz = sm - (int64_t)u;   dy = (int64_t)v - sm;   dx = -sm;
    } else if (d < 3 * a2 - 2 * a) {                //cara 3, z-
        spiral3dFaceCoords(d - (2 * a2 - 2 * a), a, u, v);
        dx = (int64_t)u - sm;   dy = (int64_t)v - sm;   dz = -sm;
    } else if (d < 4 * a2 - 4 * a) {                //cara 4, x+
        spiral3dFaceCoords(d - (3 * a2 - 3 * a), a, u, v);
        dz = (int64_t)u - sm;   dy = (int64_t)v - sm;   dx =  sm;
    } else if (d < 4 * a2 - 4 * a + b * b) {        //cara 5, y-
        spiral3dFaceCoords(d - (4 * a2 - 4 * a), b, u, v);
        dx = im - (int64_t)u;   dz = (int64_t)v - im;   dy = -sm;
    } else {                                        //cara 6, y+
        spiral3dFaceCoords(d - (4 * a2 - 4 * a) - b * b, b, u, v);
        dx = (int64_t)u - im;   dz = (int64_t)v - im;   dy =  sm;
    }

    return Point3D(parent.x + dx, parent.y + dy, parent.z + dz);
}

#endif
