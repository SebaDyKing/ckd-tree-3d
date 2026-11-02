#ifndef SPIRAL3DUD_HPP
#define SPIRAL3DUD_HPP

#include "Spiral3d.hpp"

#include <cstdint>
#include <cmath>

/**
 * Espiral 3D unidireccional (semiespacio), 64 bits.
 * El hijo vive de un solo lado del plano de corte del padre.
 * El lado no se guarda en el código: lo da el subárbol.
 */

/**
 * PASO 1: dominio de búsqueda.
 * Marco canónico (x', y', z'): z' es el eje de corte plegado (z' >= 0).
 * Capa m = max(|x'|, |y'|, z').
 */

enum Spiral3dAxis { SPIRAL3D_AXIS_X = 0, SPIRAL3D_AXIS_Y = 1, SPIRAL3D_AXIS_Z = 2 };

/**
 * Eje de corte según la profundidad: x, y, z, x, ...
 * @param depth profundidad del nodo
 * @return eje de corte
 */
inline Spiral3dAxis spiral3dUDAxisAt(uint64_t depth) {
    return (Spiral3dAxis)(depth % 3);
}

/** Capa más grande cuyo count cabe en 64 bits. */
static const uint64_t SPIRAL3DUD_MAX_LAYER = 1664509ULL;

/**
 * Celdas hasta la capa m: (2m+1)^2 (m+1).
 * @param m capa
 * @return cantidad de códigos en [0, count)
 */
inline uint64_t spiral3dUDCount(uint64_t m) {
    const uint64_t a = 2 * m + 1;
    return a * a * (m + 1);
}

/**
 * Primer código de la capa m: (2m-1)^2 m.
 * @param m capa
 * @return base(m); 0 si m = 0
 */
inline uint64_t spiral3dUDBase(uint64_t m) {
    if (m == 0) return 0;
    const uint64_t s = 2 * m - 1;
    return s * s * m;
}

/**
 * Celdas de la capa m: 12m^2 + 4m + 1.
 * @param m capa
 * @return celdas propias de la capa
 */
inline uint64_t spiral3dUDShell(uint64_t m) {
    return 12 * m * m + 4 * m + 1;
}

/**
 * Capa desde el código, O(1): cbrt(code/4) más una corrección de +-1.
 * @param code código, menor que count(SPIRAL3DUD_MAX_LAYER)
 * @return capa m
 */
inline uint64_t spiral3dUDLayer(uint64_t code) {
    uint64_t t = (uint64_t)std::cbrt((double)code / 4.0);
    if (t > 0 && spiral3dUDBase(t) > code) --t;
    if (t < SPIRAL3DUD_MAX_LAYER && spiral3dUDBase(t + 1) <= code) ++t;
    return t;
}

/**
 * Capa de un punto canónico: max(|cx|, |cy|, cz).
 * @param cx x'
 * @param cy y'
 * @param cz z', >= 0
 * @return capa m
 */
inline uint64_t spiral3dUDLayerOf(int64_t cx, int64_t cy, int64_t cz) {
    int64_t m = cx < 0 ? -cx : cx;
    const int64_t ay = cy < 0 ? -cy : cy;
    if (ay > m) m = ay;
    if (cz > m) m = cz;
    return (uint64_t)m;
}

/**
 * Desplazamiento real -> canónico: rota el eje de corte a z' y lo pliega con |.|.
 * @param dx, dy, dz desplazamiento hijo - padre
 * @param axis eje de corte del padre
 * @param cx, cy, cz salida: punto canónico
 */
inline void spiral3dUDToCanonical(int64_t dx, int64_t dy, int64_t dz, Spiral3dAxis axis,
                                  int64_t &cx, int64_t &cy, int64_t &cz) {
    switch (axis) {
    case SPIRAL3D_AXIS_X: cx = dy; cy = dz; cz = dx < 0 ? -dx : dx; break;
    case SPIRAL3D_AXIS_Y: cx = dz; cy = dx; cz = dy < 0 ? -dy : dy; break;
    default:              cx = dx; cy = dy; cz = dz < 0 ? -dz : dz; break;
    }
}

/**
 * Canónico -> real.
 * @param cx, cy, cz punto canónico
 * @param axis eje de corte del padre
 * @param left true si el hijo está del lado izquierdo (signo negativo)
 * @param dx, dy, dz salida: desplazamiento hijo - padre
 */
inline void spiral3dUDFromCanonical(int64_t cx, int64_t cy, int64_t cz, Spiral3dAxis axis, bool left,
                                    int64_t &dx, int64_t &dy, int64_t &dz) {
    const int64_t cut = left ? -cz : cz;
    switch (axis) {
    case SPIRAL3D_AXIS_X: dx = cut; dy = cx; dz = cy; break;
    case SPIRAL3D_AXIS_Y: dy = cut; dz = cx; dx = cy; break;
    default:              dz = cut; dx = cx; dy = cy; break;
    }
}

/**
 * PASO 2: orden dentro de la capa.
 * Se barre de x' = -m a x' = +m, sin saltos:
 * cara x' = -m (zigzag) -> rebanadas "U" -> cara x' = +m (zigzag).
 */

/**
 * Posición dentro de una U: sube por y' = +m, cruza por z' = m, baja por y' = -m.
 * @param cy y'
 * @param cz z'
 * @param m capa
 * @return q en [0, 4m]
 */
inline uint64_t spiral3dUDSlicePos(int64_t cy, int64_t cz, uint64_t m) {
    const int64_t sm = (int64_t)m;
    if (cy == sm) return (uint64_t)cz;
    if (cz == sm && cy > -sm) return (uint64_t)(2 * sm - cy);
    return 4 * m - (uint64_t)cz;
}

/**
 * Inverso de spiral3dUDSlicePos.
 * @param q posición en [0, 4m]
 * @param m capa
 * @param cy, cz salida: (y', z')
 */
inline void spiral3dUDSliceCell(uint64_t q, uint64_t m, int64_t &cy, int64_t &cz) {
    const int64_t sm = (int64_t)m;
    if (q <= m)          { cy = sm;                  cz = (int64_t)q; }
    else if (q < 3 * m)  { cy = 2 * sm - (int64_t)q; cz = sm; }
    else                 { cy = -sm;                 cz = 4 * sm - (int64_t)q; }
}

/**
 * Celda canónica -> posición dentro de la capa.
 * @param cx, cy, cz celda canónica de la capa m
 * @param m capa, >= 1
 * @return d en [0, shell(m))
 */
inline uint64_t spiral3dUDShellIndex(int64_t cx, int64_t cy, int64_t cz, uint64_t m) {
    const int64_t  sm = (int64_t)m;
    const uint64_t h  = m + 1;               //alto de una columna
    const uint64_t F  = (2 * m + 1) * h;     //celdas de una cara
    const uint64_t U  = 4 * m + 1;           //celdas de una U
    const uint64_t j  = (uint64_t)(cy + sm); //columna
    //cara x' = -m
    if (cx == -sm) return j * h + ((j & 1ULL) ? (uint64_t)cz : m - (uint64_t)cz);
    //cara x' = +m
    if (cx ==  sm) return F + (2 * m - 1) * U + j * h + ((j & 1ULL) ? m - (uint64_t)cz : (uint64_t)cz);
    //rebanada U número t; las impares van al revés
    const uint64_t t = (uint64_t)(cx + sm - 1);
    uint64_t q = spiral3dUDSlicePos(cy, cz, m);
    if (t & 1ULL) q = 4 * m - q;
    return F + t * U + q;
}

/**
 * Inverso: posición dentro de la capa -> celda canónica.
 * @param d posición en [0, shell(m))
 * @param m capa, >= 1
 * @param cx, cy, cz salida: celda canónica
 */
inline void spiral3dUDShellCell(uint64_t d, uint64_t m, int64_t &cx, int64_t &cy, int64_t &cz) {
    const int64_t  sm = (int64_t)m;
    const uint64_t h  = m + 1;
    const uint64_t F  = (2 * m + 1) * h;
    const uint64_t U  = 4 * m + 1;
    const uint64_t inner = (2 * m - 1) * U;
    //cara x' = -m
    if (d < F) {
        const uint64_t j = d / h, i = d % h;
        cx = -sm; cy = (int64_t)j - sm; cz = (int64_t)((j & 1ULL) ? i : m - i);
        return;
    }
    //cara x' = +m
    if (d >= F + inner) {
        const uint64_t r = d - F - inner;
        const uint64_t j = r / h, i = r % h;
        cx = sm; cy = (int64_t)j - sm; cz = (int64_t)((j & 1ULL) ? m - i : i);
        return;
    }
    //rebanada U
    const uint64_t r = d - F;
    const uint64_t t = r / U;
    uint64_t q = r % U;
    if (t & 1ULL) q = 4 * m - q;
    cx = (int64_t)t - sm + 1;
    spiral3dUDSliceCell(q, m, cy, cz);
}

/**
 * PASO 3: codificar y decodificar.
 * código = base(m) + d; el padre es el código 0.
 */

/**
 * Codifica un punto respecto de su padre.
 * @param current punto hijo; |current - parent| <= SPIRAL3DUD_MAX_LAYER por eje
 * @param parent punto padre
 * @param axis eje de corte del padre (spiral3dUDAxisAt(depth))
 * @return código; 0 si current == parent
 */
inline uint64_t spiral3dUDEncode(const Point3D &current, const Point3D &parent, Spiral3dAxis axis) {
    int64_t cx, cy, cz;
    spiral3dUDToCanonical(current.x - parent.x, current.y - parent.y, current.z - parent.z, axis, cx, cy, cz);
    const uint64_t m = spiral3dUDLayerOf(cx, cy, cz);
    if (m == 0) return 0;
    return spiral3dUDBase(m) + spiral3dUDShellIndex(cx, cy, cz, m);
}

/**
 * Decodifica un código respecto de su padre.
 * @param code código, menor que count(SPIRAL3DUD_MAX_LAYER)
 * @param parent punto padre
 * @param axis eje de corte del padre (spiral3dUDAxisAt(depth))
 * @param left lado del hijo, lo da el subárbol
 * @return punto hijo
 */
inline Point3D spiral3dUDDecode(uint64_t code, const Point3D &parent, Spiral3dAxis axis, bool left) {
    if (code == 0) return parent;
    const uint64_t m = spiral3dUDLayer(code);
    int64_t cx, cy, cz, dx, dy, dz;
    spiral3dUDShellCell(code - spiral3dUDBase(m), m, cx, cy, cz);
    spiral3dUDFromCanonical(cx, cy, cz, axis, left, dx, dy, dz);
    return Point3D(parent.x + dx, parent.y + dy, parent.z + dz);
}

#endif
