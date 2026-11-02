#ifndef CKDTREE3D_HPP
#define CKDTREE3D_HPP

#include "Spiral3dUD.hpp"
#include "Dac64.hpp"
#include "Dac64Extras.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <vector>

/**
 * kd-tree 3D compacto: port de ref/cKDtreeUD/cKDtreeUD.hpp a 3D y 64 bits.
 * Guarda la raíz y, por cada nodo, el código de la espiral respecto de su padre,
 * todos en un Dac64. La espiral es un parámetro: UD (semiespacio) o BI (cubo completo).
 */

/**
 * PASO A: construcción y recorrido compartidos.
 * Árbol implícito: el nodo de [start, end] es mid = (start + end) / 2, como en la referencia.
 * Eje de corte por profundidad: x, y, z, x, ...
 */

/** Espiral unidireccional: el lado del hijo lo da la rama. */
struct Spiral3dUDCodec {
    static const uint64_t maxDelta = SPIRAL3DUD_MAX_LAYER;
    static uint64_t encode(const Point3D &child, const Point3D &parent, Spiral3dAxis parentAxis) {
        return spiral3dUDEncode(child, parent, parentAxis);
    }
    static Point3D decode(uint64_t code, const Point3D &parent, Spiral3dAxis parentAxis, bool left) {
        return spiral3dUDDecode(code, parent, parentAxis, left);
    }
};

/** Espiral bidireccional: no usa el eje ni el lado. */
struct Spiral3dBICodec {
    static const uint64_t maxDelta = SPIRAL3D_MAX_LAYER;
    static uint64_t encode(const Point3D &child, const Point3D &parent, Spiral3dAxis) {
        return spiral3dEncode(child, parent);
    }
    static Point3D decode(uint64_t code, const Point3D &parent, Spiral3dAxis, bool) {
        return spiral3dDecode(code, parent);
    }
};

inline int64_t kdtree3dCoordinate(const Point3D &point, Spiral3dAxis axis) {
    return axis == SPIRAL3D_AXIS_X ? point.x : axis == SPIRAL3D_AXIS_Y ? point.y : point.z;
}

/**
 * Orden total del eje y luego los otros dos en orden cíclico (en 2D: el eje y luego el otro).
 * @return true si a va antes que b
 */
inline bool kdtree3dLess(const Point3D &a, const Point3D &b, Spiral3dAxis axis) {
    for (int step = 0; step < 3; ++step) {
        const Spiral3dAxis current = (Spiral3dAxis)((axis + step) % 3);
        const int64_t ca = kdtree3dCoordinate(a, current), cb = kdtree3dCoordinate(b, current);
        if (ca != cb) return ca < cb;
    }
    return false;
}

/**
 * Deja los puntos en el orden del árbol implícito (buildTree de la referencia).
 * @param points puntos, se reordenan
 * @param start primer índice del subárbol
 * @param end último índice del subárbol, incluido
 * @param depth profundidad del nodo de [start, end]
 */
inline void kdtree3dBuildOrder(std::vector<Point3D> &points, int64_t start, int64_t end, uint64_t depth) {
    if (start > end) return;
    const Spiral3dAxis axis = spiral3dUDAxisAt(depth);
    const int64_t mid = (start + end) / 2;
    std::nth_element(points.begin() + start, points.begin() + mid, points.begin() + end + 1,
                     [axis](const Point3D &a, const Point3D &b) { return kdtree3dLess(a, b, axis); });
    kdtree3dBuildOrder(points, start, mid - 1, depth + 1);
    kdtree3dBuildOrder(points, mid + 1, end, depth + 1);
}

/**
 * Busca un punto bajando desde el nodo de [start, end] (Link de la referencia).
 * Si la coordenada empata con el nodo, se buscan los dos lados, igual que en la referencia.
 * @param nodeAt nodeAt(index, parent, parentAxis, left) devuelve el punto del nodo
 * @param query punto buscado
 * @param parent padre del nodo de [start, end]; en la raíz, la raíz misma
 * @param start primer índice del subárbol
 * @param end último índice del subárbol, incluido
 * @param depth profundidad del nodo de [start, end]
 * @param left si el subárbol es hijo izquierdo
 * @return true si el punto está
 */
template <typename NodeAt>
inline bool kdtree3dSearch(const NodeAt &nodeAt, const Point3D &query, const Point3D &parent,
                           int64_t start, int64_t end, uint64_t depth, bool left) {
    if (start > end) return false;
    const int64_t mid = (start + end) / 2;
    const Spiral3dAxis parentAxis = spiral3dUDAxisAt(depth ? depth - 1 : 0);
    const Point3D current = depth ? nodeAt((uint64_t)mid, parent, parentAxis, left) : parent;
    if (query == current) return true;

    const Spiral3dAxis axis = spiral3dUDAxisAt(depth);
    const int64_t q = kdtree3dCoordinate(query, axis), c = kdtree3dCoordinate(current, axis);
    if (q < c) return kdtree3dSearch(nodeAt, query, current, start, mid - 1, depth + 1, true);
    if (q > c) return kdtree3dSearch(nodeAt, query, current, mid + 1, end, depth + 1, false);
    return kdtree3dSearch(nodeAt, query, current, start, mid - 1, depth + 1, true)
        || kdtree3dSearch(nodeAt, query, current, mid + 1, end, depth + 1, false);
}

/** Caja alineada a los ejes, con los bordes incluidos. */
struct Box3D {
    Point3D low, high;

    Box3D() {}
    Box3D(const Point3D &low, const Point3D &high) : low(low), high(high) {}

    /** Todo el espacio de 64 bits: la región de la raíz. */
    static Box3D everything() {
        return Box3D(Point3D(INT64_MIN, INT64_MIN, INT64_MIN), Point3D(INT64_MAX, INT64_MAX, INT64_MAX));
    }
    bool contains(const Point3D &p) const {
        return low.x <= p.x && p.x <= high.x && low.y <= p.y && p.y <= high.y && low.z <= p.z && p.z <= high.z;
    }
    bool intersects(const Box3D &o) const {
        return low.x <= o.high.x && o.low.x <= high.x && low.y <= o.high.y && o.low.y <= high.y
            && low.z <= o.high.z && o.low.z <= high.z;
    }
};

inline void kdtree3dSetCoordinate(Point3D &point, Spiral3dAxis axis, int64_t value) {
    if (axis == SPIRAL3D_AXIS_X) point.x = value;
    else if (axis == SPIRAL3D_AXIS_Y) point.y = value;
    else point.z = value;
}

/**
 * Cuenta los puntos dentro de una caja (RangeQuery de la referencia [2], algoritmo 4).
 * Cada subárbol conoce la región que cubre: si no toca la caja, se descarta sin decodificar.
 * Con empates el hijo izquierdo puede tener la coordenada del nodo, así que ambas regiones la incluyen.
 * @param nodeAt nodeAt(index, parent, parentAxis, left) devuelve el punto del nodo
 * @param box caja consultada
 * @param parent padre del nodo de [start, end]; en la raíz, la raíz misma
 * @param start primer índice del subárbol
 * @param end último índice del subárbol, incluido
 * @param depth profundidad del nodo de [start, end]
 * @param left si el subárbol es hijo izquierdo
 * @param region parte del espacio que cubre el subárbol
 * @return puntos del subárbol dentro de la caja, repetidos incluidos
 */
template <typename NodeAt>
inline uint64_t kdtree3dRangeCount(const NodeAt &nodeAt, const Box3D &box, const Point3D &parent,
                                   int64_t start, int64_t end, uint64_t depth, bool left, const Box3D &region) {
    if (start > end || !box.intersects(region)) return 0;
    const int64_t mid = (start + end) / 2;
    const Spiral3dAxis parentAxis = spiral3dUDAxisAt(depth ? depth - 1 : 0);
    const Point3D current = depth ? nodeAt((uint64_t)mid, parent, parentAxis, left) : parent;
    uint64_t count = box.contains(current) ? 1 : 0;

    const Spiral3dAxis axis = spiral3dUDAxisAt(depth);
    const int64_t c = kdtree3dCoordinate(current, axis);
    Box3D leftRegion = region, rightRegion = region;
    kdtree3dSetCoordinate(leftRegion.high, axis, c);
    kdtree3dSetCoordinate(rightRegion.low, axis, c);
    count += kdtree3dRangeCount(nodeAt, box, current, start, mid - 1, depth + 1, true, leftRegion);
    count += kdtree3dRangeCount(nodeAt, box, current, mid + 1, end, depth + 1, false, rightRegion);
    return count;
}

/**
 * Código de cada nodo respecto de su padre (generateCodes de la referencia).
 * El lado no se guarda: al decodificar lo da la rama. La raíz queda en 0.
 * @param points puntos en el orden de kdtree3dBuildOrder
 * @param codes un código por punto, mismo orden
 * @param parent padre del nodo de [start, end]; en la raíz, la raíz misma
 * @param start primer índice del subárbol
 * @param end último índice del subárbol, incluido
 * @param depth profundidad del nodo de [start, end]
 */
template <typename Codec>
inline void kdtree3dGenerateCodes(const std::vector<Point3D> &points, std::vector<uint64_t> &codes,
                                  const Point3D &parent, int64_t start, int64_t end, uint64_t depth) {
    if (start > end) return;
    const int64_t mid = (start + end) / 2;
    const Point3D &current = points[mid];
    if (depth) codes[mid] = Codec::encode(current, parent, spiral3dUDAxisAt(depth - 1));
    kdtree3dGenerateCodes<Codec>(points, codes, current, start, mid - 1, depth + 1);
    kdtree3dGenerateCodes<Codec>(points, codes, current, mid + 1, end, depth + 1);
}

/**
 * Revisa que el cubo envolvente quepa en la espiral: padre e hijo están dentro de él.
 * @throws std::invalid_argument si algún eje supera Codec::maxDelta
 */
template <typename Codec>
inline void kdtree3dCheckExtent(const std::vector<Point3D> &points) {
    Point3D low = points[0], high = points[0];
    for (const Point3D &p : points) {
        low  = Point3D(std::min(low.x, p.x), std::min(low.y, p.y), std::min(low.z, p.z));
        high = Point3D(std::max(high.x, p.x), std::max(high.y, p.y), std::max(high.z, p.z));
    }
    const uint64_t extent = std::max({ (uint64_t)(high.x - low.x), (uint64_t)(high.y - low.y),
                                       (uint64_t)(high.z - low.z) });
    if (extent > Codec::maxDelta)
        throw std::invalid_argument("CompactKdTree3d: el cubo envolvente supera la capa máxima de la espiral");
}

/**
 * Raíz y códigos de un conjunto de puntos, listos para guardar en un DAC.
 * @param points puntos, se reordenan al orden del árbol
 * @param root raíz del árbol
 * @return un código por punto en el orden del árbol
 */
template <typename Codec>
inline std::vector<uint64_t> kdtree3dEncodeTree(std::vector<Point3D> &points, Point3D &root) {
    kdtree3dCheckExtent<Codec>(points);
    kdtree3dBuildOrder(points, 0, (int64_t)points.size() - 1, 0);
    root = points[((int64_t)points.size() - 1) / 2];
    std::vector<uint64_t> codes(points.size(), 0);
    kdtree3dGenerateCodes<Codec>(points, codes, root, 0, (int64_t)points.size() - 1, 0);
    return codes;
}

/**
 * PASO B: árbol comprimido.
 * Cada nodo guarda Codec::encode(nodo, padre, eje del padre); la raíz se guarda aparte.
 * Cada paso del recorrido es O(1): un acceso al Dac64 y un decode de la espiral.
 */
template <typename Codec>
class CompactKdTree3d {
public:
    typedef Codec Spiral;

    /**
     * Construye el árbol. Las diferencias padre-hijo deben caber en la espiral.
     * @param points puntos; se aceptan repetidos
     * @throws std::invalid_argument si algún eje del cubo envolvente supera Codec::maxDelta
     */
    explicit CompactKdTree3d(std::vector<Point3D> points) {
        if (points.empty()) return;
        const std::vector<uint64_t> codes = kdtree3dEncodeTree<Codec>(points, root_);
        codes_ = dac64Build(codes.data(), codes.size());
    }

    CompactKdTree3d() {}

    /** true si el punto está en el árbol (checkLink de la referencia). */
    bool contains(const Point3D &query) const {
        if (codes_.count == 0) return false;
        const Dac64 &codes = codes_;
        const auto nodeAt = [&codes](uint64_t index, const Point3D &parent, Spiral3dAxis parentAxis, bool left) {
            return Codec::decode(dac64Access(codes, index), parent, parentAxis, left);
        };
        return kdtree3dSearch(nodeAt, query, root_, 0, (int64_t)codes_.count - 1, 0, true);
    }

    /** Puntos dentro de la caja, repetidos incluidos (RangeQuery de la referencia [2]). */
    uint64_t countInRange(const Box3D &box) const {
        if (codes_.count == 0) return 0;
        const Dac64 &codes = codes_;
        const auto nodeAt = [&codes](uint64_t index, const Point3D &parent, Spiral3dAxis parentAxis, bool left) {
            return Codec::decode(dac64Access(codes, index), parent, parentAxis, left);
        };
        return kdtree3dRangeCount(nodeAt, box, root_, 0, (int64_t)codes_.count - 1, 0, true, Box3D::everything());
    }

    uint64_t size() const { return codes_.count; }
    const Dac64 &codes() const { return codes_; }

    /** Bits en memoria: el Dac64 más la raíz. */
    uint64_t sizeInBits() const { return dac64SizeInBits(codes_) + 3 * 64; }

    /** Escribe el Dac64 y después la raíz (x, y, z), como cKDsave de la referencia. */
    bool write(FILE *file) const {
        const int64_t root[3] = { root_.x, root_.y, root_.z };
        return dac64Write(codes_, file) && fwrite(root, sizeof root[0], 3, file) == 3;
    }

    /** Lee un árbol escrito con write; queda vacío si falla. */
    bool read(FILE *file) {
        int64_t root[3] = { 0, 0, 0 };
        const bool ok = dac64Read(file, codes_) && fread(root, sizeof root[0], 3, file) == 3;
        if (!ok) { codes_ = Dac64(); root_ = Point3D(); return false; }
        root_ = Point3D(root[0], root[1], root[2]);
        return true;
    }

private:
    Point3D root_;
    Dac64 codes_;
};

/** kd-tree compacto con la espiral unidireccional (el de la tesis). */
typedef CompactKdTree3d<Spiral3dUDCodec> cKDtree3dUD;

/** kd-tree compacto con la espiral bidireccional, para comparar. */
typedef CompactKdTree3d<Spiral3dBICodec> cKDtree3dBI;

#endif
