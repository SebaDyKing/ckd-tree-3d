/* Pruebas del dominio y del recorrido de la espiral 3D unidireccional
   (semiespacio).

     ./spiral3dud_verify [capas]     por defecto 40

   Se verifica el dominio (fórmulas de conteo, capa desde el código y marco
   canónico que pliega el eje de corte), el recorrido dentro de una capa
   (desplazamiento de cada celda y su inverso) y la codificación completa
   spiral3dUDEncode / spiral3dUDDecode respecto de un padre, un eje y un lado. */

#include "Spiral3dUD.hpp"

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

//máximo de los tres valores absolutos, la capa en el cubo completo
static int64_t capaCubo(int64_t x, int64_t y, int64_t z) {
    int64_t m = llabs(x);
    if (llabs(y) > m) m = llabs(y);
    if (llabs(z) > m) m = llabs(z);
    return m;
}

struct Triple {
    int64_t a, b, c;
    bool operator<(const Triple &o) const {
        if (a != o.a) return a < o.a;
        if (b != o.b) return b < o.b;
        return c < o.c;
    }
    bool operator==(const Triple &o) const { return a == o.a && b == o.b && c == o.c; }
};

//componente del desplazamiento sobre el eje de corte
static int64_t enEje(int64_t dx, int64_t dy, int64_t dz, Spiral3dAxis eje) {
    return eje == SPIRAL3D_AXIS_X ? dx : eje == SPIRAL3D_AXIS_Y ? dy : dz;
}

int main(int argc, char **argv) {
    const int64_t N = (argc > 1) ? atoll(argv[1]) : 40;

    printf("\nVerificación del dominio y recorrido de la espiral 3D unidireccional, capas 1..%lld\n",
           (long long)N);
    printf("%llu celdas hasta la capa %lld\n\n",
           (unsigned long long)spiral3dUDCount((uint64_t)N), (long long)N);

    //capacidad: la capa máxima debe caber y la siguiente no
    {
        typedef unsigned __int128 u128;
        const u128 limite = (u128)1 << 64;
        const u128 m      = SPIRAL3DUD_MAX_LAYER;
        const u128 a      = 2 * m + 1;
        const u128 a2     = 2 * (m + 1) + 1;
        const bool cabe     = a * a * (m + 1) <= limite;
        const bool ajustado = a2 * a2 * (m + 2) > limite;
        const bool ordenado = spiral3dUDBase(SPIRAL3DUD_MAX_LAYER)
                            < spiral3dUDCount(SPIRAL3DUD_MAX_LAYER);

        char buf[128];
        snprintf(buf, sizeof buf, "capa máxima m = %llu, count = %llu",
                 (unsigned long long)SPIRAL3DUD_MAX_LAYER,
                 (unsigned long long)spiral3dUDCount(SPIRAL3DUD_MAX_LAYER));
        reportar("la capa más grande cabe en un código de 64 bits",
                 cabe && ajustado && ordenado, buf);
    }

    /* Fórmulas contra conteo directo: se recorre la caja de la capa m y se
       cuentan las celdas por capa, que deben coincidir con shell, count y base. */
    {
        bool shellOk = true, countOk = true, baseOk = true, restaOk = true;
        int64_t capaMala = 0;

        for (int64_t m = 0; m <= N && shellOk && countOk && baseOk && restaOk; ++m) {
            uint64_t propias = 0, hasta = 0;
            for (int64_t x = -m; x <= m; ++x)
            for (int64_t y = -m; y <= m; ++y)
            for (int64_t z = 0; z <= m; ++z) {
                const uint64_t capa = spiral3dUDLayerOf(x, y, z);
                if (capa == (uint64_t)m) ++propias;
                if (capa <= (uint64_t)m) ++hasta;
            }
            const uint64_t um = (uint64_t)m;
            if (propias != spiral3dUDShell(um)) { shellOk = false; capaMala = m; }
            if (hasta   != spiral3dUDCount(um)) { countOk = false; capaMala = m; }
            const uint64_t previo = m == 0 ? 0 : spiral3dUDCount(um - 1);
            if (spiral3dUDBase(um) != previo) { baseOk = false; capaMala = m; }
            if (spiral3dUDCount(um) - spiral3dUDBase(um) != spiral3dUDShell(um)) {
                restaOk = false; capaMala = m;
            }
        }

        char buf[96];
        snprintf(buf, sizeof buf, "capa m=%lld", (long long)capaMala);
        reportar("shell(m) cuenta las celdas propias de la capa", shellOk, shellOk ? "" : buf);
        reportar("count(m) cuenta las celdas hasta la capa m", countOk, countOk ? "" : buf);
        reportar("base(m) == count(m-1), y 0 en la capa 0", baseOk, baseOk ? "" : buf);
        reportar("count(m) - base(m) == shell(m)", restaOk, restaOk ? "" : buf);
    }

    /* spiral3dUDLayer no tiene bucles: corrige el cbrt con un solo if en cada
       sentido. Acá se recorren TODOS los bordes de capa, que es donde el
       redondeo de coma flotante se da vuelta, para confirmar que un paso basta
       en todo el rango de 64 bits. */
    {
        bool ok = spiral3dUDLayer(0) == 0;
        uint64_t mMala = 0;
        uint64_t bordes = 1;

        for (uint64_t m = 1; m <= SPIRAL3DUD_MAX_LAYER; ++m) {
            const uint64_t lo = spiral3dUDBase(m);
            const uint64_t hi = spiral3dUDCount(m) - 1;
            bool bien = true;
            if (spiral3dUDLayer(lo)     != m)     bien = false;   //primer código
            if (spiral3dUDLayer(lo - 1) != m - 1) bien = false;   //último de la anterior
            if (spiral3dUDLayer(hi)     != m)     bien = false;   //último código
            if (!bien && ok) { ok = false; mMala = m; }
            bordes += 3;
        }

        char buf[128];
        if (ok) snprintf(buf, sizeof buf, "%llu bordes, sin bucles",
                         (unsigned long long)bordes);
        else    snprintf(buf, sizeof buf, "falla en m=%llu", (unsigned long long)mMala);
        reportar("spiral3dUDLayer acierta en todos los bordes de capa", ok, buf);
    }

    /* Códigos al azar en todo el rango: la capa devuelta tiene que contener
       al código, base(l) <= c < count(l). */
    {
        std::mt19937_64 gen(20260905);
        std::uniform_int_distribution<uint64_t> codigo(0, spiral3dUDCount(SPIRAL3DUD_MAX_LAYER) - 1);

        bool ok = true;
        uint64_t codigoMalo = 0;
        for (int i = 0; i < 1000000 && ok; ++i) {
            const uint64_t c = codigo(gen);
            const uint64_t l = spiral3dUDLayer(c);
            if (l > SPIRAL3DUD_MAX_LAYER || spiral3dUDBase(l) > c || c >= spiral3dUDCount(l)) {
                ok = false; codigoMalo = c;
            }
        }

        char buf[128];
        if (ok) snprintf(buf, sizeof buf, "1000000 códigos al azar");
        else    snprintf(buf, sizeof buf, "código %llu", (unsigned long long)codigoMalo);
        reportar("base(l) <= código < count(l) con l = spiral3dUDLayer", ok, buf);
    }

    /* Marco canónico, ida y vuelta: por cada eje de corte y cada lado, todo
       punto del cubo que viva de ese lado tiene que plegarse con z' >= 0,
       conservar su capa, volver intacto y no chocar con otro punto del lado. */
    {
        bool plegado = true, capaOk = true, vuelta = true, unico = true;
        int ejeMalo = -1; bool ladoMalo = false;
        int64_t bx = 0, by = 0, bz = 0;

        for (int e = 0; e < 3; ++e)
        for (int l = 0; l < 2; ++l) {
            const Spiral3dAxis eje = (Spiral3dAxis)e;
            const bool izq = l == 0;
            std::vector<Triple> vistos;
            vistos.reserve((size_t)((2 * N + 1) * (2 * N + 1) * (N + 1)));

            for (int64_t x = -N; x <= N; ++x)
            for (int64_t y = -N; y <= N; ++y)
            for (int64_t z = -N; z <= N; ++z) {
                const int64_t d = enEje(x, y, z, eje);
                if (izq ? d > 0 : d < 0) continue;      //no está de este lado

                int64_t cx, cy, cz;
                spiral3dUDToCanonical(x, y, z, eje, cx, cy, cz);
                const Triple t = { cx, cy, cz };
                vistos.push_back(t);

                bool bien = true;
                if (cz < 0) { plegado = false; bien = false; }
                if (spiral3dUDLayerOf(cx, cy, cz) != (uint64_t)capaCubo(x, y, z)) {
                    capaOk = false; bien = false;
                }
                int64_t rx, ry, rz;
                spiral3dUDFromCanonical(cx, cy, cz, eje, izq, rx, ry, rz);
                if (rx != x || ry != y || rz != z) { vuelta = false; bien = false; }
                if (!bien && ejeMalo < 0) { ejeMalo = e; ladoMalo = izq; bx = x; by = y; bz = z; }
            }

            std::sort(vistos.begin(), vistos.end());
            if (std::adjacent_find(vistos.begin(), vistos.end()) != vistos.end()) {
                unico = false;
                if (ejeMalo < 0) { ejeMalo = e; ladoMalo = izq; }
            }
        }

        char buf[160];
        snprintf(buf, sizeof buf, "eje %c, lado %s, punto (%lld,%lld,%lld)",
                 ejeMalo < 0 ? '-' : "xyz"[ejeMalo], ladoMalo ? "izquierdo" : "derecho",
                 (long long)bx, (long long)by, (long long)bz);
        reportar("el marco canónico pliega el eje de corte a z' >= 0", plegado, plegado ? "" : buf);
        reportar("la capa canónica es max(|dx|,|dy|,|dz|)", capaOk, capaOk ? "" : buf);
        reportar("FromCanonical(ToCanonical(p)) == p en cada lado", vuelta, vuelta ? "" : buf);
        reportar("ningún par de puntos del lado comparte marco canónico", unico, unico ? "" : buf);
    }

    /* Plano de empate: los puntos sobre el plano del padre valen para ambos
       lados, así que FromCanonical tiene que devolverlos con left y sin left. */
    {
        bool ok = true;
        int ejeMalo = -1;
        int64_t bx = 0, by = 0, bz = 0;

        for (int e = 0; e < 3 && ok; ++e) {
            const Spiral3dAxis eje = (Spiral3dAxis)e;
            for (int64_t u = -N; u <= N && ok; ++u)
            for (int64_t v = -N; v <= N && ok; ++v) {
                //u, v recorren los dos ejes libres; el eje de corte va en 0
                int64_t p[3] = { u, v, 0 };
                if (e == 0)      { p[0] = 0; p[1] = u; p[2] = v; }
                else if (e == 1) { p[1] = 0; p[0] = u; p[2] = v; }

                int64_t cx, cy, cz;
                spiral3dUDToCanonical(p[0], p[1], p[2], eje, cx, cy, cz);
                int64_t lx, ly, lz, rx, ry, rz;
                spiral3dUDFromCanonical(cx, cy, cz, eje, true,  lx, ly, lz);
                spiral3dUDFromCanonical(cx, cy, cz, eje, false, rx, ry, rz);
                const bool izq = lx == p[0] && ly == p[1] && lz == p[2];
                const bool der = rx == p[0] && ry == p[1] && rz == p[2];
                if (cz != 0 || !izq || !der) { ok = false; ejeMalo = e; bx = p[0]; by = p[1]; bz = p[2]; }
            }
        }

        char buf[160];
        if (ok) snprintf(buf, sizeof buf, "%lld puntos por eje",
                         (long long)((2 * N + 1) * (2 * N + 1)));
        else    snprintf(buf, sizeof buf, "eje %c, punto (%lld,%lld,%lld)", "xyz"[ejeMalo],
                         (long long)bx, (long long)by, (long long)bz);
        reportar("en el plano de empate ambos lados devuelven el punto", ok, buf);
    }

    /* Recorrido de la capa, biyección: cada celda propia de la capa m recibe un
       desplazamiento distinto en [0, shell(m)) y el inverso la devuelve. */
    {
        bool rango = true, unico = true, vuelta = true;
        int64_t mMala = 0, bx = 0, by = 0, bz = 0;
        uint64_t celdas = 0;

        for (int64_t m = 1; m <= N && rango && unico && vuelta; ++m) {
            const uint64_t um = (uint64_t)m;
            const uint64_t shell = spiral3dUDShell(um);
            std::vector<unsigned char> visto(shell, 0);
            uint64_t propias = 0;

            for (int64_t x = -m; x <= m; ++x)
            for (int64_t y = -m; y <= m; ++y)
            for (int64_t z = 0; z <= m; ++z) {
                if (spiral3dUDLayerOf(x, y, z) != um) continue;
                ++propias;
                const uint64_t d = spiral3dUDShellIndex(x, y, z, um);
                bool bien = true;
                if (d >= shell) { rango = false; bien = false; }
                else if (visto[d]) { unico = false; bien = false; }
                else visto[d] = 1;
                if (bien) {
                    int64_t rx, ry, rz;
                    spiral3dUDShellCell(d, um, rx, ry, rz);
                    if (rx != x || ry != y || rz != z) { vuelta = false; bien = false; }
                }
                if (!bien && mMala == 0) { mMala = m; bx = x; by = y; bz = z; }
            }
            //propias == shell ya se probó arriba; con rango y unico, el mapa
            //de bits queda lleno y la asignación es sobreyectiva
            if (propias != shell && mMala == 0) { rango = false; mMala = m; }
            celdas += propias;
        }

        char buf[160];
        if (mMala == 0) snprintf(buf, sizeof buf, "%llu celdas en %lld capas",
                                 (unsigned long long)celdas, (long long)N);
        else snprintf(buf, sizeof buf, "capa m=%lld, celda (%lld,%lld,%lld)",
                      (long long)mMala, (long long)bx, (long long)by, (long long)bz);
        reportar("cada celda de la capa recibe un desplazamiento distinto en [0, shell)",
                 rango && unico, buf);
        reportar("ShellCell(ShellIndex(c)) == c en toda la capa", vuelta, buf);
    }

    /* Continuidad: dos desplazamientos consecutivos son celdas vecinas cara
       con cara, es decir a distancia Manhattan 1. */
    {
        uint64_t pares = 0, saltos = 0;
        int64_t mMala = 0; uint64_t dMalo = 0;

        for (int64_t m = 1; m <= N; ++m) {
            const uint64_t um = (uint64_t)m;
            const uint64_t shell = spiral3dUDShell(um);
            int64_t px, py, pz;
            spiral3dUDShellCell(0, um, px, py, pz);
            for (uint64_t d = 1; d < shell; ++d) {
                int64_t x, y, z;
                spiral3dUDShellCell(d, um, x, y, z);
                const int64_t dist = llabs(x - px) + llabs(y - py) + llabs(z - pz);
                if (dist != 1) { if (saltos == 0) { mMala = m; dMalo = d - 1; } ++saltos; }
                px = x; py = y; pz = z;
                ++pares;
            }
        }

        char buf[160];
        if (saltos == 0) snprintf(buf, sizeof buf, "%llu pares consecutivos, 0 saltos",
                                  (unsigned long long)pares);
        else snprintf(buf, sizeof buf, "%llu saltos, el primero en m=%lld entre d=%llu y d+1",
                      (unsigned long long)saltos, (long long)mMala, (unsigned long long)dMalo);
        reportar("el recorrido es continuo", saltos == 0, buf);
    }

    /* Empalmes del recorrido: la capa empieza en (-m,-m,m), la cara x' = -m
       termina en (-m,+m,0) y la primera rebanada U empieza justo al lado en
       (-m+1,+m,0); la última rebanada termina en (m-1,-m,0) y la cara
       x' = +m empieza al lado en (+m,-m,0); la capa termina en (m,m,m). */
    {
        bool ok = true;
        int64_t mMala = 0;

        for (int64_t m = 1; m <= N && ok; ++m) {
            const uint64_t um = (uint64_t)m;
            const uint64_t F = (2 * um + 1) * (um + 1);
            const uint64_t inner = (2 * um - 1) * (4 * um + 1);
            int64_t x, y, z;
            spiral3dUDShellCell(0, um, x, y, z);
            if (x != -m || y != -m || z != m) ok = false;
            spiral3dUDShellCell(F - 1, um, x, y, z);
            if (x != -m || y != m || z != 0) ok = false;
            spiral3dUDShellCell(F, um, x, y, z);
            if (x != -m + 1 || y != m || z != 0) ok = false;
            spiral3dUDShellCell(F + inner - 1, um, x, y, z);
            if (x != m - 1 || y != -m || z != 0) ok = false;
            spiral3dUDShellCell(F + inner, um, x, y, z);
            if (x != m || y != -m || z != 0) ok = false;
            spiral3dUDShellCell(spiral3dUDShell(um) - 1, um, x, y, z);
            if (x != m || y != m || z != m) ok = false;
            if (!ok) mMala = m;
        }

        char buf[96];
        if (ok) snprintf(buf, sizeof buf, "capas 1..%lld", (long long)N);
        else    snprintf(buf, sizeof buf, "capa m=%lld", (long long)mMala);
        reportar("las caras y las rebanadas U empalman en los extremos", ok, buf);
    }

    /* Orden esperado de la capa 1, celda por celda (códigos 1..17). */
    {
        static const int64_t orden[17][3] = {
            { -1, -1, 1 }, { -1, -1, 0 }, { -1, 0, 0 }, { -1, 0, 1 }, { -1, 1, 1 }, { -1, 1, 0 },
            { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 }, { 0, -1, 1 }, { 0, -1, 0 },
            { 1, -1, 0 }, { 1, -1, 1 }, { 1, 0, 1 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 },
        };
        bool ok = true;
        for (uint64_t d = 0; d < 17; ++d) {
            int64_t x, y, z;
            spiral3dUDShellCell(d, 1, x, y, z);
            if (x != orden[d][0] || y != orden[d][1] || z != orden[d][2]) ok = false;
        }
        reportar("el orden de la capa 1 coincide con el documentado", ok, ok ? "17 celdas" : "");
    }

    /* Numeración esperada de las capas 1 y 2:
       filas x' = +m..-m, columnas y' = -m..+m, -1 = celda de otra capa. */
    {
        static const int capa1[2][3][3] = {
            { { 13, 14, 17 }, { 10, 9, 8 }, { 1, 4, 5 } },       //z' = 1
            { { 12, 15, 16 }, { 11, -1, 7 }, { 2, 3, 6 } },      //z' = 0
        };
        static const int capa2[3][5][5] = {
            { { 62, 63, 68, 69, 74 }, { 57, 56, 55, 54, 53 }, { 44, 45, 46, 47, 48 },
              { 39, 38, 37, 36, 35 }, { 18, 23, 24, 29, 30 } },  //z' = 2
            { { 61, 64, 67, 70, 73 }, { 58, -1, -1, -1, 52 }, { 43, -1, -1, -1, 49 },
              { 40, -1, -1, -1, 34 }, { 19, 22, 25, 28, 31 } },  //z' = 1
            { { 60, 65, 66, 71, 72 }, { 59, -1, -1, -1, 51 }, { 42, -1, -1, -1, 50 },
              { 41, -1, -1, -1, 33 }, { 20, 21, 26, 27, 32 } },  //z' = 0
        };

        bool ok = true;
        int64_t mMala = 0, bx = 0, by = 0, bz = 0;
        uint64_t celdas = 0;
        for (int64_t m = 1; m <= 2 && ok; ++m) {
            const uint64_t um = (uint64_t)m;
            for (int64_t z = m; z >= 0 && ok; --z)
            for (int64_t x = m; x >= -m && ok; --x)
            for (int64_t y = -m; y <= m && ok; ++y) {
                const int esperado = m == 1 ? capa1[m - z][m - x][y + m]
                                            : capa2[m - z][m - x][y + m];
                const bool propia = spiral3dUDLayerOf(x, y, z) == um;
                int64_t codigo = -1;
                if (propia) codigo = (int64_t)(spiral3dUDBase(um) + spiral3dUDShellIndex(x, y, z, um));
                if (codigo != esperado) { ok = false; mMala = m; bx = x; by = y; bz = z; }
                ++celdas;
            }
        }

        char buf[160];
        if (ok) snprintf(buf, sizeof buf, "%llu celdas de las cajas, códigos 1..74",
                         (unsigned long long)celdas);
        else    snprintf(buf, sizeof buf, "capa m=%lld, celda (%lld,%lld,%lld)", (long long)mMala,
                         (long long)bx, (long long)by, (long long)bz);
        reportar("la numeración de las capas 1 y 2 coincide con la documentada", ok, buf);
    }

    /* Ida y vuelta en capas grandes: celdas al azar de la cáscara, con una
       coordenada fijada al borde para que caigan en la capa. */
    {
        static const uint64_t capas[3] = { 813, 100000, SPIRAL3DUD_MAX_LAYER };
        std::mt19937_64 gen(20260905);

        bool ok = true;
        uint64_t mMala = 0, dMax = 0;
        int64_t bx = 0, by = 0, bz = 0;
        for (int c = 0; c < 3 && ok; ++c) {
            const uint64_t m  = capas[c];
            const int64_t  sm = (int64_t)m;
            const uint64_t shell = spiral3dUDShell(m);
            std::uniform_int_distribution<int64_t> libre(-sm, sm);
            std::uniform_int_distribution<int64_t> alto(0, sm);
            std::uniform_int_distribution<int>     cara(0, 4);

            for (int i = 0; i < 20000 && ok; ++i) {
                int64_t x = libre(gen), y = libre(gen), z = alto(gen);
                switch (cara(gen)) {
                case 0:  z = sm;  break;
                case 1:  x = sm;  break;
                case 2:  x = -sm; break;
                case 3:  y = sm;  break;
                default: y = -sm; break;
                }
                const uint64_t d = spiral3dUDShellIndex(x, y, z, m);
                int64_t rx, ry, rz;
                spiral3dUDShellCell(d < shell ? d : 0, m, rx, ry, rz);
                if (spiral3dUDLayerOf(x, y, z) != m || d >= shell
                    || rx != x || ry != y || rz != z
                    || spiral3dUDBase(m) + d >= spiral3dUDCount(m)) {
                    ok = false; mMala = m; bx = x; by = y; bz = z;
                }
                if (d > dMax) dMax = d;
            }
        }

        char buf[160];
        if (ok) snprintf(buf, sizeof buf, "3 capas, 20000 celdas por capa, desplazamiento máximo %llu",
                         (unsigned long long)dMax);
        else    snprintf(buf, sizeof buf, "capa m=%llu, celda (%lld,%lld,%lld)", (unsigned long long)mMala,
                         (long long)bx, (long long)by, (long long)bz);
        reportar("ida y vuelta en capas grandes", ok, buf);
    }

    /* Codificación completa, por cada eje y cada lado, con el padre en el
       origen: se recorre todo el cubo hasta la capa N y se queda con los
       puntos del lado (eje de corte <= 0 a la izquierda, >= 0 a la derecha). */
    {
        bool vuelta = true, biyeccion = true, cero = true, empate = true, compone = true;
        int ejeMalo = -1; bool ladoMalo = false;
        int64_t bx = 0, by = 0, bz = 0;
        uint64_t puntos = 0;
        const Point3D origen;
        const uint64_t total = spiral3dUDCount((uint64_t)N);

        for (int e = 0; e < 3; ++e)
        for (int l = 0; l < 2; ++l) {
            const Spiral3dAxis eje = (Spiral3dAxis)e;
            const bool izq = l == 0;
            std::vector<unsigned char> visto(total, 0);
            uint64_t cubiertos = 0;

            for (int64_t x = -N; x <= N; ++x)
            for (int64_t y = -N; y <= N; ++y)
            for (int64_t z = -N; z <= N; ++z) {
                const int64_t d = enEje(x, y, z, eje);
                if (izq ? d > 0 : d < 0) continue;      //no está de este lado

                const Point3D p(x, y, z);
                const uint64_t code = spiral3dUDEncode(p, origen, eje);
                ++puntos;
                bool bien = true;

                if (code >= total) { biyeccion = false; bien = false; }
                else if (visto[code]) { biyeccion = false; bien = false; }
                else { visto[code] = 1; ++cubiertos; }

                if (!(spiral3dUDDecode(code, origen, eje, izq) == p)) { vuelta = false; bien = false; }
                if ((code == 0) != (x == 0 && y == 0 && z == 0)) { cero = false; bien = false; }

                if (d == 0 && !(spiral3dUDDecode(code, origen, eje, !izq) == p)) {
                    empate = false; bien = false;
                }

                //composición manual: base(m) + ShellIndex
                int64_t cx, cy, cz;
                spiral3dUDToCanonical(x, y, z, eje, cx, cy, cz);
                const uint64_t m = spiral3dUDLayerOf(cx, cy, cz);
                const uint64_t manual = m == 0 ? 0 : spiral3dUDBase(m) + spiral3dUDShellIndex(cx, cy, cz, m);
                if (code != manual) { compone = false; bien = false; }

                if (!bien && ejeMalo < 0) { ejeMalo = e; ladoMalo = izq; bx = x; by = y; bz = z; }
            }
            if (cubiertos != total) {
                biyeccion = false;
                if (ejeMalo < 0) { ejeMalo = e; ladoMalo = izq; }
            }
        }

        char buf[160], ok[160];
        snprintf(buf, sizeof buf, "eje %c, lado %s, punto (%lld,%lld,%lld)",
                 ejeMalo < 0 ? '-' : "xyz"[ejeMalo], ladoMalo ? "izquierdo" : "derecho",
                 (long long)bx, (long long)by, (long long)bz);
        snprintf(ok, sizeof ok, "3 ejes x 2 lados, %llu puntos, capas 0..%lld",
                 (unsigned long long)puntos, (long long)N);
        reportar("Decode(Encode(p)) == p en cada eje y lado", vuelta, vuelta ? ok : buf);
        snprintf(ok, sizeof ok, "6 combinaciones, cada una cubre [0, %llu)", (unsigned long long)total);
        reportar("Encode es biyección sobre [0, count(N)) en cada eje y lado", biyeccion,
                 biyeccion ? ok : buf);
        reportar("código 0 si y solo si el hijo es el padre", cero, cero ? "" : buf);
        reportar("en el plano de empate ambos lados decodifican el mismo punto", empate,
                 empate ? "puntos con desplazamiento 0 en el eje de corte" : buf);
        reportar("Encode == base(m) + ShellIndex (composición manual)", compone, compone ? "" : buf);
    }

    /* Padres lejos del origen (hasta +-1e12) con hijos al azar de todas las
       capas, de la capa máxima y sus esquinas, en cada eje y lado. */
    {
        std::mt19937_64 gen(20260907);
        std::uniform_int_distribution<int64_t> lejos(-1000000000000LL, 1000000000000LL);

        bool ok = true;
        uint64_t mMala = 0, hijos = 0, maxCode = 0;
        int ejeMalo = -1; bool ladoMalo = false;

        //capas: 1 a 20000 al azar, más la máxima con esquinas y celdas al azar
        for (int e = 0; e < 3 && ok; ++e)
        for (int l = 0; l < 2 && ok; ++l) {
            const Spiral3dAxis eje = (Spiral3dAxis)e;
            const bool izq = l == 0;
            for (int i = 0; i < 40000 && ok; ++i) {
                const bool maxima = i >= 20000;
                const uint64_t m = maxima ? SPIRAL3DUD_MAX_LAYER
                                          : std::uniform_int_distribution<uint64_t>(1, 20000)(gen);
                const int64_t sm = (int64_t)m;
                std::uniform_int_distribution<int64_t> libre(-sm, sm), alto(0, sm);
                int64_t cx = libre(gen), cy = libre(gen), cz = alto(gen);
                switch (std::uniform_int_distribution<int>(0, 4)(gen)) {
                case 0:  cz = sm;  break;
                case 1:  cx = sm;  break;
                case 2:  cx = -sm; break;
                case 3:  cy = sm;  break;
                default: cy = -sm; break;
                }
                if (i >= 39992) {                  //las 8 esquinas de la capa máxima
                    const int k = i - 39992;
                    cx = (k & 1) ? sm : -sm;
                    cy = (k & 2) ? sm : -sm;
                    cz = (k & 4) ? sm : 0;
                }

                int64_t dx, dy, dz;
                spiral3dUDFromCanonical(cx, cy, cz, eje, izq, dx, dy, dz);
                const Point3D padre(lejos(gen), lejos(gen), lejos(gen));
                const Point3D hijo(padre.x + dx, padre.y + dy, padre.z + dz);

                const uint64_t code = spiral3dUDEncode(hijo, padre, eje);
                const uint64_t esperado = spiral3dUDBase(m) + spiral3dUDShellIndex(cx, cy, cz, m);
                if (code != esperado || code >= spiral3dUDCount(SPIRAL3DUD_MAX_LAYER)
                    || !(spiral3dUDDecode(code, padre, eje, izq) == hijo)) {
                    ok = false; mMala = m; ejeMalo = e; ladoMalo = izq;
                }
                if (code > maxCode) maxCode = code;
                ++hijos;
            }
        }

        char buf[160];
        if (ok) snprintf(buf, sizeof buf, "%llu hijos, padres hasta +-1e12, código máximo %llu",
                         (unsigned long long)hijos, (unsigned long long)maxCode);
        else    snprintf(buf, sizeof buf, "eje %c, lado %s, capa m=%llu", "xyz"[ejeMalo],
                         ladoMalo ? "izquierdo" : "derecho", (unsigned long long)mMala);
        reportar("ida y vuelta con padres lejanos y hasta la capa máxima", ok, buf);
    }

    if (errores) printf("\nVERIFICACIÓN FALLIDA (%d prueba%s con error)\n\n",
                        errores, errores == 1 ? "" : "s");
    else         printf("\nTODAS LAS PRUEBAS PASARON\n\n");
    return errores ? 1 : 0;
}
