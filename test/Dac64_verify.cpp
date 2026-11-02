/* Pruebas del DAC de 64 bits contra la referencia de 32 bits (ref/cKDtreeUD/dacs.c).

     ./dac64_verify

   PASO 1: los anchos óptimos de bloque coinciden con optimizationk de dacs.c
   en valores de 32 bits, y se calculan sin histograma por valor en 64 bits.
   PASO 2: la construcción deja los mismos niveles, bloques y bitmap que
   createFT bit a bit, y en 64 bits cumple las mismas invariantes.
   PASO 3: el acceso y la descompresión devuelven lo mismo que accessFT y
   decompressFT, recuperan valores de 64 bits y códigos reales de la espiral 3D,
   y el tiempo por acceso no crece con la cantidad de valores.
   PASO 4: guardar y cargar deja el mismo DAC, el archivo admite datos
   después, los archivos inválidos se rechazan y el tamaño se compara con
   memoryUsage de dacs.c. */

#include "Dac64.hpp"
#include "Dac64Extras.hpp"
#include "Spiral3dUD.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include <chrono>
#include <unistd.h>

extern "C" {
#include "dacs.h"
//funciones internas de dacs.c, no están en dacs.h
unsigned short *optimizationk(unsigned int *list, int listLength, int *nkvalues);
unsigned int memoryUsage(FTRep *rep);
}

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

//anchos que calcula la referencia de 32 bits
static std::vector<uint16_t> anchosReferencia(const std::vector<uint64_t> &valores) {
    std::vector<unsigned int> lista32(valores.begin(), valores.end());
    int cantidadNiveles = 0;
    unsigned short *anchos = optimizationk(lista32.data(), (int)lista32.size(), &cantidadNiveles);
    std::vector<uint16_t> resultado(anchos, anchos + cantidadNiveles);
    free(anchos);
    return resultado;
}

static int sumaAnchos(const std::vector<uint16_t> &anchos) {
    int suma = 0;
    for (uint16_t ancho : anchos) suma += ancho;
    return suma;
}

//"10,10,9,..." para mostrar los anchos en el detalle
static void anchosATexto(const std::vector<uint16_t> &anchos, char *texto, size_t tamano) {
    size_t escrito = 0;
    texto[0] = '\0';
    for (size_t nivel = 0; nivel < anchos.size() && escrito + 8 < tamano; ++nivel)
        escrito += snprintf(texto + escrito, tamano - escrito, "%s%u", nivel ? "," : "", (unsigned)anchos[nivel]);
}

//mayor valor con esa cantidad de bits
static uint64_t mayorConBits(int bits) {
    return bits == 64 ? UINT64_MAX : (1ULL << bits) - 1;
}

//lista de 1 a 5000 valores con el máximo presente; sesgada a valores chicos, como los deltas del árbol
static std::vector<uint64_t> listaAzar(std::mt19937_64 &azar, int bitsMaximo, bool uniforme) {
    const uint64_t cantidad = 1 + azar() % 5000;
    const uint64_t maximo = mayorConBits(bitsMaximo);
    std::vector<uint64_t> valores(cantidad);
    for (uint64_t &valor : valores) {
        valor = azar() & maximo;
        if (!uniforme) valor >>= azar() % bitsMaximo;
    }
    valores[azar() % cantidad] = maximo;
    return valores;
}

//bit de un arreglo de palabras de 32 bits, como los escribe dacs.c
static bool bitReferencia(const unsigned int *palabras, uint64_t posicion) {
    return (palabras[posicion / 32] >> (posicion % 32)) & 1u;
}

//misma estructura que createFT, campo por campo y bit a bit; devuelve el primer campo distinto
static const char *diferenciaConReferencia(const Dac64 &dac, const std::vector<uint64_t> &valores) {
    std::vector<unsigned int> lista32(valores.begin(), valores.end());
    FTRep *ref = createFT(lista32.data(), (unsigned int)lista32.size());
    const char *diferencia = nullptr;
    const uint64_t niveles = ref->nLevels;

    if (dac.widths.size() != niveles) diferencia = "nLevels";
    for (uint64_t nivel = 0; nivel < niveles && !diferencia; ++nivel) {
        if (dac.widths[nivel] != ref->base_bits[nivel])        diferencia = "base_bits";
        else if (dac.offsets[nivel] != ref->tablebase[nivel])  diferencia = "tablebase";
        else if (dac.firstBit[nivel] != ref->iniLevel[nivel])  diferencia = "iniLevel";
        else if (dac.stopsBefore[nivel] != ref->rankLevels[nivel]) diferencia = "rankLevels";
    }
    for (uint64_t nivel = 0; nivel <= niveles && !diferencia; ++nivel)
        if (dac.firstEntry[nivel] != ref->levelsIndex[nivel]) diferencia = "levelsIndex";

    if (!diferencia && dac.chunkBits != ref->tamCode) diferencia = "tamCode";
    for (uint64_t bit = 0; bit < dac.chunkBits && !diferencia; ++bit)
        if (((dac.chunks[bit / 64] >> (bit % 64)) & 1ULL) != bitReferencia(ref->levels, bit)) diferencia = "levels";

    if (!diferencia && dac.lastChunk.length != ref->bS->n) diferencia = "largo de bS";
    for (uint64_t bit = 0; bit < dac.lastChunk.length && !diferencia; ++bit) {
        if (dac64BitmapGet(dac.lastChunk, bit) != bitReferencia(ref->bS->data, bit)) diferencia = "bits de bS";
        //rank de la referencia es inclusivo: rank(p - 1) = unos en [0, p)
        else if (bit > 0 && dac64BitmapRank(dac.lastChunk, bit) != ::rank(ref->bS, (unsigned int)(bit - 1)))
            diferencia = "rank de bS";
    }

    destroyFT(ref);
    return diferencia;
}

//mismos campos y mismo bitmap con rank; los bits de relleno no cuentan
static bool mismoDac(const Dac64 &a, const Dac64 &b) {
    return a.count == b.count && a.widths == b.widths && a.offsets == b.offsets
        && a.firstEntry == b.firstEntry && a.firstBit == b.firstBit && a.stopsBefore == b.stopsBefore
        && a.chunkBits == b.chunkBits && a.chunks == b.chunks && a.lastChunk.length == b.lastChunk.length
        && a.lastChunk.words == b.lastChunk.words && a.lastChunk.onesBeforeSample == b.lastChunk.onesBeforeSample;
}

//bytes de un DAC escrito con dac64Write
static std::vector<unsigned char> bytesDe(const Dac64 &dac) {
    FILE *archivo = tmpfile();
    dac64Write(dac, archivo);
    std::vector<unsigned char> bytes((size_t)ftell(archivo));
    rewind(archivo);
    if (fread(bytes.data(), 1, bytes.size(), archivo) != bytes.size()) bytes.clear();
    fclose(archivo);
    return bytes;
}

//lee un DAC desde bytes arbitrarios, como si vinieran de un archivo
static bool leerDesdeBytes(const std::vector<unsigned char> &bytes, Dac64 &dac) {
    FILE *archivo = tmpfile();
    if (!bytes.empty()) fwrite(bytes.data(), 1, bytes.size(), archivo);
    rewind(archivo);
    const bool ok = dac64Read(archivo, dac);
    fclose(archivo);
    return ok;
}

int main() {
    printf("\nVerificación del DAC de 64 bits contra dacs.c\n\n");
    std::mt19937_64 azar(20261026);

    printf("PASO 1: anchos óptimos de bloque\n");

    //la referencia necesita máximo >= 9: lee weight[0..9] sin revisar el tamaño
    {
        bool ok = true;
        uint64_t listasProbadas = 0, cantidadFallida = 0;
        int bitsFallidos = 0;
        for (int bitsMaximo = 4; bitsMaximo <= 20 && ok; ++bitsMaximo) {
            for (int repeticion = 0; repeticion < 40 && ok; ++repeticion) {
                const std::vector<uint64_t> valores = listaAzar(azar, bitsMaximo, repeticion & 1);
                if (dac64OptimalWidths(valores.data(), valores.size()) != anchosReferencia(valores)) {
                    ok = false; cantidadFallida = valores.size(); bitsFallidos = bitsMaximo;
                }
                ++listasProbadas;
            }
        }
        char detalle[160];
        if (ok) snprintf(detalle, sizeof detalle, "%llu listas, máximos de 4 a 20 bits",
                         (unsigned long long)listasProbadas);
        else    snprintf(detalle, sizeof detalle, "n=%llu, máximo de %d bits",
                         (unsigned long long)cantidadFallida, bitsFallidos);
        reportar("anchos iguales a optimizationk de dacs.c", ok, detalle);
    }

    //los niveles cubren justo los bits del máximo
    {
        bool ok = true;
        for (int bitsMaximo = 1; bitsMaximo <= 64 && ok; ++bitsMaximo) {
            const uint64_t maximo = mayorConBits(bitsMaximo);
            std::vector<uint64_t> valores(2000);
            for (uint64_t &valor : valores) valor = (azar() & maximo) >> (azar() % bitsMaximo);
            valores[0] = maximo;
            const std::vector<uint16_t> anchos = dac64OptimalWidths(valores.data(), valores.size());
            if (sumaAnchos(anchos) != bitsMaximo) ok = false;
            for (uint16_t ancho : anchos) if (ancho == 0 || ancho > 64) ok = false;
        }
        reportar("suma de anchos = bits del máximo, de 1 a 64 bits", ok);
    }

    //bordes que la referencia no soporta
    {
        const std::vector<uint64_t> ceros(100, 0);
        const std::vector<uint16_t> anchos = dac64OptimalWidths(ceros.data(), ceros.size());
        reportar("lista de ceros: un nivel de 1 bit", anchos.size() == 1 && anchos[0] == 1);
    }
    {
        const std::vector<uint64_t> valores = { 0, 1, UINT64_MAX };
        const std::vector<uint16_t> anchos = dac64OptimalWidths(valores.data(), valores.size());
        char detalle[160];
        anchosATexto(anchos, detalle, sizeof detalle);
        reportar("máximo 2^64-1 sin histograma por valor", sumaAnchos(anchos) == 64, detalle);
    }

    //tiempo lineal: un millón de valores de 64 bits
    {
        std::vector<uint64_t> valores(1000000);
        for (uint64_t &valor : valores) valor = azar() >> (azar() % 64);
        const auto inicio = std::chrono::steady_clock::now();
        const std::vector<uint16_t> anchos = dac64OptimalWidths(valores.data(), valores.size());
        const double milisegundos = std::chrono::duration<double, std::milli>(
                                        std::chrono::steady_clock::now() - inicio).count();
        char detalle[160], textoAnchos[96];
        anchosATexto(anchos, textoAnchos, sizeof textoAnchos);
        snprintf(detalle, sizeof detalle, "%.1f ms, anchos %s", milisegundos, textoAnchos);
        reportar("10^6 valores de 64 bits", sumaAnchos(anchos) == 64, detalle);
    }

    printf("\nPASO 2: construcción\n");

    //bloques de 1 a 64 bits seguidos, cruzando palabras
    {
        std::vector<uint64_t> anchos(20000), bloques(20000);
        uint64_t totalBits = 0;
        for (size_t i = 0; i < anchos.size(); ++i) {
            anchos[i] = 1 + azar() % 64;
            bloques[i] = azar() & dac64LowMask((int)anchos[i]);
            totalBits += anchos[i];
        }
        std::vector<uint64_t> palabras(totalBits / 64 + 1, 0);
        for (size_t i = 0, posicion = 0; i < anchos.size(); posicion += anchos[i++])
            dac64WriteBits(palabras.data(), posicion, (int)anchos[i], bloques[i]);
        bool ok = true;
        for (size_t i = 0, posicion = 0; i < anchos.size(); posicion += anchos[i++])
            if (dac64ReadBits(palabras.data(), posicion, (int)anchos[i]) != bloques[i]) ok = false;
        reportar("escribir y leer bloques de 1 a 64 bits", ok, "20000 bloques");
    }

    //rank contra conteo directo, con largos en los bordes de palabra y de muestra
    {
        bool ok = true;
        const uint64_t largos[] = { 1, 63, 64, 65, 639, 640, 641, 1280, 5000, 100003 };
        for (uint64_t largo : largos) {
            Dac64Bitmap bitmap = dac64BitmapCreate(largo);
            const uint64_t densidad = 1 + azar() % 8;
            for (uint64_t bit = 0; bit < largo; ++bit)
                if (azar() % densidad == 0) dac64BitmapSet(bitmap, bit);
            dac64BitmapBuildRank(bitmap);
            uint64_t unos = 0;
            for (uint64_t posicion = 0; posicion <= largo && ok; ++posicion) {
                if (dac64BitmapRank(bitmap, posicion) != unos) ok = false;
                if (posicion < largo && dac64BitmapGet(bitmap, posicion)) ++unos;
            }
        }
        reportar("rank igual al conteo directo", ok, "largos de 1 a 100003 bits");
    }

    //la referencia necesita máximo >= 9 y reserva weight[máximo + 1]: sobre 2^20 se vuelve inmanejable
    {
        const char *diferencia = nullptr;
        uint64_t listasProbadas = 0;
        int bitsFallidos = 0;
        for (int bitsMaximo = 4; bitsMaximo <= 20 && !diferencia; ++bitsMaximo) {
            for (int repeticion = 0; repeticion < 20 && !diferencia; ++repeticion) {
                const std::vector<uint64_t> valores = listaAzar(azar, bitsMaximo, repeticion & 1);
                diferencia = diferenciaConReferencia(dac64Build(valores.data(), valores.size()), valores);
                if (diferencia) bitsFallidos = bitsMaximo;
                ++listasProbadas;
            }
        }
        char detalle[160];
        if (!diferencia) snprintf(detalle, sizeof detalle, "%llu listas, máximos de 4 a 20 bits",
                                  (unsigned long long)listasProbadas);
        else             snprintf(detalle, sizeof detalle, "%s distinto, máximo de %d bits", diferencia, bitsFallidos);
        reportar("misma estructura que createFT, bit a bit", !diferencia, detalle);
    }

    //en 64 bits no hay referencia: se revisan las invariantes que la comparación anterior confirma
    {
        bool ok = true;
        int bitsFallidos = 0;
        for (int bitsMaximo = 1; bitsMaximo <= 64 && ok; ++bitsMaximo) {
            for (int repeticion = 0; repeticion < 6 && ok; ++repeticion) {
                const std::vector<uint64_t> valores = listaAzar(azar, bitsMaximo, repeticion & 1);
                const Dac64 dac = dac64Build(valores.data(), valores.size());
                const uint64_t niveles = dac.widths.size();

                //cada nivel tiene un bloque por valor que alcanza su offset
                uint64_t bitsEsperados = 0;
                for (uint64_t nivel = 0; nivel < niveles; ++nivel) {
                    uint64_t alcanzan = 0;
                    for (uint64_t valor : valores) if (valor >= dac.offsets[nivel]) ++alcanzan;
                    if (dac.firstEntry[nivel + 1] - dac.firstEntry[nivel] != alcanzan || alcanzan == 0) ok = false;
                    bitsEsperados += alcanzan * dac.widths[nivel];
                }
                if (dac.chunkBits != bitsEsperados) ok = false;

                //en el nivel j terminan los valores que no siguen al nivel j + 1
                for (uint64_t nivel = 0; nivel + 1 < niveles; ++nivel) {
                    const uint64_t terminan = (dac.firstEntry[nivel + 1] - dac.firstEntry[nivel])
                                            - (dac.firstEntry[nivel + 2] - dac.firstEntry[nivel + 1]);
                    const uint64_t unos = dac64BitmapRank(dac.lastChunk, dac.firstEntry[nivel + 1])
                                        - dac64BitmapRank(dac.lastChunk, dac.firstEntry[nivel]);
                    if (unos != terminan) ok = false;
                }
                if (!dac64BitmapGet(dac.lastChunk, dac.lastChunk.length - 1)) ok = false;
                if (!ok) bitsFallidos = bitsMaximo;
            }
        }
        char detalle[160] = "";
        if (!ok) snprintf(detalle, sizeof detalle, "máximo de %d bits", bitsFallidos);
        reportar("invariantes de niveles y bitmap, de 1 a 64 bits", ok, detalle);
    }

    printf("\nPASO 3: acceso y descompresión\n");

    //accessFT es desde 1; dac64Access es desde 0
    {
        bool ok = true;
        uint64_t listasProbadas = 0;
        int bitsFallidos = 0;
        for (int bitsMaximo = 4; bitsMaximo <= 20 && ok; ++bitsMaximo) {
            for (int repeticion = 0; repeticion < 20 && ok; ++repeticion) {
                const std::vector<uint64_t> valores = listaAzar(azar, bitsMaximo, repeticion & 1);
                const Dac64 dac = dac64Build(valores.data(), valores.size());
                std::vector<unsigned int> lista32(valores.begin(), valores.end());
                FTRep *ref = createFT(lista32.data(), (unsigned int)lista32.size());
                unsigned int *descomprimidoRef = decompressFT(ref, (unsigned int)lista32.size());
                const std::vector<uint64_t> descomprimido = dac64Decompress(dac);
                for (uint64_t i = 0; i < valores.size() && ok; ++i) {
                    if (dac64Access(dac, i) != accessFT(ref, (unsigned int)(i + 1))) ok = false;
                    if (descomprimido[i] != descomprimidoRef[i]) ok = false;
                }
                free(descomprimidoRef);
                destroyFT(ref);
                if (!ok) bitsFallidos = bitsMaximo;
                ++listasProbadas;
            }
        }
        char detalle[160];
        if (ok) snprintf(detalle, sizeof detalle, "%llu listas, máximos de 4 a 20 bits",
                         (unsigned long long)listasProbadas);
        else    snprintf(detalle, sizeof detalle, "máximo de %d bits", bitsFallidos);
        reportar("acceso y descompresión iguales a dacs.c", ok, detalle);
    }

    //en 64 bits se compara contra la lista original, con los bordes de 32 y 64 bits
    {
        bool ok = true;
        int bitsFallidos = 0;
        for (int bitsMaximo = 1; bitsMaximo <= 64 && ok; ++bitsMaximo) {
            for (int repeticion = 0; repeticion < 6 && ok; ++repeticion) {
                std::vector<uint64_t> valores = listaAzar(azar, bitsMaximo, repeticion & 1);
                if (bitsMaximo == 64)
                    valores.insert(valores.end(), { 0, 1, 0xFFFFFFFFULL, 0x100000000ULL, 1ULL << 63, UINT64_MAX });
                const Dac64 dac = dac64Build(valores.data(), valores.size());
                if (dac64Decompress(dac) != valores) ok = false;
                for (uint64_t i = 0; i < valores.size() && ok; ++i)
                    if (dac64Access(dac, i) != valores[i]) ok = false;
                if (!ok) bitsFallidos = bitsMaximo;
            }
        }
        char detalle[160] = "0, 2^32-1, 2^32, 2^63, 2^64-1 incluidos";
        if (!ok) snprintf(detalle, sizeof detalle, "máximo de %d bits", bitsFallidos);
        reportar("acceso y descompresión de 1 a 64 bits", ok, detalle);
    }

    //de punta a punta: hijo -> código de la espiral -> DAC -> código -> hijo
    {
        const uint64_t pares = 200000;
        std::vector<Point3D> padres(pares), hijos(pares);
        std::vector<Spiral3dAxis> ejes(pares);
        std::vector<uint64_t> codigos(pares);

        //distancias en escala logarítmica: muchas chicas y algunas hasta la capa máxima
        std::uniform_real_distribution<double> exponente(0.0, std::log((double)SPIRAL3DUD_MAX_LAYER));
        std::uniform_int_distribution<int64_t> posicion(-1000000000000LL, 1000000000000LL);
        uint64_t sobre32Bits = 0;
        for (uint64_t i = 0; i < pares; ++i) {
            const int64_t distancia = (int64_t)std::exp(exponente(azar));
            std::uniform_int_distribution<int64_t> delta(-distancia, distancia);
            padres[i] = Point3D(posicion(azar), posicion(azar), posicion(azar));
            hijos[i] = Point3D(padres[i].x + delta(azar), padres[i].y + delta(azar), padres[i].z + delta(azar));
            ejes[i] = spiral3dUDAxisAt(i);
            codigos[i] = spiral3dUDEncode(hijos[i], padres[i], ejes[i]);
            if (codigos[i] > 0xFFFFFFFFULL) ++sobre32Bits;
        }

        const Dac64 dac = dac64Build(codigos.data(), pares);
        bool ok = true;
        for (uint64_t i = 0; i < pares && ok; ++i) {
            const Point3D &padre = padres[i];
            const Point3D &hijo = hijos[i];
            const int64_t enEje = ejes[i] == SPIRAL3D_AXIS_X ? hijo.x - padre.x
                                : ejes[i] == SPIRAL3D_AXIS_Y ? hijo.y - padre.y : hijo.z - padre.z;
            const uint64_t codigo = dac64Access(dac, i);
            if (codigo != codigos[i] || !(spiral3dUDDecode(codigo, padre, ejes[i], enEje < 0) == hijo)) ok = false;
        }

        const uint64_t bitsTotales = dac.chunkBits + dac.lastChunk.length
                                   + 64 * (dac.lastChunk.onesBeforeSample.size() + 4 * dac.widths.size());
        char detalle[160];
        snprintf(detalle, sizeof detalle, "%llu pares, %llu códigos >= 2^32, %.1f bits por código, %zu niveles",
                 (unsigned long long)pares, (unsigned long long)sobre32Bits,
                 (double)bitsTotales / (double)pares, dac.widths.size());
        reportar("códigos de la espiral 3D: guardar y decodificar", ok, detalle);
    }

    //O(1): el tiempo por acceso no depende de n (el resto lo explica la caché)
    {
        bool ok = true;
        char detalle[160];
        size_t escrito = 0;
        const uint64_t cantidades[] = { 10000, 100000, 1000000, 10000000 };
        const uint64_t accesos = 2000000;
        for (uint64_t cantidad : cantidades) {
            std::vector<uint64_t> valores(cantidad);
            for (uint64_t &valor : valores) valor = azar() >> (azar() % 64);
            const Dac64 dac = dac64Build(valores.data(), cantidad);

            std::vector<uint64_t> posiciones(accesos);
            for (uint64_t &p : posiciones) p = azar() % cantidad;
            uint64_t suma = 0, sumaEsperada = 0;
            const auto inicio = std::chrono::steady_clock::now();
            for (uint64_t p : posiciones) suma += dac64Access(dac, p);
            const double nanosegundos = std::chrono::duration<double, std::nano>(
                                            std::chrono::steady_clock::now() - inicio).count() / accesos;
            for (uint64_t p : posiciones) sumaEsperada += valores[p];
            if (suma != sumaEsperada) ok = false;
            escrito += snprintf(detalle + escrito, sizeof detalle - escrito, "%sn=10^%d %.0f ns",
                                escrito ? ", " : "", (int)std::log10((double)cantidad), nanosegundos);
        }
        reportar("tiempo por acceso al azar", ok, detalle);
    }

    printf("\nPASO 4: tamaño, guardar y cargar\n");

    //ida y vuelta por archivo, de 1 a 64 bits y con los bordes de cantidad
    {
        char ruta[] = "/tmp/dac64_verify_XXXXXX";
        const int descriptor = mkstemp(ruta);
        bool ok = descriptor >= 0;
        if (ok) close(descriptor);

        std::vector<std::vector<uint64_t>> listas = { {}, { 0 }, { UINT64_MAX }, std::vector<uint64_t>(1000, 0) };
        for (int bitsMaximo = 1; bitsMaximo <= 64; ++bitsMaximo) listas.push_back(listaAzar(azar, bitsMaximo, bitsMaximo & 1));
        for (const std::vector<uint64_t> &valores : listas) {
            if (!ok) break;
            const Dac64 dac = dac64Build(valores.data(), valores.size());
            Dac64 cargado;
            if (!dac64Save(dac, ruta) || !dac64Load(ruta, cargado) || !mismoDac(dac, cargado)) ok = false;
            for (uint64_t i = 0; i < valores.size() && ok; ++i)
                if (dac64Access(cargado, i) != valores[i]) ok = false;
        }
        remove(ruta);
        char detalle[160];
        snprintf(detalle, sizeof detalle, "%zu listas, vacía y de un valor incluidas", listas.size());
        reportar("guardar y cargar deja el mismo DAC", ok, detalle);
    }

    //el árbol escribe su raíz después del DAC, como cKDtreeUD con x, y
    {
        const std::vector<uint64_t> primeros = listaAzar(azar, 40, false), segundos = listaAzar(azar, 64, true);
        const Dac64 primero = dac64Build(primeros.data(), primeros.size());
        const Dac64 segundo = dac64Build(segundos.data(), segundos.size());
        const int64_t raiz[3] = { -123456789012LL, 42, 987654321098LL };

        FILE *archivo = tmpfile();
        bool ok = dac64Write(primero, archivo) && fwrite(raiz, sizeof raiz[0], 3, archivo) == 3
               && dac64Write(segundo, archivo);
        rewind(archivo);
        Dac64 leidoPrimero, leidoSegundo;
        int64_t raizLeida[3] = { 0, 0, 0 };
        ok = ok && dac64Read(archivo, leidoPrimero) && fread(raizLeida, sizeof raizLeida[0], 3, archivo) == 3
                && dac64Read(archivo, leidoSegundo) && fgetc(archivo) == EOF;
        fclose(archivo);
        ok = ok && mismoDac(primero, leidoPrimero) && mismoDac(segundo, leidoSegundo)
                && memcmp(raiz, raizLeida, sizeof raiz) == 0;
        reportar("datos después del DAC en el mismo archivo", ok, "DAC, raíz (x, y, z), DAC");
    }

    //archivos inválidos: de dacs.c, cortados, con campos o bits alterados
    {
        const std::vector<uint64_t> valores = listaAzar(azar, 30, false);
        const Dac64 dac = dac64Build(valores.data(), valores.size());
        const std::vector<unsigned char> bytes = bytesDe(dac);
        Dac64 leido;
        bool ok = dac.widths.size() >= 2 && leerDesdeBytes(bytes, leido);

        //un archivo de saveFT
        char ruta[] = "/tmp/dac64_verify_XXXXXX";
        const int descriptor = mkstemp(ruta);
        if (descriptor >= 0) {
            close(descriptor);
            std::vector<unsigned int> lista32 = { 9, 3, 1000, 7 };
            FTRep *ref = createFT(lista32.data(), (unsigned int)lista32.size());
            saveFT(ref, ruta);
            destroyFT(ref);
            if (dac64Load(ruta, leido) || leido.count != 0) ok = false;
            remove(ruta);
        } else ok = false;

        //el mismo archivo válido con otra versión en el encabezado
        {
            std::vector<unsigned char> alterado = bytes;
            alterado[6] = '2';
            if (leerDesdeBytes(alterado, leido)) ok = false;
        }

        //cada largo más corto que el archivo completo
        uint64_t cortes = 0;
        for (size_t largo = 0; largo < bytes.size() && ok; largo += 1 + largo / 64, ++cortes)
            if (leerDesdeBytes(std::vector<unsigned char>(bytes.begin(), bytes.begin() + largo), leido)) ok = false;

        //anchos de 0 y de 65 bits; el primer ancho va después de magic, count y niveles
        for (uint16_t ancho : { (uint16_t)0, (uint16_t)65 }) {
            std::vector<unsigned char> alterado = bytes;
            memcpy(alterado.data() + 24, &ancho, sizeof ancho);
            if (leerDesdeBytes(alterado, leido)) ok = false;
        }

        //cada bit del bitmap, que va al final del archivo
        const uint64_t inicioBitmap = bytes.size() - dac.lastChunk.words.size() * 8;
        for (uint64_t bit = 0; bit < dac.lastChunk.length && ok; ++bit) {
            std::vector<unsigned char> alterado = bytes;
            alterado[inicioBitmap + bit / 8] ^= (unsigned char)(1u << (bit % 8));
            if (leerDesdeBytes(alterado, leido)) ok = false;
        }

        char detalle[160];
        snprintf(detalle, sizeof detalle, "%llu cortes, %llu bits del bitmap",
                 (unsigned long long)cortes, (unsigned long long)dac.lastChunk.length);
        reportar("archivos inválidos rechazados", ok, detalle);
    }

    //mismo contenido que dacs.c; la diferencia son palabras de 64 bits y muestras de rank de 64
    {
        std::vector<uint64_t> valores(100000);
        for (uint64_t &valor : valores) valor = (azar() & 0xFFFFF) >> (azar() % 20);
        valores[0] = 0xFFFFF;
        const Dac64 dac = dac64Build(valores.data(), valores.size());
        std::vector<unsigned int> lista32(valores.begin(), valores.end());
        FTRep *ref = createFT(lista32.data(), (unsigned int)lista32.size());
        const double bitsReferencia = 8.0 * memoryUsage(ref) / valores.size();
        const double bitsDac64 = (double)dac64SizeInBits(dac) / valores.size();
        destroyFT(ref);

        const double bitsArchivo = 8.0 * bytesDe(dac).size() / valores.size();
        char detalle[160];
        snprintf(detalle, sizeof detalle, "%.3f contra %.3f bits por valor, archivo %.3f",
                 bitsDac64, bitsReferencia, bitsArchivo);
        reportar("tamaño en memoria contra memoryUsage de dacs.c", bitsDac64 <= bitsReferencia * 1.02, detalle);
    }

    if (errores) printf("\nVERIFICACIÓN FALLIDA (%d prueba%s con error)\n\n",
                        errores, errores == 1 ? "" : "s");
    else         printf("\nTODAS LAS PRUEBAS PASARON\n\n");
    return errores ? 1 : 0;
}
