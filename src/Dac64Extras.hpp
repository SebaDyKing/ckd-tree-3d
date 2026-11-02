#ifndef DAC64_EXTRAS_HPP
#define DAC64_EXTRAS_HPP

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "Dac64.hpp"

/**
 * Extras del DAC de 64 bits que las búsquedas no usan.
 * La lógica del DAC (anchos, construcción y acceso) está en Dac64.hpp.
 */

/**
 * Descompresión O(n): recorre los niveles en orden, sin rank.
 * La usan las pruebas para compararse con decompressFT de dacs.c.
 */

/**
 * Todos los valores en orden (decompressFT de dacs.c).
 * @param dac DAC armado
 * @return los count valores
 */
inline std::vector<uint64_t> dac64Decompress(const Dac64 &dac) {
    if (dac.count == 0) return {};
    const int lastLevel = (int)dac.widths.size() - 1;
    std::vector<uint64_t> values(dac.count);
    std::vector<uint64_t> nextBit(dac.firstBit);
    std::vector<uint64_t> nextEntry(dac.firstEntry.begin(), dac.firstEntry.end() - 1);
    for (uint64_t i = 0; i < dac.count; ++i) {
        uint64_t value = 0;
        int shift = 0;
        int level = 0;
        for (;; ++level) {
            const int width = dac.widths[level];
            value += dac64ReadBits(dac.chunks.data(), nextBit[level], width) << shift;
            nextBit[level] += width;
            if (level == lastLevel || dac64BitmapGet(dac.lastChunk, nextEntry[level]++)) break;
            shift += width;
        }
        values[i] = value + dac.offsets[level];
    }
    return values;
}

/**
 * Guardar y cargar.
 * El archivo guarda solo count, widths, firstEntry, los bloques y el bitmap.
 * Lo demás se recalcula al cargar con las mismas funciones de la construcción,
 * así un archivo no puede quedar inconsistente. Carga O(n / 64).
 */

/** Encabezado del archivo; distingue este formato del de dacs.c. */
static const char DAC64_FILE_MAGIC[8] = { 'D', 'A', 'C', '6', '4', 'v', '1', '\0' };

/** Tope de valores: count * 64 bits debe caber en 64 bits. */
static const uint64_t DAC64_MAX_COUNT = 1ULL << 57;

template <typename T>
inline bool dac64WriteArray(FILE *file, const T *data, uint64_t length) {
    return length == 0 || fwrite(data, sizeof(T), length, file) == length;
}

template <typename T>
inline bool dac64ReadArray(FILE *file, std::vector<T> &data, uint64_t length) {
    data.assign(length, 0);
    return length == 0 || fread(data.data(), sizeof(T), length, file) == length;
}

/**
 * Escribe el DAC desde la posición actual del archivo; se puede seguir escribiendo después.
 * @param dac DAC armado
 * @param file archivo abierto en modo binario
 * @return true si se escribió completo
 */
inline bool dac64Write(const Dac64 &dac, FILE *file) {
    const uint64_t levels = dac.widths.size();
    bool ok = dac64WriteArray(file, DAC64_FILE_MAGIC, sizeof DAC64_FILE_MAGIC)
           && dac64WriteArray(file, &dac.count, 1);
    if (!ok || dac.count == 0) return ok;
    return dac64WriteArray(file, &levels, 1)
        && dac64WriteArray(file, dac.widths.data(), levels)
        && dac64WriteArray(file, dac.firstEntry.data(), levels + 1)
        && dac64WriteArray(file, dac.chunks.data(), dac.chunks.size())
        && dac64WriteArray(file, dac.lastChunk.words.data(), dac.lastChunk.words.size());
}

/**
 * Lee un DAC escrito con dac64Write desde la posición actual del archivo.
 * Revisa el encabezado y la coherencia de los niveles antes de reservar los bloques.
 * @param file archivo abierto en modo binario
 * @param dac DAC leído; queda vacío si falla
 * @return true si el archivo es válido
 */
inline bool dac64Read(FILE *file, Dac64 &dac) {
    dac = Dac64();
    char magic[sizeof DAC64_FILE_MAGIC];
    uint64_t levels = 0;
    if (fread(magic, 1, sizeof magic, file) != sizeof magic
        || memcmp(magic, DAC64_FILE_MAGIC, sizeof magic) != 0
        || fread(&dac.count, sizeof dac.count, 1, file) != 1) return false;
    if (dac.count == 0) return true;
    if (dac.count > DAC64_MAX_COUNT || fread(&levels, sizeof levels, 1, file) != 1
        || levels == 0 || levels > (uint64_t)DAC64_WORD_BITS) { dac = Dac64(); return false; }

    //anchos de 1 a 64 que suman a lo más 64; niveles que se achican sin vaciarse
    bool valid = dac64ReadArray(file, dac.widths, levels) && dac64ReadArray(file, dac.firstEntry, levels + 1);
    int totalWidth = 0;
    for (uint64_t level = 0; valid && level < levels; ++level) {
        totalWidth += dac.widths[level];
        const uint64_t reaching = dac.firstEntry[level + 1] - dac.firstEntry[level];
        const uint64_t above = level ? dac.firstEntry[level] - dac.firstEntry[level - 1] : dac.count;
        valid = dac.widths[level] >= 1 && totalWidth <= DAC64_WORD_BITS
             && dac.firstEntry[level + 1] > dac.firstEntry[level] && reaching <= above;
    }
    valid = valid && dac.firstEntry[0] == 0 && dac.firstEntry[1] == dac.count;
    if (!valid) { dac = Dac64(); return false; }

    dac.offsets = dac64LevelOffsets(dac.widths);
    dac64LayoutLevels(dac);
    dac.lastChunk.length = dac.firstEntry[levels - 1] + 1;
    valid = dac64ReadArray(file, dac.chunks, dac.chunkBits / DAC64_WORD_BITS + 1)
         && dac64ReadArray(file, dac.lastChunk.words, dac.lastChunk.length / DAC64_WORD_BITS + 1);
    if (!valid) { dac = Dac64(); return false; }

    //en el nivel j terminan los que no siguen al j + 1; si no, el acceso saldría del nivel
    dac64IndexStops(dac);
    valid = dac64BitmapGet(dac.lastChunk, dac.lastChunk.length - 1);
    for (uint64_t level = 0; valid && level + 1 < levels; ++level) {
        const uint64_t ending = (dac.firstEntry[level + 1] - dac.firstEntry[level])
                              - (dac.firstEntry[level + 2] - dac.firstEntry[level + 1]);
        valid = dac.stopsBefore[level + 1] - dac.stopsBefore[level] == ending;
    }
    if (!valid) dac = Dac64();
    return valid;
}

/**
 * Guarda el DAC en un archivo nuevo (saveFT de dacs.c).
 * @param dac DAC armado
 * @param path ruta del archivo
 * @return true si se guardó completo
 */
inline bool dac64Save(const Dac64 &dac, const char *path) {
    FILE *file = fopen(path, "wb");
    if (!file) return false;
    const bool ok = dac64Write(dac, file);
    return fclose(file) == 0 && ok;
}

/**
 * Carga un DAC guardado con dac64Save (loadFT de dacs.c).
 * @param path ruta del archivo
 * @param dac DAC cargado; queda vacío si falla
 * @return true si el archivo es válido
 */
inline bool dac64Load(const char *path, Dac64 &dac) {
    FILE *file = fopen(path, "rb");
    if (!file) { dac = Dac64(); return false; }
    const bool ok = dac64Read(file, dac);
    fclose(file);
    return ok;
}

#endif
