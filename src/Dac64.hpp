#ifndef DAC64_HPP
#define DAC64_HPP

#include <cstdint>
#include <vector>

/**
 * DAC (Directly Addressable Codes) de 64 bits.
 * Port de ref/cKDtreeUD/dacs.c (Brisaboa, Ladra, Navarro) con valores uint64_t.
 * Descompresión completa, guardar y cargar están en Dac64Extras.hpp.
 */

/**
 * PASO 1: anchos óptimos de bloque por nivel.
 * Misma programación dinámica que optimizationk de dacs.c, pero el histograma
 * se cuenta por cantidad de bits (65 contadores) y no por valor (maxInt + 1).
 * Construcción O(n + 64^2) en tiempo y O(1) en memoria extra.
 */

/** Costo del rank por bit del bitmap, igual que FACT_RANK en dacs.c. */
static const uint64_t DAC64_RANK_COST_FACTOR = 20;

/** Bits de un uint64_t: tope de bits por valor y de niveles. */
static const int DAC64_WORD_BITS = 64;

/**
 * Bits necesarios para representar un valor.
 * @param value valor
 * @return 0 si value = 0, si no 64 - clz(value)
 */
inline int dac64BitLength(uint64_t value) {
    return value ? DAC64_WORD_BITS - __builtin_clzll(value) : 0;
}

/**
 * Anchos de bloque de cada nivel, del menos al más significativo.
 * La suma de los anchos es bits(max), con 1 como mínimo (lista de ceros).
 * @param values valores
 * @param count cantidad de valores
 * @return anchos en bits, uno por nivel
 */
inline std::vector<uint16_t> dac64OptimalWidths(const uint64_t *values, uint64_t count) {
    uint64_t maxValue = 0;
    for (uint64_t i = 0; i < count; ++i)
        if (maxValue < values[i]) maxValue = values[i];

    //nBits de dacs.c; con maxValue = 0 la referencia se indefine, aquí vale 0
    const int highestBit = (maxValue ? dac64BitLength(maxValue) : 1) - 1;

    //countBelow[b] = valores menores que 2^b; el 0 cuenta con 1 bit, así countBelow[0] = 0
    uint64_t countBelow[DAC64_WORD_BITS + 2] = { 0 };
    for (uint64_t i = 0; i < count; ++i) {
        const int bitLength = dac64BitLength(values[i]);
        ++countBelow[bitLength ? bitLength : 1];
    }
    for (int bit = 1; bit <= highestBit + 1; ++bit) countBelow[bit] += countBelow[bit - 1];

    //para guardar los bits firstBit..highestBit: costo, niveles y ancho del primer nivel
    uint64_t bestCost[DAC64_WORD_BITS + 1], levelCount[DAC64_WORD_BITS + 1], firstWidth[DAC64_WORD_BITS + 1];
    for (int firstBit = highestBit; firstBit >= 0; --firstBit) {
        const uint64_t reaching = count - countBelow[firstBit]; //valores con bits desde firstBit

        //cortar el primer nivel antes de nextBit; los empates quedan con el nextBit más alto
        uint64_t splitCost = UINT64_MAX;
        int bestNextBit = highestBit;
        for (int nextBit = highestBit; nextBit >= firstBit + 1; --nextBit) {
            const uint64_t cost = bestCost[nextBit]
                                + reaching * (uint64_t)(nextBit - firstBit + 1) //bloque + bit del bitmap
                                + reaching / DAC64_RANK_COST_FACTOR;
            if (splitCost > cost) { splitCost = cost; bestNextBit = nextBit; }
        }

        const uint64_t singleLevelCost = reaching * (uint64_t)(highestBit - firstBit + 1);
        if (splitCost < singleLevelCost) {
            bestCost[firstBit]   = splitCost;
            levelCount[firstBit] = levelCount[bestNextBit] + 1;
            firstWidth[firstBit] = bestNextBit - firstBit;
        } else {
            bestCost[firstBit]   = singleLevelCost;
            levelCount[firstBit] = 1;
            firstWidth[firstBit] = highestBit - firstBit + 1;
        }
    }

    std::vector<uint16_t> widths;
    for (uint64_t level = 0, bit = 0; level < levelCount[0]; ++level) {
        widths.push_back((uint16_t)firstWidth[bit]);
        bit += firstWidth[bit];
    }
    return widths;
}

/**
 * PASO 2: construcción.
 * El nivel j guarda value - offsets[j] en bloques de widths[j] bits (tablebase de dacs.c),
 * así un valor ocupa solo los niveles que necesita.
 * El bitmap marca con 1 el último bloque de cada valor; el último nivel no lo necesita.
 * Construcción O(n * niveles), con niveles <= 64.
 */

/** Palabras entre muestras del rank: 640 bits, igual que dacs.c (factor 20 en palabras de 32). */
static const uint64_t DAC64_RANK_SAMPLE_WORDS = 10;

/**
 * Máscara con los width bits bajos en 1.
 * @param width ancho, de 0 a 64
 * @return máscara
 */
inline uint64_t dac64LowMask(int width) {
    return width >= DAC64_WORD_BITS ? UINT64_MAX : (1ULL << width) - 1;
}

/**
 * Lee width bits desde position; el bloque puede cruzar dos palabras.
 * @param words arreglo de bits
 * @param position primer bit
 * @param width ancho, de 1 a 64
 * @return bloque leído
 */
inline uint64_t dac64ReadBits(const uint64_t *words, uint64_t position, int width) {
    const uint64_t word = position / DAC64_WORD_BITS;
    const int offset = (int)(position % DAC64_WORD_BITS);
    uint64_t chunk = words[word] >> offset;
    if (offset + width > DAC64_WORD_BITS) chunk |= words[word + 1] << (DAC64_WORD_BITS - offset);
    return chunk & dac64LowMask(width);
}

/**
 * Escribe width bits desde position. El arreglo parte en cero y cada bit se escribe una vez.
 * @param words arreglo de bits
 * @param position primer bit
 * @param width ancho, de 1 a 64
 * @param chunk bloque, se usan sus width bits bajos
 */
inline void dac64WriteBits(uint64_t *words, uint64_t position, int width, uint64_t chunk) {
    const uint64_t word = position / DAC64_WORD_BITS;
    const int offset = (int)(position % DAC64_WORD_BITS);
    chunk &= dac64LowMask(width);
    words[word] |= chunk << offset;
    if (offset + width > DAC64_WORD_BITS) words[word + 1] |= chunk >> (DAC64_WORD_BITS - offset);
}

/** Bitmap con rank O(1): conteo acumulado cada DAC64_RANK_SAMPLE_WORDS palabras más popcount. */
struct Dac64Bitmap {
    uint64_t length = 0;
    std::vector<uint64_t> words;
    std::vector<uint64_t> onesBeforeSample;
};

inline bool dac64BitmapGet(const Dac64Bitmap &bitmap, uint64_t position) {
    return (bitmap.words[position / DAC64_WORD_BITS] >> (position % DAC64_WORD_BITS)) & 1ULL;
}

inline void dac64BitmapSet(Dac64Bitmap &bitmap, uint64_t position) {
    bitmap.words[position / DAC64_WORD_BITS] |= 1ULL << (position % DAC64_WORD_BITS);
}

/**
 * Bitmap de length bits en cero.
 * @param length cantidad de bits
 * @return bitmap sin índice de rank (ver dac64BitmapBuildRank)
 */
inline Dac64Bitmap dac64BitmapCreate(uint64_t length) {
    Dac64Bitmap bitmap;
    bitmap.length = length;
    bitmap.words.assign(length / DAC64_WORD_BITS + 1, 0);
    return bitmap;
}

/** Arma las muestras del rank; se llama después de fijar todos los bits. */
inline void dac64BitmapBuildRank(Dac64Bitmap &bitmap) {
    const uint64_t samples = bitmap.words.size() / DAC64_RANK_SAMPLE_WORDS + 1;
    bitmap.onesBeforeSample.assign(samples, 0);
    uint64_t ones = 0;
    for (uint64_t word = 0; word < bitmap.words.size(); ++word) {
        if (word % DAC64_RANK_SAMPLE_WORDS == 0) bitmap.onesBeforeSample[word / DAC64_RANK_SAMPLE_WORDS] = ones;
        ones += (uint64_t)__builtin_popcountll(bitmap.words[word]);
    }
}

/**
 * Unos en [0, position). Equivale a rank(position - 1) de bitrankw32int.
 * @param bitmap bitmap con rank armado
 * @param position de 0 a length
 * @return cantidad de unos antes de position
 */
inline uint64_t dac64BitmapRank(const Dac64Bitmap &bitmap, uint64_t position) {
    const uint64_t lastWord = position / DAC64_WORD_BITS;
    const uint64_t sample = lastWord / DAC64_RANK_SAMPLE_WORDS;
    uint64_t ones = bitmap.onesBeforeSample[sample];
    for (uint64_t word = sample * DAC64_RANK_SAMPLE_WORDS; word < lastWord; ++word)
        ones += (uint64_t)__builtin_popcountll(bitmap.words[word]);
    const int bitsInLastWord = (int)(position % DAC64_WORD_BITS);
    if (bitsInLastWord) ones += (uint64_t)__builtin_popcountll(bitmap.words[lastWord] & dac64LowMask(bitsInLastWord));
    return ones;
}

/** DAC armado. Entre paréntesis, el campo equivalente de FTRep en dacs.c. */
struct Dac64 {
    uint64_t count = 0;                  //valores guardados (listLength)
    std::vector<uint16_t> widths;        //ancho de bloque por nivel (base_bits)
    std::vector<uint64_t> offsets;       //valor restado en cada nivel (tablebase)
    std::vector<uint64_t> firstEntry;    //primer bloque de cada nivel, niveles + 1 (levelsIndex)
    std::vector<uint64_t> firstBit;      //primer bit de cada nivel en chunks (iniLevel)
    std::vector<uint64_t> stopsBefore;   //unos del bitmap antes de cada nivel (rankLevels)
    uint64_t chunkBits = 0;              //bits usados en chunks (tamCode)
    std::vector<uint64_t> chunks;        //bloques de todos los niveles, seguidos (levels)
    Dac64Bitmap lastChunk;               //1 si el bloque es el último del valor (bS)
};

/**
 * Valor que resta cada nivel: offsets[j] = suma de 2^(w0 + ... + wi) para i < j.
 * El exponente no llega a 64 porque los anchos suman a lo más 64.
 * @param widths ancho de bloque por nivel
 * @return un offset por nivel
 */
inline std::vector<uint64_t> dac64LevelOffsets(const std::vector<uint16_t> &widths) {
    const int levels = (int)widths.size();
    std::vector<uint64_t> offsets(levels);
    uint64_t offset = 0;
    int usedBits = 0;
    for (int level = 0; level < levels; ++level) {
        offsets[level] = offset;
        usedBits += widths[level];
        if (level + 1 < levels) offset += 1ULL << usedBits;
    }
    return offsets;
}

/** Primer bit de cada nivel y bits usados, desde widths y firstEntry. */
inline void dac64LayoutLevels(Dac64 &dac) {
    const int levels = (int)dac.widths.size();
    dac.firstBit.assign(levels, 0);
    dac.chunkBits = 0;
    for (int level = 0; level < levels; ++level) {
        dac.firstBit[level] = dac.chunkBits;
        dac.chunkBits += (dac.firstEntry[level + 1] - dac.firstEntry[level]) * dac.widths[level];
    }
}

/** Arma el rank del bitmap y los unos antes de cada nivel; el bitmap ya debe estar completo. */
inline void dac64IndexStops(Dac64 &dac) {
    dac64BitmapBuildRank(dac.lastChunk);
    const int levels = (int)dac.widths.size();
    dac.stopsBefore.assign(levels, 0);
    for (int level = 0; level < levels; ++level)
        dac.stopsBefore[level] = dac64BitmapRank(dac.lastChunk, dac.firstEntry[level]);
}

/**
 * Construye el DAC de una lista de valores.
 * @param values valores
 * @param count cantidad de valores
 * @return DAC con los valores en el mismo orden
 */
inline Dac64 dac64Build(const uint64_t *values, uint64_t count) {
    Dac64 dac;
    dac.count = count;
    if (count == 0) return dac;

    const std::vector<uint16_t> widths = dac64OptimalWidths(values, count);
    const int plannedLevels = (int)widths.size();
    const std::vector<uint64_t> offsets = dac64LevelOffsets(widths);

    //nivel más alto de cada valor; los offsets crecen, así que se cuentan de abajo hacia arriba
    std::vector<uint64_t> reaching(plannedLevels, 0);
    for (uint64_t i = 0; i < count; ++i)
        for (int level = 0; level < plannedLevels && values[i] >= offsets[level]; ++level)
            ++reaching[level];

    int levels = 0;
    while (levels < plannedLevels && reaching[levels] != 0) ++levels;

    dac.widths.assign(widths.begin(), widths.begin() + levels);
    dac.offsets.assign(offsets.begin(), offsets.begin() + levels);
    dac.firstEntry.assign(levels + 1, 0);
    for (int level = 0; level < levels; ++level)
        dac.firstEntry[level + 1] = dac.firstEntry[level] + reaching[level];
    dac64LayoutLevels(dac);
    dac.chunks.assign(dac.chunkBits / DAC64_WORD_BITS + 1, 0);

    //un bit por bloque de los niveles 0..L-2, más uno final en 1 (bits_BS_len de dacs.c)
    dac.lastChunk = dac64BitmapCreate(dac.firstEntry[levels - 1] + 1);

    std::vector<uint64_t> nextEntry(dac.firstEntry.begin(), dac.firstEntry.end() - 1);
    for (uint64_t i = 0; i < count; ++i) {
        int topLevel = levels - 1;
        while (values[i] < dac.offsets[topLevel]) --topLevel;

        //del bloque menos significativo al más significativo; antes del último el ancho es < 64
        uint64_t remaining = values[i] - dac.offsets[topLevel];
        for (int level = 0; level <= topLevel; ++level) {
            const uint64_t entry = nextEntry[level]++;
            const uint64_t position = dac.firstBit[level] + (entry - dac.firstEntry[level]) * dac.widths[level];
            dac64WriteBits(dac.chunks.data(), position, dac.widths[level], remaining);
            if (level < topLevel) remaining >>= dac.widths[level];
            else if (level < levels - 1) dac64BitmapSet(dac.lastChunk, entry);
        }
    }
    dac64BitmapSet(dac.lastChunk, dac.lastChunk.length - 1);
    dac64IndexStops(dac);
    return dac;
}

/**
 * PASO 3: acceso.
 * Acceso O(1): a lo más un bloque por nivel (niveles <= 64) y un rank O(1) por nivel;
 * nada depende de la cantidad de valores.
 */

/**
 * Valor en la posición index (accessFT(rep, index + 1) de dacs.c).
 * @param dac DAC armado
 * @param index posición, desde 0
 * @return valor guardado
 */
inline uint64_t dac64Access(const Dac64 &dac, uint64_t index) {
    const int lastLevel = (int)dac.widths.size() - 1;
    uint64_t value = 0;
    int shift = 0;
    int level = 0;
    uint64_t entry = index; //posición del valor dentro del nivel
    for (;; ++level) {
        const int width = dac.widths[level];
        value += dac64ReadBits(dac.chunks.data(), dac.firstBit[level] + entry * width, width) << shift;

        const uint64_t bitmapPosition = dac.firstEntry[level] + entry;
        if (level == lastLevel || dac64BitmapGet(dac.lastChunk, bitmapPosition)) break;

        //en el nivel siguiente solo están los valores que no terminaron en este
        entry -= dac64BitmapRank(dac.lastChunk, bitmapPosition) - dac.stopsBefore[level];
        shift += width;
    }
    return value + dac.offsets[level];
}

/**
 * PASO 4: tamaño en memoria.
 */

/**
 * Bits que ocupa el DAC en memoria (memoryUsage de dacs.c, pero en bits).
 * @param dac DAC armado
 * @return bits de todos los campos
 */
inline uint64_t dac64SizeInBits(const Dac64 &dac) {
    const uint64_t scalars = 3; //count, chunkBits y largo del bitmap
    const uint64_t words = scalars + dac.offsets.size() + dac.firstEntry.size() + dac.firstBit.size()
                         + dac.stopsBefore.size() + dac.chunks.size() + dac.lastChunk.words.size()
                         + dac.lastChunk.onesBeforeSample.size();
    return words * DAC64_WORD_BITS + dac.widths.size() * 16;
}

#endif
