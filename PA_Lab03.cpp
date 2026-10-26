#include <array>
#include <chrono>
#include <clocale>
#include <cstdint>
#include <cstdio>
#include <cwchar>
#include <random>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

constexpr int      K = 1000;     // довжина повідомлення, біти
constexpr int      N = 16;       // степінь полінома (розрядність CRC)
constexpr uint16_t POLY = 0x8005;   // x^16 + x^15 + x^2 + 1
constexpr uint16_t POLY_REFL = 0xA001;   // дзеркальне відображення POLY
constexpr int      BYTES = K / 8;    // повідомлення в байтах
constexpr int      BATCH_SIZE = 1000;     // повідомлень у серії
constexpr int      BATCHES = 100;      // серій (усього BATCH_SIZE*BATCHES дослідів)

static_assert(K % 8 == 0, "K має бути кратним 8, щоб упакувати біти в байти");

using Bits = std::vector<uint8_t>;              
using Table = std::array<uint16_t, 256>;

uint16_t reflect16(uint16_t v) {
    uint16_t r = 0;
    for (int i = 0; i < 16; ++i) {
        r = static_cast<uint16_t>((r << 1) | (v & 1));
        v >>= 1;
    }
    return r;
}

uint8_t reflect8(uint8_t v) {
    uint8_t r = 0;
    for (int i = 0; i < 8; ++i) {
        r = static_cast<uint8_t>((r << 1) | (v & 1));
        v >>= 1;
    }
    return r;
}

void initConsole() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_U16TEXT);
#else
    std::setlocale(LC_ALL, "");
#endif
}

// Упакування бітів у байти: перший біт повідомлення стає СТАРШИМ бітом байта
// (потрібно для прямого табличного алгоритму)
void packMsbFirst(const Bits& bits, uint8_t* out) {
    for (int j = 0; j < BYTES; ++j) {
        uint8_t b = 0;
        for (int i = 0; i < 8; ++i) b = static_cast<uint8_t>((b << 1) | bits[8 * j + i]);
        out[j] = b;
    }
}

// Упакування бітів у байти: перший біт повідомлення стає МОЛОДШИМ бітом байта
// (потрібно для дзеркального табличного алгоритму)
void packLsbFirst(const Bits& bits, uint8_t* out) {
    for (int j = 0; j < BYTES; ++j) {
        uint8_t b = 0;
        for (int i = 7; i >= 0; --i) b = static_cast<uint8_t>((b << 1) | bits[8 * j + i]);
        out[j] = b;
    }
}

// 1. Простий послідовний алгоритм
uint16_t crcSimple(const Bits& msg) {
    uint16_t reg = 0;
    for (uint8_t bit : msg) {                 // біти повідомлення
        const bool out = (reg & 0x8000) != 0;
        reg = static_cast<uint16_t>((reg << 1) | bit);
        if (out) reg ^= POLY;
    }
    for (int i = 0; i < N; ++i) {             // N дописаних нулів
        const bool out = (reg & 0x8000) != 0;
        reg = static_cast<uint16_t>(reg << 1);
        if (out) reg ^= POLY;
    }
    return reg;
}

// 2. Табличний алгоритм (прямий)
Table makeTable() {
    Table t{};
    for (int i = 0; i < 256; ++i) {
        uint16_t r = static_cast<uint16_t>(i << 8);
        for (int k = 0; k < 8; ++k)
            r = (r & 0x8000) ? static_cast<uint16_t>((r << 1) ^ POLY)
            : static_cast<uint16_t>(r << 1);
        t[i] = r;
    }
    return t;
}

uint16_t crcTable(const uint8_t* data, const Table& T) {
    uint16_t reg = 0;
    for (int i = 0; i < BYTES; ++i) {
        const uint8_t top = static_cast<uint8_t>(reg >> 8);
        reg = static_cast<uint16_t>((reg << 8) | data[i]);
        reg ^= T[top];
    }
    for (int i = 0; i < N / 8; ++i) {         // два дописаних нульових байти
        const uint8_t top = static_cast<uint8_t>(reg >> 8);
        reg = static_cast<uint16_t>(reg << 8);
        reg ^= T[top];
    }
    return reg;
}

// 3. Дзеркальний послідовний алгоритм
uint16_t crcSimpleRefl(const Bits& msg) {
    uint16_t reg = 0;
    for (uint8_t bit : msg) {
        const bool x = ((reg & 1) ^ bit) != 0;
        reg >>= 1;
        if (x) reg ^= POLY_REFL;
    }
    return reg;
}

// 4. Дзеркальний табличний алгоритм
Table makeTableRefl() {
    Table t{};
    for (int i = 0; i < 256; ++i) {
        uint16_t c = static_cast<uint16_t>(i);
        for (int k = 0; k < 8; ++k)
            c = (c & 1) ? static_cast<uint16_t>((c >> 1) ^ POLY_REFL)
            : static_cast<uint16_t>(c >> 1);
        t[i] = c;
    }
    return t;
}

uint16_t crcTableRefl(const uint8_t* data, const Table& T) {
    uint16_t reg = 0;
    for (int i = 0; i < BYTES; ++i)
        reg = static_cast<uint16_t>((reg >> 8) ^ T[(reg ^ data[i]) & 0xFF]);
    return reg;
}

// 5. Стандартизований алгоритм CRC-16/ARC:
//    width=16, poly=0x8005, init=0x0000, RefIn=true, RefOut=true,
//    XorOut=0x0000, check("123456789")=0xBB3D.
constexpr uint16_t ARC_INIT = 0x0000;
constexpr bool     ARC_REFIN = true;
constexpr bool     ARC_REFOUT = true;
constexpr uint16_t ARC_XOROUT = 0x0000;

uint16_t crc16Arc(const uint8_t* data, size_t n) {
    uint16_t crc = ARC_INIT;
    for (size_t i = 0; i < n; ++i) {
        uint8_t b = ARC_REFIN ? reflect8(data[i]) : data[i];
        crc ^= static_cast<uint16_t>(b << 8);
        for (int k = 0; k < 8; ++k)
            crc = (crc & 0x8000) ? static_cast<uint16_t>((crc << 1) ^ POLY)
            : static_cast<uint16_t>(crc << 1);
    }
    if (ARC_REFOUT) crc = reflect16(crc);
    return static_cast<uint16_t>(crc ^ ARC_XOROUT);
}

using Clock = std::chrono::steady_clock;

template <class F>
double timeNs(F&& f) {
    const auto t0 = Clock::now();
    f();
    const auto t1 = Clock::now();
    return std::chrono::duration<double, std::nano>(t1 - t0).count();
}

int main() {
    initConsole();
    const Table T = makeTable();
    const Table TR = makeTableRefl();

    std::random_device rd;
    const uint64_t seed = (static_cast<uint64_t>(rd()) << 32) ^ rd();
    std::mt19937_64 rng(seed);
    std::wprintf(L"Поліном: x^16 + x^15 + x^2 + 1 (0x8005), K = %d біт, seed = %llu\n",
        K, static_cast<unsigned long long>(seed));
    std::wprintf(L"Дослідів: %d серій x %d повідомлень = %d\n\n", BATCHES, BATCH_SIZE,
        BATCHES * BATCH_SIZE);

    std::vector<Bits> msgs(BATCH_SIZE, Bits(K));
    std::vector<std::array<uint8_t, BYTES>> msb(BATCH_SIZE), lsb(BATCH_SIZE);
    std::vector<uint16_t> r1(BATCH_SIZE), r2(BATCH_SIZE), r3(BATCH_SIZE), r4(BATCH_SIZE),
        r5(BATCH_SIZE);

    // Сума часів за всіма серіями, нс: індекси 0..4 - алгоритми 1..5
    double total[5] = { 0, 0, 0, 0, 0 };
    long long mismatches = 0;

    for (int batch = -1; batch < BATCHES; ++batch) {     // batch = -1 -- розігрів, не враховується
        for (int j = 0; j < BATCH_SIZE; ++j) {
            for (int i = 0; i < K; i += 64) {
                uint64_t w = rng();
                for (int b = 0; b < 64 && i + b < K; ++b) msgs[j][i + b] = (w >> b) & 1;
            }
            packMsbFirst(msgs[j], msb[j].data());
            packLsbFirst(msgs[j], lsb[j].data());
        }

        double t[5];
        t[0] = timeNs([&] { for (int j = 0; j < BATCH_SIZE; ++j)
            r1[j] = crcSimple(msgs[j]); });
        t[1] = timeNs([&] { for (int j = 0; j < BATCH_SIZE; ++j)
            r2[j] = crcTable(msb[j].data(), T); });
        t[2] = timeNs([&] { for (int j = 0; j < BATCH_SIZE; ++j)
            r3[j] = crcSimpleRefl(msgs[j]); });
        t[3] = timeNs([&] { for (int j = 0; j < BATCH_SIZE; ++j)
            r4[j] = crcTableRefl(lsb[j].data(), TR); });
        t[4] = timeNs([&] { for (int j = 0; j < BATCH_SIZE; ++j)
            r5[j] = crc16Arc(lsb[j].data(), BYTES); });

        if (batch >= 0) {
            for (int a = 0; a < 5; ++a) total[a] += t[a];
            // Перехресна перевірка: усі реалізації мають збігатися
            for (int j = 0; j < BATCH_SIZE; ++j) {
                const bool ok = r1[j] == r2[j] && r3[j] == r4[j] && r3[j] == r5[j] &&
                    r1[j] == reflect16(r3[j]);
                if (!ok) ++mismatches;
            }
        }
    }

    // Кількість ітерацій основного циклу (базових операцій) на одну контрольну суму
    const int iters[5] = { K + N, BYTES + N / 8, K, BYTES, K };
    const wchar_t* names[5] = {
        L"1. Простий послідовний",
        L"2. Табличний",
        L"3. Дзеркальний послідовний",
        L"4. Дзеркальний табличний",
        L"5. Стандартизований CRC-16/ARC" };
    const double experiments = static_cast<double>(BATCHES) * BATCH_SIZE;

    std::wprintf(L"Перехресна перевірка: розбіжностей %lld з %d\n\n", mismatches,
        BATCHES * BATCH_SIZE);
    std::wprintf(L"Ресурсна складність та середній час\n");
    std::wprintf(L"%-40ls %9ls %14ls %12ls %10ls\n", L"Алгоритм", L"ітерацій",
        L"пам'ять табл.", L"нс на CRC", L"нс/ітер.");
    for (int a = 0; a < 5; ++a) {
        const double avg = total[a] / experiments;
        const int it = iters[a];
        const bool usesTable = (a == 1 || a == 3);
        std::wprintf(L"%-40ls %9d %12u Б %12.1f %10.3f\n", names[a], it,
            usesTable ? static_cast<unsigned>(sizeof(Table)) : 0u, avg, avg / it);
    }
    std::wprintf(L"\nПрискорення відносно простого послідовного: табличний x%.1f, "
        L"дзеркальний табличний x%.1f (дзеркальний послідовний x%.1f)\n",
        total[0] / total[1], total[0] / total[3], total[0] / total[2]);
    return mismatches == 0 ? 0 : 2;
}