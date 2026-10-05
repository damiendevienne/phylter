// Standalone lifetime/bounds check for src/mc.cpp; compile with ASan/LSan.
#include <cmath>
#include <vector>
#ifdef PHYLTER_TRACK_ARRAYS
#include <cstdlib>
#include <cstdio>
#include <new>
static long live_arrays = 0;
void* operator new[](std::size_t size) {
    void* p = std::malloc(size ? size : 1);
    if (!p) throw std::bad_alloc();
    ++live_arrays;
    return p;
}
void operator delete[](void* p) noexcept {
    if (p) { --live_arrays; std::free(p); }
}
void operator delete[](void* p, std::size_t) noexcept { operator delete[](p); }
#endif
extern double mlmccN(const double[], int, int);
int main() {
    for (int repeat = 0; repeat < 20; ++repeat) {
        for (int n : {1, 2, 3, 20, 99, 100, 101, 200, 1001}) {
            std::vector<double> x(n, 1.0);
            for (int reflect : {0, 1}) {
                if (!std::isfinite(mlmccN(x.data(), n, reflect))) return 1;
            }
            for (int i = 0; i < n; ++i) x[i] = std::sin(i * 1.234) + i * 0.001;
            for (int reflect : {0, 1}) {
                if (!std::isfinite(mlmccN(x.data(), n, reflect))) return 2;
            }
        }
    }
#ifdef PHYLTER_TRACK_ARRAYS
    std::printf("Unreleased raw arrays: %ld\n", live_arrays);
    return live_arrays == 0 ? 0 : 3;
#endif
}
