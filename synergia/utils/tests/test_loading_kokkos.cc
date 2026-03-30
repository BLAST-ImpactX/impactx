#include <dlfcn.h>
#include <stdio.h>

#include <Kokkos_Core.hpp>

int
main()
{
    void* lib_kokkoscore = dlopen("libkokkoscore.so", RTLD_NOW);
    if (!lib_kokkoscore) {
        printf("Failed to load libkokkoscore.so\n");
        return 1;
    }

    void* p_initialize =
        dlsym(lib_kokkoscore,
              "__ZN6Kokkos10initializeERKNS_22InitializationSettingsE");

    if (!p_initialize) {
        printf("Failed to find Kokkos::initialize\n");
        dlclose(lib_kokkoscore);
        return 2;
    }

    void* p_finalize = dlsym(lib_kokkoscore, "__ZN6Kokkos8finalizeEv");
    if (!p_finalize) {
        printf("Failed to find Kokkos::finalize\n");
        dlclose(lib_kokkoscore);
        return 3;
    }

    // Try calling the Kokkos::initialize() function used in the Python code
    // from a C++ program that is not linked to the library. Make sure to call
    // Kokkos::finalize() to clean up after ourselves.
    using initialization_fun_t = void (*)(Kokkos::InitializationSettings);
    using finalization_fun_t = void (*)();

    auto call_ctor = reinterpret_cast<initialization_fun_t>(p_initialize);
    auto call_finalize = reinterpret_cast<finalization_fun_t>(p_finalize);

    (*call_ctor)(Kokkos::InitializationSettings());
    (*call_finalize)();

    dlclose(lib_kokkoscore);
    printf("dlclose(lib_kokkoscore) completed\n");
    return 0;
}
