#include <iostream>
#include <numbers>

#include "include/blacklemon.h"
#include "include/global.h"

Signal generateSignal(const PSsk& secretKey) {
    Signal signal;
    auto publicKey = PSpkGen(PSparam, secretKey);
    PSsignal(signal, publicKey, PSparam);
    return signal;
}

Signal generateSignal() {
    auto secretKey = PSskGen(PSparam);
    return generateSignal(secretKey);
}

Signal generateSnakeEye() {
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> dgg(std::numbers::pi);
    Signal signal;
    // signal.a = NativeVector(PSparam.n, PSparam.q);
    // signal.b = NativeVector(PSparam.ell, PSparam.q);
    signal.a = dgg.GenerateVector(PSparam.n, PSparam.q);
    signal.b = dgg.GenerateVector(PSparam.ell, PSparam.q);
    return signal;
}

int main() {
    std::size_t detected{0};
    std::size_t trials{100};
    for (std::size_t i = 0; i < trials; ++i) {
        auto secretKey = PSskGen(PSparam);
        // auto signal = generateSignal();
        auto signal = generateSnakeEye();
        if (PSdetect(signal, secretKey, PSparam))
            ++detected;
    }
    if (detected == 0)
        std::cout << "No snakes today" << std::endl;
    else if (detected == trials)
        std::cout << "Hiss" << std::endl;
    else
        std::cout << "Misdetected " << detected << " of " << trials << std::endl;
    return 0;
}
