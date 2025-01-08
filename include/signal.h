#ifndef SIGNAL_H
#define SIGNAL_H

#include <cstdint>
#include <vector>

#include <openfhe.h>

struct param{
    int n; int q; int ell; int h; double sigma; int r;
    param(int n, int q, int ell, int h, double sigma, int r)
    : n(n), q(q), ell(ell), h(h), sigma(sigma), r(r)
    {}
};

typedef NativeVector PSsk;

struct Signal{
    NativeVector a;
    NativeVector b;
};

typedef Signal PSpk;

PSsk PSskGen(const param& param) 
{
    lbcrypto::TernaryUniformGeneratorImpl<NativeVector> tug;
    PSsk sk = tug.GenerateVector(param.n, param.q, param.h);
    return sk;
}

std::vector<std::vector<uint64_t>> expand(NativeVector a, int n, int q) {
    std::vector<std::vector<uint64_t>> res(n);

    for (int cnt = 0; cnt < n; cnt++) {
        res[cnt].resize(n);
        int ind = 0;

        for (int i = cnt; i >= 0 && ind < n; i--) {
            res[cnt][ind] = a[i].ConvertToInt();
            ind++;
        }

        for (int i = n-1; i > cnt && ind < n; i--) {
            res[cnt][ind] = q - a[i].ConvertToInt();
            ind++;
        }
    }

    return res;
}

NativeVector ringMult(NativeVector a, NativeVector b, int n, int q, int ell) 
{
    NativeVector res = NativeVector(ell, q);

    std::vector<std::vector<uint64_t>> a_expanded = expand(a, n, q);
    for (int i = 0; i < ell; i++) {
        uint64_t temp = 0;
        for (int j = 0; j < n; j++) {
            temp = (temp + a_expanded[i][j] * b[j].ConvertToInt()) % q;
        }
        res[i] = temp;
    }

    return res;
}


PSpk PSpkGen(const param& param, const PSsk& sk) 
{   
    PSpk pk;

    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(param.q);

    pk.a = dug.GenerateVector(param.n);
    pk.b = ringMult(pk.a, sk, param.n, param.q, param.n);

    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> dgg(param.sigma);
    for (int i = 0; i < param.n; i++) {
        pk.b[i].ModAddFastEq(dgg.GenerateInteger(param.q), param.q);
    }

    return pk;
}

void PSsignal(Signal& sig, const PSpk& pk, const param& param) 
{
    lbcrypto::TernaryUniformGeneratorImpl<NativeVector> tug;
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> dgg(param.sigma);

    NativeVector r = tug.GenerateVector(param.n, param.q, param.h);
    sig.a = ringMult(pk.a, r, param.n, param.q, param.n);
    sig.b = ringMult(pk.b, r, param.n, param.q, param.ell);

    for(int i = 0; i < param.n; i++){
        sig.a[i].ModAddFastEq(dgg.GenerateInteger(param.q), param.q);
    }
    for(int i = 0; i < param.ell; i++){
        sig.b[i].ModAddFastEq(dgg.GenerateInteger(param.q), param.q);
    }
}

bool PSdetect(const Signal& sig, const PSsk& sk, const param& param)
{
    NativeVector d{sig.b.ModSub(ringMult(sig.a, sk, param.n, param.q, param.ell))};
    for (std::size_t i = 0; i < d.GetLength(); ++i) {
        int n = d[i].ConvertToInt();
        if (n > param.r && n < param.q - param.r)
            return false;
    }
    return true;
}

#endif
