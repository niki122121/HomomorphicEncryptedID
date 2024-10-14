#include <cmath>

#include "ipcl/ipcl.hpp"

#include <openssl/opensslv.h>
#if OPENSSL_VERSION_NUMBER < 0x10100000L
#error "OpenSSL version 1.1.0 or later is required.
#endif

//CODE STOPPER
#if OPENSSL_VERSION_NUMBER > 0x30000020L
#error "The operation EC_POINTs_mul is used in this code and for future versions of OpenSSL it could be removed. Make sure your version still has this operation and remove this code stopper if so"
#endif


#include <openssl/ec.h>
#include <openssl/evp.h>
#include <openssl/bn.h>

#include <iostream>
#include <random>
#include <chrono>
using namespace std::chrono;


#define EMB_LENG_0 128
#define EMB_LENG_1 256
#define EMB_LENG_2 512
#define EMB_LENG_3 4096

#define BIT_PRECISION_0 8
#define BIT_PRECISION_1 10
#define BIT_PRECISION_2 12

typedef unsigned char byte;
typedef unsigned int uint;
typedef long long slong;
typedef unsigned long long unslong;

//seed for random generation
uint resultSeed;

//EC global params
EC_GROUP* curve;
BN_CTX* ctx;
BIGNUM* curve_p;
BIGNUM* curve_a;
BIGNUM* curve_b;
BIGNUM* orderMinus1;
BIGNUM* aux_x_bn;
byte bnBuffer[32];

//timers
slong paillierCompTimer = 0, paillierDecTimer = 0;
slong ecegFindyTimer = 0, ecegCompTimer = 0, ecegDecTimer = 0, ecegUnmapTimer = 0;
double paillierPrecision = 0, ecegPrecision = 0;
std::chrono::high_resolution_clock::time_point start, stop;
std::chrono::nanoseconds timer_duration;

template <const int leng, typename T>
void printVec(T* vec){
    for(int i=0; i<leng; i++){
        std::cout << vec[i] << ", ";
    }
    std::cout << std::endl;
}

inline unslong hashECP(EC_POINT* point){
    EC_POINT_get_affine_coordinates(curve, point, aux_x_bn, NULL, ctx);
    BN_bn2bin(aux_x_bn, bnBuffer);
    unslong* pointHash_p = reinterpret_cast<unslong*>(bnBuffer);
    return *pointHash_p;
}
//test the chosen hash function for storing ec points into unordered_map
//in this case, we take the first 8 bytes of the x coordenate
template <const int bitPrecision>
void testBSGSHashFunction() {
    std::cout  << "TESTING HASH FUNCTION WITH ALL CASES" << std::endl;
    std::unordered_map<unslong, uint> babySteps;

    BIGNUM* i_bn = BN_new();
    EC_POINT* auxP = EC_POINT_new(curve);
    for(uint i=1; i<(1<<bitPrecision); i++){
        BN_set_word(i_bn, i);
        EC_POINT_mul(curve, auxP, i_bn, NULL, NULL, ctx);
        unslong ECP_hash = hashECP(auxP);
        if (babySteps.find(ECP_hash) == babySteps.end()) {
            babySteps[ECP_hash] = (uint)123;
        } else {
            std::cout << "BABY STEP COLLISION IN HASH FUNCTION (use a different message range or hash function)" << std::endl;
            std::cerr << "ERROR: BABY STEP COLLISION IN HASH FUNCTION (use a different message range or hash function)" << std::endl;
            BN_free(i_bn);
            EC_POINT_free(auxP);

            exit(EXIT_FAILURE);
        }
    }
    BN_free(i_bn);
    EC_POINT_free(auxP);
    std::cout  << "TEST SUCCESS\n";
    std::cout << "------------------------------------------------------------\n\n";
}
template <const int bsgs_n>
void computeBabySteps(std::unordered_map<unslong, uint> &baby_steps){
    BIGNUM* j_bn = BN_new();
    EC_POINT* P_j = EC_POINT_new(curve);
    for (uint j = 0; j <= bsgs_n; j++) {
        BN_set_word(j_bn, j);
        EC_POINT_mul(curve, P_j, j_bn, NULL, NULL, ctx);
        baby_steps[hashECP(P_j)] = j;
    }
    BN_free(j_bn);
    EC_POINT_free(P_j);
}
template <const int bsgs_n>
inline uint computeGiantSteps(EC_POINT* point, EC_POINT* neg_nG, std::unordered_map<unslong, uint> &baby_steps){
    BIGNUM* i_bn = BN_new();
    EC_POINT* P_i = EC_POINT_new(curve);
    for (uint i = 0; i <= bsgs_n; i++){
        BN_set_word(i_bn, i);
        EC_POINT_mul(curve, P_i, NULL, neg_nG, i_bn, ctx);
        EC_POINT_add(curve, P_i, P_i, point, ctx);
        auto it = baby_steps.find(hashECP(P_i));
        if(it != baby_steps.end()){
            BN_free(i_bn);
            EC_POINT_free(P_i);
            return i*bsgs_n + it->second;
        }
    }
    std::cout << "GIANT STEP NOT FOUND,   SEED: " << resultSeed << "\n";
    std::cout << "ATTEMPTING BRUTE FORCE SEARCH OF " << (bsgs_n*bsgs_n)<< " DEPTH \n";
    for(uint i=0; i<= (bsgs_n*bsgs_n); i++){
        BN_set_word(i_bn, i);
        EC_POINT_mul(curve, P_i, i_bn, NULL, NULL, ctx);
        int foundPoint = EC_POINT_cmp(curve, P_i, point, ctx);
        if(foundPoint == 0){
            std::cout << "FOUND POINT:  " << i << "\n";
            return i;
        }
    }
    std::cout << "DIDN'T FIND POINT, RETURNING 0\n";
    BN_free(i_bn);
    EC_POINT_free(P_i);
    return 0;
}

//from coordinate x, get coordinate y decompressing point
//This code can be optimized for curves where p mod 4 = 3  (no need for modular square root)
inline void decompress(BIGNUM* y, BIGNUM* x){
    BIGNUM *y2 = BN_new();
    BIGNUM *tmp = BN_new();
    BIGNUM *tmp2 = BN_new();

    //compute y^2 = x^3 + ax + b (mod p)
    BN_mod_sqr(tmp, x, curve_p, ctx);
    BN_mod_mul(tmp, tmp, x, curve_p, ctx);
    BN_mod_mul(tmp2, curve_a, x, curve_p, ctx);
    BN_mod_add(y2, tmp, tmp2, curve_p, ctx);
    BN_mod_add(y2, y2, curve_b, curve_p, ctx);

    if (!BN_mod_sqrt(y, y2, curve_p, ctx)) {
        std::cout << "INVALID EC POINT, NO Y COORDINATE FOUND" << std::endl;
        std::cerr << "ERROR: INVALID EC POINT, NO Y COORDINATE FOUND" << std::endl;
        exit(EXIT_FAILURE);
        BN_free(y2);
        BN_free(tmp);
        BN_free(tmp2);
    }
    BN_free(y2);
    BN_free(tmp);
    BN_free(tmp2);
}


template <const int bitPrecision>
inline uint floatToUint(float val) {
    return static_cast<uint>(std::round(static_cast<double>(val)*(static_cast<unslong>(1) << bitPrecision)));;
}
template <const int adjustedBitPrecision>
inline float uintToFloat(uint val) {
    return static_cast<float>(static_cast<double>(val)/(static_cast<unslong>(1) << adjustedBitPrecision));
}

//cosine distance for NORMALIZED embeddings
template <const int embSize>
float cosDist_normalized(float embA[embSize], float embB[embSize]) {
    float cDist = 0;
    for (int i = 0; i < embSize; i++) {
        cDist += (embA[i] * embB[i]);
    }
    return cDist;
}
//euclidean distance for NORMALIZED embeddings
float eucDist_normalized(float cDist) {
    return 2-2*cDist;
}


//paillier cosine distance for NORMALIZED embeddings
template <const int embSize, const int bitPrecision>
inline float pallierCosDis_normalized(float embA[embSize], float embB[embSize], const ipcl::KeyPair& keyPaillier) {
    ipcl::PlainText embA_plaintext;
    ipcl::PlainText embB_plaintext;
    ipcl::CipherText embA_ciphertext;
    ipcl::CipherText cDist_ciphertext = keyPaillier.pub_key.encrypt(ipcl::PlainText(0));
    for(int i=0; i<embSize; i++){
        uint embA_uint = floatToUint<bitPrecision>(embA[i]);
        uint embB_uint = floatToUint<bitPrecision>(embB[i]);

        embA_plaintext = ipcl::PlainText(embA_uint);
        embA_ciphertext = keyPaillier.pub_key.encrypt(embA_plaintext);

        embB_plaintext = ipcl::PlainText(embB_uint);

        start = std::chrono::high_resolution_clock::now();
        cDist_ciphertext = cDist_ciphertext + (embA_ciphertext * embB_plaintext);
        stop = std::chrono::high_resolution_clock::now();
        timer_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        paillierCompTimer+= timer_duration.count();
    }
    start = std::chrono::high_resolution_clock::now();
    uint result_uint = keyPaillier.priv_key.decrypt(cDist_ciphertext).getElementVec(0)[0];
    float result_float = uintToFloat<2*bitPrecision>(result_uint);
    stop = std::chrono::high_resolution_clock::now();
    timer_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    paillierDecTimer+= timer_duration.count();
    return result_float;
}

//EC elGamal cosine distance for NORMALIZED embeddings
//For EC elGamal, unlike Paillier instead of ecnrypting, multiplying and adding as we go, we put all the points on an array to utilize the openssl built in Shamir's Trick
template <const int embSize, const int bitPrecision, const int bsgs_n>
inline float ecegCosDis_normalized(float embA[embSize], float embB[embSize], BIGNUM* secKeyECEG, EC_POINT* pubKeyECEG, EC_POINT* neg_nG, std::unordered_map<unslong, uint> &baby_steps) {
    BIGNUM* embA_bigNum = BN_new();
    BIGNUM* ephemeralKey = BN_new();
    BIGNUM** embB_bigNums = new BIGNUM*[embSize];
    EC_POINT** embA_ciphertexts_1 = new EC_POINT*[embSize];
    EC_POINT** embA_ciphertexts_2 = new EC_POINT*[embSize];
    BIGNUM* comp_x_1 = BN_new();
    BIGNUM* decomp_y_1 = BN_new();
    BIGNUM* comp_x_2 = BN_new();
    BIGNUM* decomp_y_2 = BN_new();
    for(int i=0; i<embSize; i++){
        uint embA_uint = floatToUint<bitPrecision>(embA[i]);
        uint embB_uint = floatToUint<bitPrecision>(embB[i]);

        BN_set_word(embA_bigNum, embA_uint);
        BN_rand_range(ephemeralKey, orderMinus1); //rand number in range [0, order-2]
        BN_add_word(ephemeralKey, 1);             //rand number in range [1, order-1]
        embA_ciphertexts_1[i] = EC_POINT_new(curve);
        embA_ciphertexts_2[i] = EC_POINT_new(curve);
        EC_POINT_mul(curve, embA_ciphertexts_1[i], ephemeralKey, NULL, NULL, ctx);
        EC_POINT_mul(curve, embA_ciphertexts_2[i], embA_bigNum, pubKeyECEG, ephemeralKey, ctx); // map and encrypt => m*G + k*Q (from paper notation)

        embB_bigNums[i] = BN_new();
        BN_set_word(embB_bigNums[i], embB_uint);

        //decompression (this part is a bit redundant here but serves to measure decompression times in a real application)
        EC_POINT_get_affine_coordinates(curve, embA_ciphertexts_1[i], comp_x_1, NULL, ctx);
        EC_POINT_get_affine_coordinates(curve, embA_ciphertexts_2[i], comp_x_2, NULL, ctx);
        start = std::chrono::high_resolution_clock::now();
        decompress(decomp_y_1, comp_x_1);
        decompress(decomp_y_2, comp_x_2);
        stop = std::chrono::high_resolution_clock::now();
        timer_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        ecegFindyTimer+= timer_duration.count();
    }
    //computation
    start = std::chrono::high_resolution_clock::now();
    EC_POINT* ciphertext_result_1 = EC_POINT_new(curve);
    EC_POINT* ciphertext_result_2 = EC_POINT_new(curve);
    EC_POINTs_mul(curve, ciphertext_result_1, NULL, embSize, const_cast<const EC_POINT**>(embA_ciphertexts_1), const_cast<const BIGNUM**>(embB_bigNums), ctx);
    EC_POINTs_mul(curve, ciphertext_result_2, NULL, embSize, const_cast<const EC_POINT**>(embA_ciphertexts_2), const_cast<const BIGNUM**>(embB_bigNums), ctx);
    stop = std::chrono::high_resolution_clock::now();
    timer_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    ecegCompTimer+= timer_duration.count();
    //decryption
    start = std::chrono::high_resolution_clock::now();
    EC_POINT* map_result = EC_POINT_new(curve);
    EC_POINT_mul(curve, map_result, NULL, ciphertext_result_1, secKeyECEG, ctx);
    EC_POINT_invert(curve, map_result, ctx);
    EC_POINT_add(curve, map_result, map_result, ciphertext_result_2, ctx);
    stop = std::chrono::high_resolution_clock::now();
    timer_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    ecegDecTimer+= timer_duration.count();

    //unmapping
    start = std::chrono::high_resolution_clock::now();
    uint result_uint = computeGiantSteps<bsgs_n>(map_result, neg_nG, baby_steps);
    float result_float = uintToFloat<2*bitPrecision>(result_uint);
    stop = std::chrono::high_resolution_clock::now();
    timer_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    ecegUnmapTimer+= timer_duration.count();

    //cleanup
    BN_free(embA_bigNum);
    BN_free(ephemeralKey);
    for(int i=0; i<embSize; i++){
        BN_free(embB_bigNums[i]);
        EC_POINT_free(embA_ciphertexts_1[i]);
        EC_POINT_free(embA_ciphertexts_2[i]);
    }
    delete[] embB_bigNums;
    delete[] embA_ciphertexts_1;
    delete[] embA_ciphertexts_2;
    EC_POINT_free(ciphertext_result_1);
    EC_POINT_free(ciphertext_result_2);
    EC_POINT_free(map_result);
    BN_free(comp_x_1);
    BN_free(comp_x_2);
    BN_free(decomp_y_1);
    BN_free(decomp_y_2);

    return result_float;
}

//this code test timing based on embSize and bitPrecision
template <const int embSize, const int bitPrecision, const int totalIterations>
void timeTest(std::mt19937 &gen, std::uniform_real_distribution<float> &unif, const ipcl::KeyPair& keyPaillier, BIGNUM* secKeyECEG, EC_POINT* pubKeyECEG) {
    float embA[embSize];
    float embB[embSize];
    float lengA = 0, lengB = 0;

    float cDist, cDist_Paillier, cDist_ECEG;
    paillierCompTimer = 0;
    paillierDecTimer = 0;
    ecegFindyTimer = 0;
    ecegCompTimer = 0;
    ecegDecTimer = 0;
    ecegUnmapTimer = 0;
    paillierPrecision = 0;
    ecegPrecision = 0;

    //Cos baby step init and precompute -nG
    constexpr int bsgs_n = 1 << bitPrecision;
    std::unordered_map<unslong, uint> baby_steps;
    computeBabySteps<bsgs_n>(baby_steps);
    BIGNUM* n_bn = BN_new();
    BN_set_word(n_bn, bsgs_n);
    EC_POINT* neg_nG = EC_POINT_new(curve);
    EC_POINT_mul(curve, neg_nG, n_bn, NULL, NULL, ctx);
    EC_POINT_invert(curve, neg_nG, ctx);

    std::cout << "EMB LENGTH: " << embSize << "    BIT PRECISION: " << bitPrecision << "    ITERATIONS: " << totalIterations << "\n";
    //main loop
    for (int n = 0; n < totalIterations; n++) {
        //Get embeddings and normalize them
        lengA = 0;
        lengB = 0;
        for (int i = 0; i < embSize; i++) {
            embA[i] = unif(gen);
            embB[i] = unif(gen);
            lengA += embA[i]*embA[i];
            lengB += embB[i]*embB[i];
        }
        lengA = sqrt(lengA);
        lengB = sqrt(lengB);
        for (int i = 0; i < embSize; i++) {
            embA[i] /= lengA;
            embB[i] /= lengB;
        }

        cDist = cosDist_normalized<embSize>(embA, embB);
        //eDist = eucDist_normalized<embSize>(embA, embB);

        cDist_Paillier = pallierCosDis_normalized<embSize, bitPrecision>(embA, embB, keyPaillier);
        paillierPrecision += std::abs((cDist - cDist_Paillier) / cDist) * 100;

        cDist_ECEG = ecegCosDis_normalized<embSize, bitPrecision, bsgs_n>(embA, embB, secKeyECEG, pubKeyECEG, neg_nG, baby_steps);
        ecegPrecision += std::abs((cDist - cDist_ECEG) / cDist) * 100;
    }
    paillierCompTimer /= totalIterations;
    paillierDecTimer /= totalIterations;
    paillierPrecision = 100 - paillierPrecision/totalIterations;
    ecegFindyTimer /= totalIterations;
    ecegCompTimer /= totalIterations;
    ecegDecTimer /= totalIterations;
    ecegUnmapTimer /= totalIterations;
    ecegPrecision = 100 - ecegPrecision/totalIterations;

    std::cout << "Paillier Dist Avrg Computation Time: " << static_cast<double>(paillierCompTimer/1000)/1000 << "ms\n";
    std::cout << "Paillier Dist Avrg Decryption Time:  " << static_cast<double>(paillierDecTimer/1000)/1000 << "ms\n";
    std::cout << "Paillier Fixed Decimal Accuracy: " << paillierPrecision << "%\n";  //ecegPrecision should is identical to paillier precision)
    std::cout << "EC ElGamal Dist Avrg Decompression Time:  " << static_cast<double>(ecegFindyTimer/1000)/1000 << "ms\n";
    std::cout << "EC ElGamal Dist Avrg Computation Time: " << static_cast<double>(ecegCompTimer/1000)/1000 << "ms\n";
    std::cout << "EC ElGamal Dist Avrg Decryption Time:  " << static_cast<double>(ecegDecTimer/1000)/1000 << "ms\n";
    std::cout << "EC ElGamal Dist Avrg Unmapping Time:  " << static_cast<double>(ecegUnmapTimer/1000)/1000 << "ms\n";
    std::cout << "EC ElGamal Fixed Decimal Accuracy: " << ecegPrecision << "%\n";
    std::cout << "------------------------------------------------------------\n" << std::endl;

    //BSGS cleanup
    BN_free(n_bn);
    EC_POINT_free(neg_nG);
}

int main()
{
    std::cout << "PROGRAM START\n";
    std::cout << "OpenSSL VERSION: " << OpenSSL_version(OPENSSL_VERSION) << "\n";
    std::cout << "PAILLIER MOD LENGTH: 2048 bits, SQUARED MOD LENGTH: 4096 bits\n";
    std::cout << "ECEG UTILIZED CURVE: secp256r1\n";
    std::cout << "------------------------------------------------------------\n\n";

    //auxiliary stuff
    aux_x_bn = BN_new();
    std::random_device rd;
    resultSeed = rd();
    std::mt19937 gen(resultSeed);
    std::uniform_real_distribution<float> unif(0.001, 0.999);   //0 and 1 need to be accounted for and slightly increased of decreased

    //Paillier key
    ipcl::KeyPair keyPaillier = ipcl::generateKeypair(2048, true);

    //ECEG init
    OPENSSL_init_crypto(OPENSSL_INIT_LOAD_CONFIG, NULL);
    int nid = NID_X9_62_prime256v1;//secp256r1
    curve = EC_GROUP_new_by_curve_name(nid);
    ctx = BN_CTX_new();
    curve_p = BN_new();
    curve_a = BN_new();
    curve_b = BN_new();
    EC_GROUP_get_curve(curve, curve_p, curve_a, curve_b, ctx);

    //ECEG order - 1 for random generation
    orderMinus1 = BN_new();
    EC_GROUP_get_order(curve, orderMinus1, ctx);
    BN_sub_word(orderMinus1, 1);

    //ECEG sk and pk
    BIGNUM* secKeyECEG = BN_new();
    EC_POINT* pubKeyECEG = EC_POINT_new(curve);
    BN_rand_range(secKeyECEG, orderMinus1); //rand number in range [0, order-2]
    BN_add_word(secKeyECEG, 1);             //rand number in range [1, order-1]
    EC_POINT_mul(curve, pubKeyECEG, secKeyECEG, NULL, NULL, ctx);

    //OPTIONAL: run hash checking function
    testBSGSHashFunction<16>();

    //the iteration number is decreased to avoid the program running forever (since embedding size increases exponentially)
    timeTest<EMB_LENG_0, BIT_PRECISION_0, 80>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_0, BIT_PRECISION_1, 80>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_0, BIT_PRECISION_2, 80>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);

    timeTest<EMB_LENG_1, BIT_PRECISION_0, 40>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_1, BIT_PRECISION_1, 40>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_1, BIT_PRECISION_2, 40>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);

    timeTest<EMB_LENG_2, BIT_PRECISION_0, 20>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_2, BIT_PRECISION_1, 20>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_2, BIT_PRECISION_2, 20>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);

    timeTest<EMB_LENG_3, BIT_PRECISION_0, 10>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_3, BIT_PRECISION_1, 10>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);
    timeTest<EMB_LENG_3, BIT_PRECISION_2, 10>(gen, unif, keyPaillier, secKeyECEG, pubKeyECEG);

    //cleanup
    BN_free(aux_x_bn);
    BN_free(curve_p);
    BN_free(curve_a);
    BN_free(curve_b);
    BN_free(orderMinus1);
    BN_free(secKeyECEG);
    EC_POINT_free(pubKeyECEG);
    EC_GROUP_free(curve);
    BN_CTX_free(ctx);


    std::cout << "PROGRAM END\n";
    char stopChar = std::getchar();
    std::cout << stopChar << std::endl;

    return 0;
}

