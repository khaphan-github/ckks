// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "seal/seal.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;
using namespace seal;

void ckks_demo()
{
    // 1. Khởi tạo tham số CKKS
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));

    // 2. Thiết lập scale
    double scale = pow(2.0, 40);

    // 3. Tạo context và các loại khóa
    SEALContext context(parms);
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);

    // 4. Chuẩn bị dữ liệu đầu vào
    size_t slot_count = encoder.slot_count();
    vector<double> input(slot_count);
    double step = 1.0 / (slot_count - 1);
    for (size_t i = 0; i < slot_count; ++i) input[i] = i * step;

    // 5. Encode và mã hóa input
    Plaintext x_plain;
    encoder.encode(input, scale, x_plain);
    Ciphertext x1_encrypted;
    encryptor.encrypt(x_plain, x1_encrypted);

    // 6. Encode các hệ số
    Plaintext plain_coeff3, plain_coeff1, plain_coeff0;
    encoder.encode(3.14159265, scale, plain_coeff3);
    encoder.encode(0.4, scale, plain_coeff1);
    encoder.encode(1.0, scale, plain_coeff0);

    // 7. Tính x^2, relinearize, rescale
    Ciphertext x2_encrypted;
    evaluator.square(x1_encrypted, x2_encrypted);
    evaluator.relinearize_inplace(x2_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x2_encrypted);

    // 8. Tính PI*x, rescale
    Ciphertext pix_encrypted;
    evaluator.multiply_plain(x1_encrypted, plain_coeff3, pix_encrypted);
    evaluator.rescale_to_next_inplace(pix_encrypted);

    // 9. Tính PI*x^3 = (PI*x)*x^2, relinearize, rescale
    evaluator.multiply_inplace(x2_encrypted, pix_encrypted);
    evaluator.relinearize_inplace(x2_encrypted, relin_keys);
    evaluator.rescale_to_next_inplace(x2_encrypted);

    // 10. Tính 0.4*x, rescale
    Ciphertext x1_scaled = x1_encrypted;
    evaluator.multiply_plain_inplace(x1_scaled, plain_coeff1);
    evaluator.rescale_to_next_inplace(x1_scaled);

    // 11. Chuẩn hóa scale về 2^40
    x2_encrypted.scale() = scale;
    x1_scaled.scale() = scale;

    // 12. Modulus switch về cùng level
    parms_id_type last_parms_id = x2_encrypted.parms_id();
    evaluator.mod_switch_to_inplace(x1_scaled, last_parms_id);
    evaluator.mod_switch_to_inplace(plain_coeff0, last_parms_id);

    // 13. Cộng các thành phần
    Ciphertext encrypted_result;
    evaluator.add(x2_encrypted, x1_scaled, encrypted_result);
    evaluator.add_plain_inplace(encrypted_result, plain_coeff0);

    // 14. Giải mã và decode kết quả
    Plaintext plain_result;
    decryptor.decrypt(encrypted_result, plain_result);
    vector<double> result;
    encoder.decode(plain_result, result);

    // In ra 5 giá trị đầu tiên để kiểm tra
    for (size_t i = 0; i < 5; ++i)
        cout << result[i] << endl;
}
