cd ./SEAL
cmake -S . -B build -DSEAL_BUILD_EXAMPLES=ON
cmake --build build -j$(nproc)
./build/bin/sealexamples

choose 9
# File
/workspaces/crypto_encrypt/SEAL/native/examples/5_ckks_basics.cpp# ckks
