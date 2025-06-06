# `MIACL` PoC Implementation, based on the implementation of https://github.com/k4m4/bm-poc/tree/main

> Proof-of-concept implementations of anonymous credentials with decentralized issuance.

# Relevant
The implementation of MIACL can be found in bm_poc/py/miacl.py, while the implementation of the blind multi snowblind can be found in bm_poc/py/bm_sb.py

To rerun the evaluation:

```sh
cd py
poetry shell
poetry install
maturin develop --release
python miacl.py
python bm_sb.py
```
