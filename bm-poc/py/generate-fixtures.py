import os
import json
from hashlib import sha256
from ec import from_ecc_py
import py_eth_pairing
from bm_bls import BM_BLS
from bm_sb import MIACL
BN128FQ, BN128Point = from_ecc_py('BN128', py_eth_pairing)


def write_bm_sb_fixture(pks, m, sigma):
    (R_bar, y_bar, z_bar) = sigma
    data = {
        "R_bar": R_bar.p,
        "m": m.n,
        "pks": list(map(lambda x: x.p, pks)),
        "y_bar": y_bar.n,
        "z_bar": z_bar.n,
    }

    with open('data/input.json', 'w') as f:
        json.dump(data, f)

def generate_bm_sb_fixture(num_signers):
    m = BN128FQ.rand()
    hash = lambda x: sha256(x).digest()
    max_int_size = 32
    bm = MIACL(BN128Point, BN128FQ, max_int_size, hash, num_signers)
    (sks, pks) = bm.keygen()
    sigma = bm.sign(pks, sks, m)
    assert bm.verify(pks, m, sigma)
    write_bm_sb_fixture(pks, m, sigma)

if __name__ == '__main__':
    num_messages = int(os.environ.get('NUM_MESSAGES', 1))
    num_signers = int(os.environ.get('NUM_SIGNERS', 1))
    generate_bm_bls_aggr_fixture(num_messages, num_signers)
    generate_bm_bls_aggr_hm_fixture(num_messages, num_signers)
    generate_bm_sb_fixture(num_signers)
