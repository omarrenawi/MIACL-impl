# This code is built on top of https://github.com/k4m4/bm-poc

from hashlib import sha256
from unittest import TestCase
from utils import serialize, controller, multi_controller
from ec import ECPoint, FQ
import time

class MIACLException(Exception):
    pass

class PoK_DL:
    def __init__(self, ec_point: ECPoint, fq: FQ, hash):
        self.g = ec_point.G1()
        self.hash = hash
        self.fq = fq
        self.ec_point = ec_point
        self.domain = self.hash(b"DOMAIN_PoK_DL")

    
    def _H(self, *args) -> int:
        return int.from_bytes(
            self.hash(serialize(args)), 
            'big'
        )
    
    def _H_FQ(self, *args) -> FQ: #TODO
        return self.fq(self._H(self.domain, *args))
    

    def H_p(self, *args) -> FQ:
        #return self._H_FQ(b"com", *args)
        return self._H_FQ(*args)

    def rand(self) -> FQ:
        return self.fq.rand()
    
    def prove(self, h,  op, basis):
        size= len(op)
        if size != len(basis):
            raise MIACLException("ABORT")
        nonce = [self.rand() for _ in range(size)]
        R = self.ec_point.sum([basis[i] * nonce[i] for i in range(size)])
        c = self.H_p(h,R)
        S = [ nonce[i] - c * op[i] for i in range(size) ]
        return (R, S)

    def verify(self, h, basis, pi):
        (R, S) = pi
        size= len(S)
        if size != len(basis):
            raise MIACLException("ABORT")
        c= self.H_p(h, R)
        return h * c +  self.ec_point.sum([basis[i]* S[i] for i in range(size)]) == R

class MIACL:
    def __init__(self, ec_point: ECPoint, fq: FQ, max_int_size, hash, num_of_signers, num_of_attrs):
        self.ec_point = ec_point
        self.fq = fq
        self.max_int_size = max_int_size
        self.hash = hash
        self.domain = self.hash(b"DOMAIN_MIACL")
        self.num_of_signers = num_of_signers
        self.num_of_attrs=num_of_attrs
        self.g = ec_point.G1()
        self.h = self.g * self.fq.rand() 
        self.t = self.g * self.fq.rand()
        self.h_i = [self.g * self.rand() for _ in range(self.num_of_attrs+1)]
        self.pok_dl = PoK_DL(self.ec_point, self.fq, self.hash)
        self.registration=list()


    def _H(self, *args) -> int: #TODO
        return int.from_bytes(
            self.hash(serialize(args)),
            'big'
        )

    def _H_FQ(self, *args) -> FQ: #TODO
        return self.fq(self._H(self.domain, *args))

    def H_sig(self, *args) -> FQ:
        return self._H_FQ(*args)

    def H_com(self, *args) -> FQ:
        return self._H_FQ(*args)

    def H_rnd(self, *args) -> FQ:
        return self._H_FQ(*args)
    
    def H_sh(self, *args) -> FQ:
        return self._H_FQ(*args)

    def rand(self) -> FQ:
        return self.fq.rand()


    def keygen(self):
        sks = []
        pks = []
        for _ in range(self.num_of_signers):
            sk = self.fq.rand()
            pk = self.g * sk
            sks.append(sk)
            pks.append(pk)

        return (sks, pks)
    
    def verify(self, pks, m, sigma):
        (mu, zeta, zeta_1, b_bar, pi, R_bar, y_bar, z_bar) = sigma

        c_i_bar = [self.H_sig(pks, pk, zeta, zeta_1, R_bar, self.h * y_bar + (zeta+self.ec_point.neg(zeta_1)) * b_bar , self.t * mu + zeta * y_bar ,m) for pk in pks]
        zeta1_check = self.pok_dl.verify(zeta_1, self.h_i + [self.g], pi)

        LHS = R_bar + self.ec_point.sum([pks[i] * (c_i_bar[i] + y_bar**3) for i in range(self.num_of_signers)])

        RHS = self.g * z_bar + self.h * y_bar + (zeta + self.ec_point.neg(zeta_1)) * b_bar
        R_bar_check = LHS == RHS

        return y_bar != 0 and R_bar_check and zeta1_check

    def U_reg(self, L, S):
        L_0 = self.rand()
        L_i = [L_0] + L
        C = self.ec_point.sum(
            [self.h_i[i] * L_i[i] for i in range(self.num_of_attrs+1)]
            )
        pi= self.pok_dl.prove(C, L_i, self.h_i)
        for s in S:
            s.send((C, pi))
        yield (C,L_0)

    def S_reg(self):
        C, pi = yield
        if self.pok_dl.verify(C, self.h_i, pi):
            self.registration.append(C)
            yield

    def reg(self, L):
        return multi_controller(
            lambda *args: self.U_reg(L, *args),
            [lambda *args: self.S_reg(*args) for _ in range(self.num_of_signers)]
        )

    def U_sign(self, pks, m, C, L_i, S):
        rnd_i= [0] * self.num_of_signers
        for i,s in enumerate(S):
            rnd_i[i]=s.send(())

        for s in S:
            s.send(rnd_i.copy())


        rnd = self.H_rnd(rnd_i)


        t_1 = self.g * rnd + C
        t_2 = self.t + self.ec_point.neg(t_1)

        A_i = []
        B_i = []
        com_i = []
        beta_i = []
        for s in S:
            (A, B, com) = s.send(())
            A_i.append(A)
            B_i.append(B)
            com_i.append(com)
            beta_i.append(self.fq.rand())

        alpha = self.rand()
        r = self.fq.rand()
        gamma = self.rand()
        gamma_cubed = gamma**3
        tau = self.rand()
        eta= self.rand()
        zeta = self.t * gamma
        zeta_1 = t_1 * gamma
        zeta_2= zeta + self.ec_point.neg(zeta_1)
        theta = self.t * tau
        A_sum = self.ec_point.sum(A_i)
        B_sum = self.ec_point.sum(B_i)
        alpha_cubed = alpha**3
        B_bar= self.ec_point.sum([B_sum * (gamma * alpha),zeta_2 * eta])
        R_bar = self.ec_point.sum(
            [
                self.g * r,
                self.ec_point.sum([pks[i] * (alpha_cubed * gamma_cubed * beta_i[i]) for i in range(self.num_of_signers)]),
                A_sum * (alpha_cubed * gamma_cubed),
                B_bar
            ]
        )

        alpha_neg_cubed = alpha**(-3)
        gamma_neg_cubed = gamma**(-3)
        c_j = []
        b_i = []
        y_i = []
        for i, s in enumerate(S):
            c_j.append(self.H_sig(pks, pks[i], zeta, zeta_1, R_bar, B_bar, theta, m) * alpha_neg_cubed * gamma_neg_cubed + beta_i[i])
            com_ic = com_i.copy()
            com_ic[i] = b''
            B_ic=B_i.copy()
            b, y = s.send((c_j[i], com_ic, B_ic))
            b_i.append(b)
            y_i.append(y)

        z_i = []
        for i, s in enumerate(S):
            y_ic = y_i.copy()
            y_ic[i] = b''
            b_ic=b_i.copy()

            z = s.send((b_ic, y_ic))
            z_i.append(z)

        b_sum = self.fq.sum(b_i)
        y_sum = self.fq.sum(y_i)
        z_sum = self.fq.sum(z_i)

        B_check = B_sum == t_2 * b_sum + self.h * y_sum
        A_check = self.g * z_sum == A_sum + self.ec_point.sum([pks[i] * (c_j[i] + y_sum**3) for i in range(self.num_of_signers)])
        if not A_check or not B_check:
            raise MIACLException("ABORT")

        z_bar = r + gamma_cubed * alpha_cubed * z_sum
        y_bar = alpha * y_sum * gamma
        b_bar = alpha * b_sum + eta
        mu = tau - gamma * y_bar
        op = [gamma * el for el in L_i + [rnd]]
        pi = self.pok_dl.prove(zeta_1, op ,self.h_i + [self.g])
        sigma = (mu, zeta, zeta_1, b_bar, pi, R_bar, y_bar, z_bar)

        if not self.verify(pks, m, sigma):
            raise MIACLException("ABORT")

        show_par = (rnd, gamma) # we get L_0 from the ourput of U_reg
        yield (sigma, show_par)

    def S_sign(self, i, pk, sk, C):
        yield
        rnd = self.rand()
        rnd_i=yield(rnd)
        rnd_i[i]=rnd
        rnd = self.H_rnd(rnd_i)

        if C not in self.registration:
            raise MIACLException("ABORT")
        a = self.fq.rand()
        b = self.fq.rand()
        y = self.rand()
        t_1 = self.g * rnd + C
        t_2 = self.ec_point.sum([self.t, self.ec_point.neg(t_1)])
        A = self.g * a
        B = t_2 * b + self.h * y
        com = self.H_com(i, b, y)
        yield
        c, com_i, B_i = yield (A, B, com)
        b_i, y_i = yield ((b, y))

        com_i[i] = com
        y_i[i] = y
        b_i[i] = b

        if any(com_i[i] != self.H_com(i, b_i[i], y_i[i]) for i in range(self.num_of_signers)):
            raise MIACLException("ABORT")
                

        if any(B_i[i] != self.h * y_i[i] + t_2 * b_i[i]   for i in range(self.num_of_signers)):
            raise MIACLException("ABORT")

        y_sum = self.fq.sum(y_i)
        z = (a + (c + y_sum**3) * sk)
        yield z

    def sign(self, pks, sks, m, C, L_i):
        return multi_controller(
            lambda *args: self.U_sign(pks, m, C, L_i, *args),
            [(lambda i, pk, sk, C: lambda *args: self.S_sign(i, pk, sk, C, *args))(i, pks[i], sks[i], C) for i in range(self.num_of_signers)]
        )
    

    def U_show(self, apk, m, sig, L, showpar):
        L_0, rnd, gamma = showpar
        L= [L_0]+ L
        Gamma = self.g * gamma
        Psi_i = [self.h_i[i] * gamma for i in range(self.num_of_attrs+1)]

        r_sdl = self.rand()
        h_sdl_i = [self.h_i[i] * r_sdl for i in range(self.num_of_attrs+1)]
        g_sdl = self.g * r_sdl
        t_sdl = self.t * r_sdl

        r_i = [self.rand() for _ in range(self.num_of_attrs+1)]
        r_g = self.rand()

        R = Gamma * r_g + self.ec_point.sum(
            [Psi_i[i] * r_i[i] for i in range(self.num_of_attrs+1)]
        )

        c = self.H_sh(apk, m, sig, L, L_0, Gamma, Psi_i, h_sdl_i, g_sdl, t_sdl, R)

        s_sdl = r_sdl - c * gamma
        s_i = [r_i[i] - c * L[i] for i in range(self.num_of_attrs+1)]
        s_Gamma = r_g - rnd * c

        pi_op = (R, s_i + [s_Gamma])
        pi_sdl = (Gamma, Psi_i, g_sdl, t_sdl, h_sdl_i, s_sdl)
        return (pi_op, pi_sdl)



    def V_show(self, apk, m, sig, L, L_0, pi):
        (pi_op, pi_sdl) = pi
        (R, S) = pi_op
        (s_i, s_Gamma) = S[:-1], S[-1]
        (Gamma, Psi_i, g_sdl, t_sdl, h_sdl_i, s_sdl) = pi_sdl
        L= [L_0]+L

        _, zeta, zeta_1, *_ = sig

        if not self.verify(apk, m, sig):
            return 0
        c = self.H_sh(apk, m, sig, L, L_0, Gamma, Psi_i, h_sdl_i, g_sdl, t_sdl, R)

        if self.g * s_sdl + Gamma * c != g_sdl:
            return 0

        if self.t * s_sdl + zeta * c != t_sdl:
            return 0
        
        for i in range(self.num_of_attrs+1):
            if Psi_i[i] * c + self.h_i[i] * s_sdl != h_sdl_i[i]:
                return 0
            
        lhs = zeta_1 * c + Gamma * s_Gamma + self.ec_point.sum(
            [Psi_i[i] * s_i[i] for i in range(self.num_of_attrs+1)]
        )

        if lhs != R:
            assert False

        return 1
    
    def show(self, apk, m, sig, L, L_0, show_par):
        pi= self.U_show(apk, m,sig, L, show_par)
        return self.V_show(apk, m, sig, L, L_0, pi)

from ec import from_ecc_py
import py_eth_pairing

BN128FQ, BN128Point = from_ecc_py('BN128', py_eth_pairing)



def sign_eval(trials, num_of_signers, num_of_attrs):
    # token generation
    trial_res=list()
    for _ in range(trials):
        ec_point, fq, max_int_size = (BN128Point, BN128FQ, 32)
        m = BN128FQ.rand()
        hash = lambda x: sha256(x).digest()
        bm = MIACL(ec_point, fq, max_int_size, hash, num_of_signers, num_of_attrs)
        (sks, pks) = bm.keygen()
        L= [bm.rand() for _ in range(num_of_attrs)]
        C, L_0 = bm.reg(L)
        L_i=[L_0]+L
        start = time.time_ns()
        bm.sign(pks, sks, m, C, L_i)
        end = time.time_ns()
        time_ms = (end-start) / 10**6
        trial_res.append(time_ms)
    return round(sum(trial_res)/ trials,2)

def vrfy_eval(trials, num_of_signers, num_of_attrs):
    # token verification 
    trial_res=list()
    for _ in range(trials):
        ec_point, fq, max_int_size = (BN128Point, BN128FQ, 32)
        m = BN128FQ.rand()
        hash = lambda x: sha256(x).digest()
        bm = MIACL(ec_point, fq, max_int_size, hash, num_of_signers, num_of_attrs)
        (sks, pks) = bm.keygen()
        L= [bm.rand() for _ in range(num_of_attrs)]
        C, L_0 = bm.reg(L)
        L_i=[L_0]+L
        (sigma, show_par) = bm.sign(pks, sks, m, C, L_i)
        start = time.time_ns()
        bm.verify(pks, m, sigma)
        end = time.time_ns()
        time_ms = (end-start) / 10**6
        trial_res.append(time_ms)
    return round(sum(trial_res)/ trials,2)

def cred_presnt_eval(trials, num_of_signers, num_of_attrs):
    # Credentials presentation
    trial_res=list()
    for _ in range(trials):
        ec_point, fq, max_int_size = (BN128Point, BN128FQ, 32)
        m = BN128FQ.rand()
        hash = lambda x: sha256(x).digest()
        bm = MIACL(ec_point, fq, max_int_size, hash, num_of_signers, num_of_attrs)
        (sks, pks) = bm.keygen()
        L= [bm.rand() for _ in range(num_of_attrs)]
        C, L_0 = bm.reg(L)
        L_i=[L_0]+L
        (sigma, show_par) = bm.sign(pks, sks, m, C, L_i)
        (rnd, gamma)= show_par
        show_par=(L_0, rnd, gamma)
        start = time.time_ns()
        pi = bm.U_show(pks, m, sigma, L, show_par)
        end = time.time_ns()
        time_ms = (end-start) / 10**6
        trial_res.append(time_ms)
    return round(sum(trial_res)/ trials,2)

def cred_ver_eval(trials, num_of_signers, num_of_attrs):
    # Credentials verification
    trial_res=list()
    for _ in range(trials):
        ec_point, fq, max_int_size = (BN128Point, BN128FQ, 32)
        m = BN128FQ.rand()
        hash = lambda x: sha256(x).digest()
        bm = MIACL(ec_point, fq, max_int_size, hash, num_of_signers, num_of_attrs)
        (sks, pks) = bm.keygen()
        L= [bm.rand() for _ in range(num_of_attrs)]
        C, L_0 = bm.reg(L)
        L_i=[L_0]+L
        (sigma, show_par) = bm.sign(pks, sks, m, C, L_i)
        (rnd, gamma)= show_par
        show_par=(L_0, rnd, gamma)
        pi = bm.U_show(pks, m, sigma, L, show_par)
        start = time.time_ns()
        bm.V_show(pks, m, sigma, L, L_0, pi)
        end = time.time_ns()
        time_ms = (end-start) / 10**6
        trial_res.append(time_ms)
    return round(sum(trial_res)/ trials,2)


class TestMIACL(TestCase):
    def test(self):
        for ec_point, fq, max_int_size in [(BN128Point, BN128FQ, 32)]:
            with self.subTest(msg=f"Testing with {ec_point.__name__} and {fq.__name__}"):
                m = BN128FQ.rand()
                num_of_signers = 3
                num_of_attrs= 3
                hash = lambda x: sha256(x).digest()
                bm = MIACL(ec_point, fq, max_int_size, hash, num_of_signers, num_of_attrs)
                L= [bm.rand() for _ in range(num_of_attrs)]
                C, L_0 = bm.reg(L)
                L_i=[L_0]+L
                (sks, pks) = bm.keygen()
                (sigma, show_par) = bm.sign(pks, sks, m, C, L_i)
                self.assertTrue(bm.verify(pks, m, sigma))
                (rnd, gamma)= show_par
                show_par=(L_0, rnd, gamma)
                pi = bm.U_show(pks, m, sigma, L, show_par)
                self.assertTrue(bm.V_show(pks, m, sigma, L, L_0, pi))


class TestMIACL(TestCase):
    def test(self):
        ec_point, fq, max_int_size = (BN128Point, BN128FQ, 32)
        m = BN128FQ.rand()
        hash = lambda x: sha256(x).digest()
        num_of_signers=2
        bm = MIACL(ec_point, fq, max_int_size, hash, num_of_signers, 3)
        (sks, pks) = bm.keygen()
        L= [bm.rand() for _ in range(3)]
        C, L_0 = bm.reg(L)
        L_i=[L_0]+L
        (sigma, show_par) = bm.sign(pks, sks, m, C, L_i)
        (rnd, gamma)= show_par
        show_par=(L_0, rnd, gamma)
        pi = bm.U_show(pks, m, sigma, L, show_par)
        start = time.time_ns()
        self.assertTrue(bm.V_show(pks, m, sigma, L, L_0, pi))
        end = time.time_ns()
        time_ms = (end-start) / 10**6
        print(time_ms)



if __name__=="__main__":
        ec_point, fq, max_int_size = (BN128Point, BN128FQ, 32)
        m = BN128FQ.rand()
        num_of_signers = [2**i for i in range(7)]
        num_of_attrs = [2 * i for i in range(1,8)]
        hash = lambda x: sha256(x).digest()
        trials=10
        eval_file= open("miacl-eval.csv", "w")
        eval_file.write("#Signers,#Attrs,Issuance,Token Verification,Credential Presentation,Credential Verification\n")
        for i in num_of_signers:
            for j in num_of_attrs:
                sign= sign_eval(trials, i, j)
                vrfy= vrfy_eval(trials, i, j)
                cred_pres= cred_presnt_eval(trials, i, j)
                cred_ver= cred_ver_eval(trials, i, j)
                out = f"{i},{j},{sign},{vrfy},{cred_pres},{cred_ver}\n"
                eval_file.write(out)
        eval_file.close()