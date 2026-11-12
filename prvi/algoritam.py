import os

_TU = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(os.getcwd(), 'blosum62.json')):
    os.chdir(_TU)

import argparse
from dataclasses import dataclass

import numpy as np

from oceni import oceni  # noqa: F401


ALFABET = 'ACDEFGHIKLMNPQRSTVWY'
ALFABET_BEZ = {s: ''.join(a for a in ALFABET if a != s) for s in ALFABET}
POCETNA_SEKVENCA = "TEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTN"


@dataclass
class Kandidat:
    seq: str
    ocene: tuple
    rang: int = 0
    razr: float = 0.0


def dominira(a, b):
    return (all(x >= y for x, y in zip(a.ocene, b.ocene))
            and any(x > y for x, y in zip(a.ocene, b.ocene)))


def razrijedenost(sloj):
    if len(sloj) <= 2:
        for k in sloj:
            k.razr = float('inf')
        return
    doprinos = {id(k): 0.0 for k in sloj}
    for i in range(3):
        sortirani = sorted(sloj, key=lambda k: k.ocene[i])
        lo = sortirani[0].ocene[i]
        hi = sortirani[-1].ocene[i]
        if hi == lo:
            continue
        doprinos[id(sortirani[0])] = float('inf')
        doprinos[id(sortirani[-1])] = float('inf')
        for j in range(1, len(sortirani) - 1):
            doprinos[id(sortirani[j])] += (
                (sortirani[j + 1].ocene[i] - sortirani[j - 1].ocene[i]) / (hi - lo)
            )
    for k in sloj:
        k.razr = doprinos[id(k)]


def _dominiran(k, svi):
    return any(dominira(d, k) for d in svi)


def pareto_slojevi(kandidati, pool=None):
    preostali = list(kandidati)
    slojevi = []
    rang = 1
    while preostali:
        if pool is not None:
            poslovi = [pool.submit(_dominiran, k, preostali) for k in preostali]
            front = [k for k, p in zip(preostali, poslovi) if not p.result()]
        else:
            front = [k for k in preostali
                      if not any(dominira(d, k) for d in preostali)]
        if not front:
            front = list(preostali)
        for k in front:
            k.rang = rang
        razrijedenost(front)
        slojevi.append(front)
        u_frontu = {id(k) for k in front}
        preostali = [k for k in preostali if id(k) not in u_frontu]
        rang += 1
    return slojevi


def azuriraj(populacija, arhiva, novi, N, rang_pool=None):
    spoj = {}
    for k in list(populacija) + list(arhiva) + list(novi):
        spoj[k.seq] = k
    if not spoj:
        return [], []
    pareto_slojevi(list(spoj.values()), pool=rang_pool)
    nova_arhiva = [k for k in spoj.values() if k.rang == 1]
    sortirani = sorted(spoj.values(), key=lambda k: (k.rang, -k.razr))
    return sortirani[:N], nova_arhiva


def izaberi_roditelje(populacija, rng):
    n = len(populacija)
    if n < 2:
        return populacija[0], populacija[0]
    i = int(rng.integers(0, n))
    j = i
    while j == i:
        j = int(rng.integers(0, n))
    return populacija[i], populacija[j]


def ukrsti(p1, p2, rng):
    L = len(p1)
    presek = int(rng.integers(1, L))
    return p1[:presek] + p2[presek:]


def mutiraj(seq, lam, rng):
    L = len(seq)
    s = int(rng.poisson(lam))
    if s == 0:
        return seq
    s = min(s, L)
    pozicije = rng.choice(L, size=s, replace=False)
    slova = list(seq)
    for p in pozicije:
        p = int(p)
        slova[p] = str(rng.choice(list(ALFABET_BEZ[seq[p]])))
    return ''.join(slova)


def generisi(populacija, lam, rng):
    p1, p2 = izaberi_roditelje(populacija, rng)
    return mutiraj(ukrsti(p1.seq, p2.seq, rng), lam, rng)


def pocetni_mutant(start, k, rng):
    L = len(start)
    br = max(1, int(rng.integers(1, max(2, k + 1))))
    br = min(br, L)
    pozicije = rng.choice(L, size=br, replace=False)
    slova = list(start)
    for p in pozicije:
        p = int(p)
        slova[p] = str(rng.choice(list(ALFABET_BEZ[start[p]])))
    return ''.join(slova)


def pocetna_populacija(start, N, k, rng):
    pop = [start]
    vidjeno = {start}
    while len(pop) < N:
        m = pocetni_mutant(start, k, rng)
        if m not in vidjeno:
            vidjeno.add(m)
            pop.append(m)
    return pop


def zajednicki_argparser(opis):
    p = argparse.ArgumentParser(opis)
    p.add_argument('--B', type=int, default=20000, help='ukupan budžet evaluacija')
    p.add_argument('--N', type=int, default=60, help='najviše kandidata u populaciji')
    p.add_argument('--k', type=int, default=5, help='max izmena u početnom mutantu')
    p.add_argument('--M', type=int, default=30, help='novih kandidata po iteraciji')
    p.add_argument('--lam', type=float, default=2.0, help='Poisson parametar mutacije')
    p.add_argument('--seme', type=int, default=42, help='seme slučajnih brojeva')
    p.add_argument('--start', type=str, default=POCETNA_SEKVENCA, help='početna sekvenca')
    p.add_argument('--izlaz', type=str, default='', help='putanja markdown izveštaja')
    p.add_argument('--json', type=str, default='', help='putanja JSON statistike (za izmeri.py)')
    return p


def proveri_parametre(B, N, k, M, lam):
    assert B >= 1 and N >= 2 and k >= 1 and M >= 1 and lam > 0, \
        'očekivano: B>=1, N>=2, k>=1, M>=1, lam>0'
