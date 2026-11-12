import numpy as np
from concurrent.futures import ThreadPoolExecutor

from algoritam import (POCETNA_SEKVENCA, ALFABET, Kandidat, azuriraj,
                       dominira, mutiraj, pareto_slojevi, pocetna_populacija,
                       razrijedenost, ukrsti)
from mere import hv2d, hv3d, normalizuj
from sekvencijalna import trazi as trazi_sekvencijalno
from nitna import NitnaPretraga
from procesna import kandidati_za_slanje, pokreni, susedi

OK = '  OK —'


def test_dominacije():
    a = Kandidat('AAA', (5, 5, 5))
    b = Kandidat('CCC', (4, 5, 5))
    c = Kandidat('GGG', (6, 4, 5))
    assert dominira(a, b) and not dominira(b, a)
    assert not dominira(a, c) and not dominira(c, a)
    assert not dominira(a, a)
    print(OK, 'dominacija')


def test_slojeva():
    kand = [Kandidat(s, o) for s, o in [
        ('A', (9, 1, 5)), ('B', (5, 5, 5)), ('C', (1, 9, 5)),
        ('D', (4, 4, 5)), ('E', (9, 2, 5)), ('F', (2, 9, 5))]]
    slojevi = pareto_slojevi(kand)
    imena1 = {k.seq for k in slojevi[0]}
    imena2 = {k.seq for k in slojevi[1]}
    assert imena1 == {'B', 'E', 'F'}, imena1
    assert imena2 == {'A', 'C', 'D'}, imena2
    assert all(k.rang == 1 for k in slojevi[0])
    assert all(k.rang == 2 for k in slojevi[1])
    print(OK, 'Pareto slojevi (E dominira A, F dominira C, B dominira D)')


def test_razrijedenosti():
    sloj = [Kandidat(s, (v, 7.0, 7.0))
            for s, v in [('k1', 1.0), ('k2', 2.0), ('k3', 5.0), ('k4', 10.0)]]
    razrijedenost(sloj)
    po_seq = {k.seq: k.razr for k in sloj}
    assert po_seq['k1'] == float('inf') and po_seq['k4'] == float('inf')
    assert abs(po_seq['k2'] - (5.0 - 1.0) / 9.0) < 1e-9
    assert abs(po_seq['k3'] - (10.0 - 2.0) / 9.0) < 1e-9
    print(OK, 'mera razređenosti (krajevi ∞, sredina (sledeći−prethodni)/raspon)')


def test_ukrstanja():
    rng = np.random.default_rng(7)
    p1 = ALFABET + ALFABET[:12]
    p2 = ALFABET[::-1] + ALFABET[:12]
    assert len(p1) == len(p2) == 32
    for _ in range(50):
        dete = ukrsti(p1, p2, rng)
        assert len(dete) == 32
        assert any(dete == p1[:i] + p2[i:] for i in range(1, 32))
    print(OK, 'ukrštanje (dužina sačuvana, prefiks+sufiks)')


def test_mutacije():
    rng = np.random.default_rng(11)
    seq = POCETNA_SEKVENCA
    for _ in range(200):
        m = mutiraj(seq, 2.0, rng)
        assert len(m) == 32
        assert all(s in ALFABET for s in m)
        razlike = [i for i, (a, b) in enumerate(zip(seq, m)) if a != b]
        for i in razlike:
            assert m[i] != seq[i]
    assert mutiraj(seq, 0.0, rng) == seq
    print(OK, 'mutacija (Poisson broj pozicija, uvek druga slova)')


def test_pocetne_populacije():
    rng = np.random.default_rng(3)
    pop = pocetna_populacija(POCETNA_SEKVENCA, 40, 5, rng)
    assert len(pop) == 40 and pop[0] == POCETNA_SEKVENCA
    assert len(set(pop)) == 40
    assert all(len(s) == 32 for s in pop)
    print(OK, 'početna populacija (start + N−1 različitih mutanata)')


def test_azuriranja():
    a = Kandidat('A', (9, 1, 5))
    b = Kandidat('B', (5, 5, 5))
    c = Kandidat('C', (1, 9, 5))
    d = Kandidat('D', (4, 4, 5))
    populacija, arhiva = azuriraj([], [], [a, b, c, d], 2)
    assert {k.seq for k in arhiva} == {'A', 'B', 'C'}
    assert len(populacija) == 2
    assert d not in populacija or d.rang > 1
    print(OK, 'ažuriranje (prvi sloj = arhiva, populacija do N po rangu)')


def test_determinizma_sekvencijalne():
    a1, _, s1 = trazi_sekvencijalno(B=500, N=20, k=4, M=10, lam=2.0, seme=42)
    a2, _, s2 = trazi_sekvencijalno(B=500, N=20, k=4, M=10, lam=2.0, seme=42)
    assert [k.seq for k in a1] == [k.seq for k in a2]
    assert s1.evaluacije == s2.evaluacije == 500
    a3, _, _ = trazi_sekvencijalno(B=500, N=20, k=4, M=10, lam=2.0, seme=43)
    assert [k.seq for k in a1] != [k.seq for k in a3]
    print(OK, 'determinizam (isto seme → ista arhiva; drugo seme → druga)')


def test_malog_budzeta():
    arhiva, evid, stat = trazi_sekvencijalno(B=5, N=20, k=4, M=10, lam=2.0, seme=1)
    assert stat.evaluacije == 5 and len(evid) == 5
    assert len(arhiva) >= 1
    print(OK, 'mali budžet (B=5 < N: tačno 5 evaluacija)')


def test_uslova_G():
    arhiva, evid, stat = trazi_sekvencijalno(
        B=500, N=8, k=2, M=5, lam=2.0, seme=1, G=10,
        generisi_fn=lambda: POCETNA_SEKVENCA)
    assert stat.evaluacije == 8
    print(OK, 'uslov završetka G (10 duplikata zaredom → kraj)')


def test_nitne_verzije():
    tr = NitnaPretraga(B=400, N=20, k=4, M=10, lam=2.0, T=4, seme=42)
    arhiva, evid, stat = tr.trazi()
    assert stat.evaluacije == tr.rezervisano == 400
    assert len(evid) == 400
    assert all(oc is not None for oc in evid.values())
    assert len(arhiva) >= 1
    print(OK, 'nitna verzija (T=4: bez duplikata, budžet tačno 400, bez None)')


def test_paralelnog_rangiranja():
    kand = [Kandidat(s, o) for s, o in [
        ('A', (9, 1, 5)), ('B', (5, 5, 5)), ('C', (1, 9, 5)),
        ('D', (4, 4, 5)), ('E', (9, 2, 5)), ('F', (2, 9, 5))]]
    with ThreadPoolExecutor(max_workers=4) as pool:
        slojevi = pareto_slojevi(kand, pool=pool)
    assert {k.seq for k in slojevi[0]} == {'B', 'E', 'F'}
    assert {k.seq for k in slojevi[1]} == {'A', 'C', 'D'}
    print(OK, 'paralelno rangiranje kroz pool (isti slojevi kao sekvencijalno)')


def test_topologija():
    assert susedi(0, 4, 'prsten') == [1]
    assert susedi(3, 4, 'prsten') == [0]
    assert susedi(0, 4, 'zvezda') == [1, 2, 3]
    assert susedi(2, 4, 'zvezda') == [0]
    assert susedi(1, 4, 'svima') == [0, 2, 3]
    rng = np.random.default_rng(5)
    assert susedi(0, 4, 'nasumican', rng) in ([1], [2], [3])
    print(OK, 'topologije (prsten/zvezda/svima/nasumican)')


def test_izbora_za_slanje():
    arhiva = [Kandidat(s, o) for s, o in [
        ('A', (9, 1, 5)), ('B', (5, 5, 5)), ('C', (1, 9, 5))]]
    pareto_slojevi(arhiva)
    poslato = kandidati_za_slanje(arhiva, 2)
    assert len(poslato) == 2
    assert all(isinstance(o, tuple) and len(o) == 3 for _, o in poslato)
    assert kandidati_za_slanje([], 3) == []
    assert kandidati_za_slanje(arhiva, 5) == [
        (k.seq, tuple(k.ocene)) for k in
        sorted(arhiva, key=lambda k: -k.razr)]
    print(OK, 'izbor za slanje (iz arhive po razređenosti, najviše K, sa ocenama)')


def test_procesne_verzije():
    finalna, rezultati = pokreni(B=240, N=16, k=3, M=8, lam=2.0, T=1, P=2,
                                 K=3, E=1, seme=42, lam_nacin='iste')
    assert len(rezultati) == 2
    assert sum(r['stat']['evaluacije'] for r in rezultati) == 240
    assert len(finalna) >= 1
    assert sum(r['stat']['primljeno'] for r in rezultati) > 0
    print(OK, 'procesna verzija (P=2 prsten: podela budžeta, spajanje arhivi, razmena)')


def test_hipervolumena():
    tacke = normalizuj([(0.8, 0.3, 0.5), (0.5, 0.6, 0.4), (0.2, 0.9, 0.7)],
                       lo=(0, 0, 0), hi=(1, 1, 1))
    hv = hv3d(tacke, ref=(-0.1, -0.1, -0.1))
    assert hv > 0.3, hv
    assert hv2d([(1.0, 1.0)], ref=(0.0, 0.0)) == 1.0
    print(OK, 'hipervolumen (hv2d/hv3d/normalizuj)')


if __name__ == '__main__':
    test_dominacije()
    test_slojeva()
    test_razrijedenosti()
    test_ukrstanja()
    test_mutacije()
    test_pocetne_populacije()
    test_azuriranja()
    test_determinizma_sekvencijalne()
    test_malog_budzeta()
    test_uslova_G()
    test_nitne_verzije()
    test_paralelnog_rangiranja()
    test_topologija()
    test_izbora_za_slanje()
    test_procesne_verzije()
    test_hipervolumena()
    print('\nSVE PROVERE PROŠLE ✓')
