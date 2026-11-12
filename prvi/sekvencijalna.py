import numpy as np

from algoritam import (POCETNA_SEKVENCA, Kandidat, azuriraj, generisi,
                       oceni, pocetna_populacija, proveri_parametre,
                       zajednicki_argparser)
from mere import Statistika, hipervolumen_arhive, sacuvaj_rezultate

G = 100


def trazi(B, N, k, M, lam, seme, start=POCETNA_SEKVENCA, G=G,
          generisi_fn=None):
    rng = np.random.default_rng(seme)
    proveri_parametre(B, N, k, M, lam)
    stat = Statistika()
    evidencija = {}
    budzet = B
    naredni = generisi_fn or (lambda: generisi(populacija, lam, rng))

    sekvence = pocetna_populacija(start, N, k, rng)
    novi = []
    for seq in sekvence:
        if budzet <= 0:
            break
        budzet -= 1
        stat.evaluacije += 1
        ocene = oceni(seq)
        evidencija[seq] = ocene
        novi.append(Kandidat(seq, tuple(ocene)))
    populacija, arhiva = azuriraj([], [], novi, N)

    neuspeli = 0
    while budzet > 0 and neuspeli < G:
        stat.iteracija_pocela()
        novi = []
        while len(novi) < M and budzet > 0 and neuspeli < G:
            seq = naredni()
            if seq in evidencija:
                neuspeli += 1
                continue
            budzet -= 1
            stat.evaluacije += 1
            ocene = oceni(seq)
            evidencija[seq] = ocene
            novi.append(Kandidat(seq, tuple(ocene)))
            neuspeli = 0
        populacija, arhiva = azuriraj(populacija, arhiva, novi, N)
        stat.iteracija_zavrsena()

    return arhiva, evidencija, stat


if __name__ == '__main__':
    args = zajednicki_argparser('P1 — sekvencijalna pretraga (jedna nit)').parse_args()
    arhiva, _evid, stat = trazi(args.B, args.N, args.k, args.M, args.lam,
                                args.seme, args.start)
    hv = hipervolumen_arhive(arhiva)
    parametri = {'B': args.B, 'N': args.N, 'k': args.k, 'M': args.M,
                 'lambda': args.lam, 'seme': args.seme, 'T': 1, 'P': 1, 'G': G}
    tekst = f'evaluacija: {stat.evaluacije}, arhiva: {len(arhiva)}, ' \
            f'hipervolumen: {hv:.5f}, trajanje: {stat.ukupno_trajanje():.2f}s'
    print(tekst)
    sacuvaj_rezultate(args.izlaz, args.json,
                      f'P1 — sekvencijalna pretraga (seme {args.seme})',
                      parametri, stat.to_dict(), arhiva, hv)
    if args.izlaz:
        print('izveštaj:', args.izlaz)
