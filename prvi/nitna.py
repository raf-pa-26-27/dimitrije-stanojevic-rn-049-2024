import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from algoritam import (POCETNA_SEKVENCA, Kandidat, azuriraj, generisi,
                       oceni, pocetna_populacija, proveri_parametre,
                       zajednicki_argparser)
from mere import Statistika, hipervolumen_arhive, sacuvaj_rezultate

G = 100


def podeli(niz, delova):
    delova = max(1, min(delova, len(niz)))
    vel = (len(niz) + delova - 1) // delova
    return [niz[i:i + vel] for i in range(0, len(niz), vel)]


class NitnaPretraga:

    def __init__(self, B, N, k, M, lam, T, seme, start=POCETNA_SEKVENCA, G=G):
        proveri_parametre(B, N, k, M, lam)
        assert T >= 1
        self.B, self.N, self.k, self.M, self.lam = B, N, k, M, lam
        self.T, self.start, self.G = T, start, G

        self.lock = threading.Condition()
        self.rng_lock = threading.Lock()
        self.rng = np.random.default_rng(seme)

        self.rang_pool = ThreadPoolExecutor(max_workers=T)

        self.evidencija = {}
        self.populacija, self.arhiva = [], []
        self.rezervisano = 0
        self.u_toku = 0
        self.prikupljeno = 0
        self.neuspeli = 0
        self.iteracija = 0
        self.stop = False
        self.stat = Statistika()
        self.novi_iteracije = []
        self.otvorena = threading.Event()

    def primi_poruke(self):
        return []

    def posalji_poruke(self, broj):
        pass

    def trazi(self):
        self._pocetna_populacija()
        niti = [threading.Thread(target=self._radnik) for _ in range(self.T)]
        for t in niti:
            t.start()
        try:
            self._glavna_petlja()
        finally:
            with self.lock:
                self.stop = True
            self.otvorena.set()
            for t in niti:
                t.join()
            self.rang_pool.shutdown()
        return self.arhiva, self.evidencija, self.stat

    def _glavna_petlja(self):
        while True:
            with self.lock:
                primljeni = self.primi_poruke()
                self.stat.iteracija_pocela()
                self.iteracija += 1
                self.prikupljeno = 0
                self.novi_iteracije = []
                self.otvorena.set()
                while not (self.stop or self.prikupljeno >= self.M):
                    self.lock.wait()
                while self.u_toku > 0:
                    self.lock.wait()
                self.otvorena.clear()
                novi = [Kandidat(s, self.evidencija[s]) for s in self.novi_iteracije]
                novi.extend(primljeni)
                kraj = self.stop
            if novi:
                self.populacija, self.arhiva = azuriraj(
                    self.populacija, self.arhiva, novi, self.N, self.rang_pool)
            with self.lock:
                self.stat.iteracija_zavrsena()
                self.posalji_poruke(self.iteracija)
                if kraj:
                    break

    def _pocetna_populacija(self):
        with self.rng_lock:
            sekvence = pocetna_populacija(self.start, self.N, self.k, self.rng)
        posao = sekvence[:self.B]
        with self.lock:
            for s in posao:
                self.evidencija[s] = None
            self.rezervisano = len(posao)
            self.stat.evaluacije = len(posao)
            self.u_toku = len(posao)
        niti = [threading.Thread(target=self._evaluiraj_deo, args=(deo,))
                for deo in podeli(posao, self.T)]
        for t in niti:
            t.start()
        for t in niti:
            t.join()
        with self.lock:
            kandidati = [Kandidat(s, self.evidencija[s]) for s in posao
                         if self.evidencija[s] is not None]
        self.populacija, self.arhiva = azuriraj([], [], kandidati, self.N, self.rang_pool)

    def _evaluiraj_deo(self, deo):
        for s in deo:
            self._evaluiraj(s)

    def _evaluiraj(self, seq):
        ocene = oceni(seq)
        with self.lock:
            self.evidencija[seq] = tuple(ocene)
            self.u_toku -= 1
            self.lock.notify_all()

    def _radnik(self):
        while True:
            self.otvorena.wait()
            with self.lock:
                if self.stop:
                    return
                if self.rezervisano >= self.B or self.neuspeli >= self.G:
                    self.stop = True
                    self.lock.notify_all()
                    return
            seq = self._generisi()
            with self.lock:
                if seq in self.evidencija:
                    self.neuspeli += 1
                    if self.neuspeli >= self.G:
                        self.stop = True
                        self.lock.notify_all()
                    continue
                if self.rezervisano >= self.B or self.prikupljeno >= self.M:
                    if self.rezervisano >= self.B:
                        self.stop = True
                        self.lock.notify_all()
                    continue
                self.evidencija[seq] = None
                self.rezervisano += 1
                self.stat.evaluacije += 1
                self.u_toku += 1
                self.neuspeli = 0
                self.novi_iteracije.append(seq)
                self.prikupljeno += 1
                if self.prikupljeno >= self.M:
                    self.otvorena.clear()
                self.lock.notify_all()
            self._evaluiraj(seq)

    def _generisi(self):
        with self.rng_lock:
            return generisi(self.populacija, self.lam, self.rng)


if __name__ == '__main__':
    pars = zajednicki_argparser('P2 — pretraga sa više niti')
    pars.add_argument('--T', type=int, default=4, help='broj niti')
    args = pars.parse_args()
    pretraga = NitnaPretraga(args.B, args.N, args.k, args.M, args.lam,
                             args.T, args.seme, args.start)
    arhiva, _evid, stat = pretraga.trazi()
    hv = hipervolumen_arhive(arhiva)
    parametri = {'B': args.B, 'N': args.N, 'k': args.k, 'M': args.M,
                 'lambda': args.lam, 'seme': args.seme,
                 'T': args.T, 'P': 1, 'G': G}
    print(f'evaluacija: {stat.evaluacije}, arhiva: {len(arhiva)}, '
          f'hipervolumen: {hv:.5f}, trajanje: {stat.ukupno_trajanje():.2f}s')
    sacuvaj_rezultate(args.izlaz, args.json,
                      f'P2 — nitna pretraga, T={args.T} (seme {args.seme})',
                      parametri, stat.to_dict(), arhiva, hv)
    if args.izlaz:
        print('izveštaj:', args.izlaz)
