import multiprocessing as mp
import queue as queue_mod
import time

from algoritam import (POCETNA_SEKVENCA, Kandidat, azuriraj,
                       pareto_slojevi, proveri_parametre, zajednicki_argparser)
from mere import Statistika, hipervolumen_arhive, sacuvaj_rezultate
from nitna import NitnaPretraga

LAMBDA_MIN, LAMBDA_MAX = 0.5, 4.0

MAX_PORUKA_U_REDU = 500
ROK_CEKANJA_REZULTATA = 1.0
ROK_ZIDA = 3600


def susedi(rank, P, topologija, rng=None):
    if topologija == 'prsten':
        return [(rank + 1) % P]
    if topologija == 'zvezda':
        return list(range(1, P)) if rank == 0 else [0]
    if topologija == 'svima':
        return [r for r in range(P) if r != rank]
    if topologija == 'nasumican':
        drugi = [r for r in range(P) if r != rank]
        return [drugi[int(rng.integers(0, len(drugi)))]]
    raise ValueError(f'nepoznata topologija: {topologija}')


def kandidati_za_slanje(arhiva, K):
    if not arhiva:
        return []
    izbor = sorted(arhiva, key=lambda k: -k.razr)[:max(0, K)]
    return [(k.seq, tuple(k.ocene)) for k in izbor]


class ProcesnaPretraga(NitnaPretraga):

    def __init__(self, rank, P, K, E, inbox, redovi, topologija='prsten', **kw):
        super().__init__(**kw)
        self.rank, self.P, self.K, self.E = rank, P, K, E
        self.inbox = inbox
        self.redovi = redovi
        self.topologija = topologija
        self.od_drugih = set()

    def primi_poruke(self):
        novi = []
        while True:
            try:
                poruka = self.inbox.get_nowait()
            except queue_mod.Empty:
                break
            for seq, ocene in poruka['kandidati']:
                self.stat.primljeno += 1
                if seq in self.evidencija:
                    self.stat.primljeno_odbaceno += 1
                    continue
                self.evidencija[seq] = tuple(ocene)
                self.od_drugih.add(seq)
                novi.append(Kandidat(seq, tuple(ocene)))
        return novi

    def posalji_poruke(self, broj):
        if broj % self.E != 0:
            return
        posalji = kandidati_za_slanje(self.arhiva, self.K)
        if not posalji:
            return
        with self.rng_lock:
            meta = susedi(self.rank, self.P, self.topologija, self.rng)
        poruka = {'od': self.rank, 'kandidati': posalji}
        for j in meta:
            try:
                self.redovi[j].put(poruka, timeout=0.5)
                self.stat.poslato += len(posalji)
            except queue_mod.Full:
                pass


def _proces(rank, P, cfg, redovi, rezultati):
    inbox = redovi[rank]
    seme_i = cfg['seme'] + 1009 * rank
    if cfg['lam_nacin'] == 'iste':
        lam_i = cfg['lam']
    else:
        korak = (LAMBDA_MAX - LAMBDA_MIN) / max(1, P - 1)
        lam_i = LAMBDA_MIN + rank * korak
    tr = ProcesnaPretraga(rank=rank, P=P, K=cfg['K'], E=cfg['E'],
                          inbox=inbox, redovi=redovi, topologija=cfg['topologija'],
                          B=cfg['B'] // P, N=cfg['N'], k=cfg['k'], M=cfg['M'],
                          lam=lam_i, T=cfg['T'], seme=seme_i, start=cfg['start'])
    arhiva, _evid, stat = tr.trazi()
    od_drugih_u_arhivi = sum(1 for kk in arhiva if kk.seq in tr.od_drugih)
    rezultati.put({
        'rank': rank,
        'lam': lam_i,
        'arhiva': [(kk.seq, tuple(kk.ocene)) for kk in arhiva],
        'stat': stat.to_dict(),
        'od_drugih_u_arhivi': od_drugih_u_arhivi,
    })


def pokreni(B, N, k, M, lam, T, P, K, E, seme, start=POCETNA_SEKVENCA,
            lam_nacin='razlicite', topologija='prsten'):
    proveri_parametre(B, N, k, M, lam)
    cfg = dict(B=B, N=N, k=k, M=M, lam=lam, T=T, K=K, E=E, seme=seme,
               start=start, lam_nacin=lam_nacin, topologija=topologija)
    redovi = [mp.Queue(maxsize=MAX_PORUKA_U_REDU) for _ in range(P)]
    rezultati = mp.Queue()
    procesi = [mp.Process(target=_proces, args=(i, P, cfg, redovi, rezultati))
               for i in range(P)]
    for p in procesi:
        p.start()

    prikupljeno = []
    zid = time.time() + ROK_ZIDA
    while len(prikupljeno) < P and time.time() < zid:
        try:
            prikupljeno.append(rezultati.get(timeout=ROK_CEKANJA_REZULTATA))
        except queue_mod.Empty:
            if not any(p.is_alive() for p in procesi):
                try:
                    prikupljeno.append(rezultati.get(timeout=2.0))
                except queue_mod.Empty:
                    break
    for p in procesi:
        p.join(timeout=5)
        if p.is_alive():
            p.terminate()

    spoj = {}
    for r in prikupljeno:
        for seq, ocene in r['arhiva']:
            spoj[seq] = Kandidat(seq, tuple(ocene))
    if spoj:
        finalna = pareto_slojevi(list(spoj.values()))[0]
    else:
        finalna = []
    return finalna, prikupljeno


if __name__ == '__main__':
    pars = zajednicki_argparser('P3–P6 — pretraga u više procesa sa razmenom')
    pars.add_argument('--T', type=int, default=1, help='niti po procesu')
    pars.add_argument('--P', type=int, default=4, help='broj procesa')
    pars.add_argument('--K', type=int, default=5, help='max kandidata po razmeni')
    pars.add_argument('--E', type=int, default=3, help='razmena posle svake E-te iteracije')
    pars.add_argument('--topologija', choices=['prsten', 'zvezda', 'svima', 'nasumican'],
                      default='prsten', help='topologija komunikacije')
    pars.add_argument('--lam-nacin', choices=['iste', 'razlicite'], default='razlicite',
                      help='da li svi procesi koriste istu λ')
    args = pars.parse_args()
    finalna, rezultati = pokreni(args.B, args.N, args.k, args.M, args.lam,
                                 args.T, args.P, args.K, args.E, args.seme,
                                 args.start, args.lam_nacin, args.topologija)
    hv = hipervolumen_arhive(finalna)

    linije = ['## Po procesu']
    for r in sorted(rezultati, key=lambda x: x['rank']):
        st = r['stat']
        linije.append(f"- proces {r['rank']} (λ = {r['lam']:.2f}): evaluacije "
                      f"{st['evaluacije']}, poslato {st['poslato']}, primljeno "
                      f"{st['primljeno']} (od toga poznato {st['primljeno_odbaceno']}), "
                      f"tačaka arhive od drugih {r['od_drugih_u_arhivi']}, "
                      f"trajanje {st['ukupno_trajanje_s']}s")
    evaluacija_ukupno = sum(r['stat']['evaluacije'] for r in rezultati)
    od_drugih_ukupno = sum(r['od_drugih_u_arhivi'] for r in rezultati)
    linije.append('')
    linije.append(f'Ukupno evaluacija: {evaluacija_ukupno}; '
                  f'tačaka finalne arhive poreklom od drugih procesa: {od_drugih_ukupno}.')

    parametri = {'B': args.B, 'N': args.N, 'k': args.k, 'M': args.M,
                 'lambda-nacin': args.lam_nacin, 'T': args.T, 'P': args.P,
                 'K': args.K, 'E': args.E, 'topologija': args.topologija,
                 'seme': args.seme}
    stat_dict = {
        'ukupno_trajanje_s': round(max(r['stat']['ukupno_trajanje_s'] for r in rezultati), 3),
        'evaluacije_ukupno': evaluacija_ukupno,
        'poslato_ukupno': sum(r['stat']['poslato'] for r in rezultati),
        'primljeno_ukupno': sum(r['stat']['primljeno'] for r in rezultati),
        'od_drugih_u_finalnoj': od_drugih_ukupno,
    }
    print(f'evaluacija ukupno: {evaluacija_ukupno}, arhiva: {len(finalna)}, '
          f'hipervolumen: {hv:.5f}')
    sacuvaj_rezultate(args.izlaz, args.json,
                      f'{args.topologija}, P={args.P}, λ {args.lam_nacin} '
                      f'(seme {args.seme})',
                      parametri, stat_dict, finalna, hv, dodatak='\n'.join(linije))
    if args.izlaz:
        print('izveštaj:', args.izlaz)
