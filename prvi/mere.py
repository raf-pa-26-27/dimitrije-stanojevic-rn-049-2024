import json
import time


def hv2d(tacke, ref):
    povrsina, najbolji_y = 0.0, ref[1]
    for x, y in sorted(tacke, reverse=True):
        if y > najbolji_y:
            povrsina += (x - ref[0]) * (y - najbolji_y)
            najbolji_y = y
    return povrsina


def hv3d(tacke, ref):
    tacke = sorted(tacke, key=lambda t: t[2], reverse=True)
    zapremina = 0.0
    for i, (_, _, z) in enumerate(tacke):
        sledeci_z = tacke[i + 1][2] if i + 1 < len(tacke) else ref[2]
        visina = z - sledeci_z
        if visina > 0:
            sloj = [(x, y) for x, y, _ in tacke[: i + 1]]
            zapremina += hv2d(sloj, ref[:2]) * visina
    return zapremina


def normalizuj(tacke, lo, hi):
    return [tuple((v - l) / (h - l) for v, l, h in zip(t, lo, hi)) for t in tacke]


OCENE_LO = (0.0, 0.0, 0.0)
OCENE_HI = (10.0, 10.0, 10.0)
HV_REF = (-0.1, -0.1, -0.1)


def hipervolumen_arhive(arhiva):
    if not arhiva:
        return 0.0
    ocene = [tuple(k.ocene) for k in arhiva]
    front = normalizuj(ocene, OCENE_LO, OCENE_HI)
    return hv3d(front, HV_REF)


class Statistika:

    def __init__(self):
        self.pocetak = time.perf_counter()
        self._t_iter = None
        self.trajanja_iteracija = []
        self.evaluacije = 0
        self.poslato = 0
        self.primljeno = 0
        self.primljeno_odbaceno = 0

    def iteracija_pocela(self):
        self._t_iter = time.perf_counter()

    def iteracija_zavrsena(self):
        if self._t_iter is not None:
            self.trajanja_iteracija.append(time.perf_counter() - self._t_iter)
            self._t_iter = time.perf_counter()

    def ukupno_trajanje(self):
        return time.perf_counter() - self.pocetak

    def prosecna_iteracija(self):
        if not self.trajanja_iteracija:
            return 0.0
        return sum(self.trajanja_iteracija) / len(self.trajanja_iteracija)

    def to_dict(self):
        return {
            'ukupno_trajanje_s': round(self.ukupno_trajanje(), 3),
            'prosecna_iteracija_s': round(self.prosecna_iteracija(), 4),
            'broj_iteracija': len(self.trajanja_iteracija),
            'evaluacije': self.evaluacije,
            'poslato': self.poslato,
            'primljeno': self.primljeno,
            'primljeno_odbaceno': self.primljeno_odbaceno,
        }


def rasponi_ocena(arhiva):
    if not arhiva:
        return None
    return [
        (min(k.ocene[i] for k in arhiva), max(k.ocene[i] for k in arhiva))
        for i in range(3)
    ]


def izvestaj_markdown(naslov, parametri, stat_dict, arhiva, hv, dodatak=''):
    linije = [f'# {naslov}', '']
    linije.append('## Parametri')
    for kljuc, vrednost in parametri.items():
        linije.append(f'- {kljuc}: {vrednost}')
    linije.append('')
    linije.append('## Statistika')
    for kljuc, vrednost in stat_dict.items():
        linije.append(f'- {kljuc}: {vrednost}')
    rr = rasponi_ocena(arhiva)
    if rr:
        for ime, (lo, hi) in zip(('zaraznost', 'virulentnost', 'imunoevazivnost'), rr):
            linije.append(f'- raspon {ime}: [{lo:.2f}, {hi:.2f}]')
    linije.append(f'- veličina arhive: {len(arhiva)}')
    linije.append(f'- hipervolumen: {hv:.5f}')
    if dodatak:
        linije.append('')
        linije.append(dodatak)
    linije.append('')
    linije.append('## Konačna arhiva')
    for k in sorted(arhiva, key=lambda x: -x.ocene[0]):
        t, s, e = k.ocene
        linije.append(f'- {k.seq} — zaraznost {t:.2f}, virulentnost {s:.2f}, '
                      f'imunoevazivnost {e:.2f}')
    linije.append('')
    return '\n'.join(linije)


def sacuvaj_rezultate(izlaz_md, izlaz_json, naslov, parametri, stat_dict,
                      arhiva, hv, dodatak=''):
    if izlaz_md:
        with open(izlaz_md, 'w', encoding='utf-8') as f:
            f.write(izvestaj_markdown(naslov, parametri, stat_dict, arhiva, hv, dodatak))
    if izlaz_json:
        with open(izlaz_json, 'w', encoding='utf-8') as f:
            json.dump({
                'naslov': naslov,
                'parametri': parametri,
                'stat': stat_dict,
                'hipervolumen': hv,
                'velicina_arhive': len(arhiva),
                'arhiva': [(k.seq, list(k.ocene)) for k in arhiva],
            }, f, ensure_ascii=False, indent=1)
