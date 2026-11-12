import json
import os
import subprocess
import sys

OVDE = os.path.dirname(os.path.abspath(__file__))
REZULTATI = os.path.join(OVDE, 'rezultati')

SEME = [42, 137, 2026]

ZAJEDNICKO = ['--B', '20000', '--N', '60', '--k', '5', '--M', '30', '--lam', '2']
PROCESNA_ARGS = ['--P', '4', '--T', '1', '--K', '5', '--E', '3']


def pokreni(oznaka, argv, semena):
    for seme in semena:
        md = os.path.join(REZULTATI, f'{oznaka}_s{seme}.md')
        js = os.path.join(REZULTATI, f'{oznaka}_s{seme}.json')
        print(f'>>> {oznaka} seme={seme} ...', flush=True)
        subprocess.run(argv + ['--seme', str(seme), '--izlaz', md, '--json', js],
                       cwd=OVDE, check=True)


def uporedi():
    linije = ['# Poređenje pokretanja (po pojedinačnim pokretanjima)', '',
              '| pokretanje | seme | trajanje (s) | evaluacije | '
              'veličina arhive | hipervolumen |',
              '|---|---|---|---|---|---|']
    redovi = []
    for f in sorted(os.listdir(REZULTATI)):
        if not (f.endswith('.json') and '_s' in f):
            continue
        with open(os.path.join(REZULTATI, f), encoding='utf-8') as fh:
            p = json.load(fh)
        oznaka, seme = f[:-5].split('_s')
        stat = p['stat']
        redovi.append((oznaka, int(seme),
                       stat.get('ukupno_trajanje_s', 0.0),
                       stat.get('evaluacije', stat.get('evaluacije_ukupno', 0)),
                       p['velicina_arhive'], p['hipervolumen']))
    for oznaka, seme, t, ev, va, hv in sorted(redovi):
        linije.append(f'| {oznaka} | {seme} | {t:.2f} | {ev:.0f} | {va} | {hv:.4f} |')
    linije.append('')
    linije.append('> Dopuniti tumačenjima poređenja iz postavke '
                  '(P2↔P1, P3↔P2, P4↔P3, P5↔P4, P6↔P4).')
    putanja = os.path.join(REZULTATI, 'uporedjeno.md')
    with open(putanja, 'w', encoding='utf-8') as f:
        f.write('\n'.join(linije))
    print('\n'.join(linije))


if __name__ == '__main__':
    os.makedirs(REZULTATI, exist_ok=True)
    piton = sys.executable

    pokreni('P1', [piton, 'sekvencijalna.py'] + ZAJEDNICKO, SEME)
    pokreni('P2', [piton, 'nitna.py', '--T', '4'] + ZAJEDNICKO, SEME)
    pokreni('P3', [piton, 'procesna.py'] + PROCESNA_ARGS
            + ['--lam-nacin', 'iste'] + ZAJEDNICKO, SEME)
    pokreni('P4', [piton, 'procesna.py'] + PROCESNA_ARGS
            + ['--lam-nacin', 'razlicite'] + ZAJEDNICKO, SEME)
    pokreni('P5', [piton, 'procesna.py'] + PROCESNA_ARGS
            + ['--lam-nacin', 'razlicite', '--topologija', 'zvezda'] + ZAJEDNICKO, SEME)
    pokreni('P6', [piton, 'procesna.py', '--P', '4', '--T', '1', '--K', '5',
                   '--E', '1', '--lam-nacin', 'razlicite'] + ZAJEDNICKO, SEME)

    uporedi()
