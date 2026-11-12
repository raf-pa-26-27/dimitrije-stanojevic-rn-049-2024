# Prvi projekat — optimizacija sekvenci genoma (Pareto, niti, procesi)

Optimizacija 32-aminokiselinske sekvence RBD domena spike proteina
SARS-CoV-2 po tri kriterijuma (zaraznost, virulentnost, imunoevazivnost)
Pareto evolutivnom pretragom, u tri verzije: sekvencijalnoj, sa više
niti i u više procesa sa razmenom kandidata.

## Pokretanje

Sve komande se izvršavaju iz direktorijuma `prvi/` (tamo su i
`oceni.py`, `blosum62.json`, `anchors.json`).

```bash
# provere ispravnosti
python testovi.py

# P1 — jednonitna pretraga
python sekvencijalna.py --B 20000 --N 60 --k 5 --M 30 --lam 2 --seme 42

# P2 — više niti
python nitna.py --B 20000 --T 4 --seme 42

# P3–P6 — više procesa sa razmenom (prsten; za P5 dodati --topologija zvezda)
python procesna.py --B 20000 --P 4 --T 1 --K 5 --E 3 --topologija prsten \
       --lam-nacin razlicite --seme 42

# sva merenja P1–P6 na tri semena + tabela po pojedinačnim pokretanjima
python izmeri.py
```

Parametri (B, N, k, M, λ, T, P, K, E, topologija, seme) zadaju se
spolja, argumentima komandne linije. Prag uzastopnih neuspelih
generisanja G = 100 je konstanta u kodu.

## Pravila pretrage

- Početna populacija: početna sekvenca i N−1 međusobno različitih
  mutanata (1..k zamena pozicija); evaluacije ulaze u budžet.
- Novi kandidat: dva različita nasumična roditelja iz populacije,
  ukrštanje u jednoj tački, mutacija na s ~ Poisson(λ) pozicija
  (zamenom drugim slovom). Sve slučajne odluke iz jednog numpy
  generatora (`default_rng(seme)`), kreiranog jednom po pretragi.
- Duplikati (sekvenca već u evidenciji ili evaluacija u toku) se
  preskačaju; ne troše budžet i ne ulaze u M.
- Iteracija prikupi M novih ocenjenih kandidata; zatim se populacija,
  arhiva i novi kandidati spoje, rangiraju u Pareto slojeve, ceo prvi
  sloj postaje nova arhiva, a populacija se puni po rangu, unutar sloja
  po većoj meri razređenosti, do N kandidata.
- Kraj: potrošen budžet B ili G uzastopnih neuspelih generisanja.
- Konačan rezultat je Pareto arhiva (sekvence i tri ocene).

## Organizacija niti (P2)

T radnih niti prikuplja kandidate unutar jedne iteracije; tokom
prikupljanja populacija i arhiva se ne menjaju. Deljene strukture:
evidencija (sekvenca → None tokom evaluacije → ocene), brojači
(započeti pozivi evaluatora, rezervacije u iteraciji, uzastopni neuspeli
pokušaji) i populacija/arhiva — sve pod jednom bravom (Condition);
pozivi generatora slučajnih brojeva pod zasebnom bravom. Rezervacija
mesta u evidenciji (None) vrši se pre evaluacije i odmah troši budžet
(računaju se započeti pozivi evaluatora), čime je isključena i dvostruka
evaluacija iste sekvence. Evaluacija se izvršava van brava; nit koja
završi kandidata odmah nastavlja sa sledećim, a posao se ne deli
unapred pa nema izgladnjivanja. Kad je M kandidata rezervisano,
iteracija se zatvara; glavna nit sačeka da se završe sve započete
evaluacije (nijedna započeta evaluacija ne ostaje van iteracije), pa
rangira spoj kroz poseban thread pool (isti broj niti kao radne niti,
odvojen od njih) i ažurira populaciju i arhivu pre sledeće iteracije,
dok radne niti čekaju. Ako pretraga prekine usred iteracije (budžet ili
prag G), prikupljeni kandidati se i dalje spajaju kao u običnoj
iteraciji.

## Organizacija procesa i protokol razmene (P3–P6)

Svaki proces pokreće celu pretragu sopstvenim stanjem (populacija,
arhiva, evidencija, slučajni tok) i budžetom B/P. λ je ista za sve
procese (--lam-nacin iste) ili linearno raspoređena od 0.5 do 4.0
(--lam-nacin razlicite), tako da neki procesi pretražuju oprezno, a
 drugi široko.

- Jedan ulazni red (`multiprocessing.Queue`, maxsize 500) po procesu.
- Šta se šalje: najviše K kandidata iz Pareto arhive sa najvećom merom
  razređenosti, sa svim ocenama; ako je arhiva manja od K, šalju se svi.
  Format poruke: `{'od': rank, 'kandidati': [(sekvenca, (t, s, e)), ...]}`.
- Kada: posle svake E-te završene lokalne iteracije.
- Kome: svim susedima po topologiji (`--topologija prsten | zvezda |
  svima | nasumican`, parametar programa).
- Prijem: na početku svake lokalne iteracije, neblokirajuće
  (`get_nowait` dok se red ne isprazni). Poznate sekvence se odbacuju;
  nove se upisuju u evidenciju sa primljenim ocenama, bez poziva
  evaluatora i bez trošenja budžeta, i učestvuju u ažuriranju
  populacije i arhive u toj iteraciji kao lokalni kandidati (ne ulaze
  u M). Prijem poruke ne resetuje brojač neuspelih pokušaja (resetuje
  ga samo lokalno pronađeni novi kandidat).
- Proces čiji brojač dostigne G se završava; drugi nastavljaju. Poruke
  upućene završenom procesu se odbacuju kad se njegov red napuni
  (slanje ide sa kratkim tajmautom pa pošiljalac nikad ne blokira).
- Prikupljanje rezultata: porukama — poseban red za rezultate ka
  glavnom procesu. Glavni proces čita sa tajmautom i proverom da li su
  procesi živi, pa nikad ne čeka zauvek na pokvaren proces; na kraju
  spaja lokalne arhive, uklanja duplikate i dominirane kandidate istom
  funkcijom dominacije kao u lokalnoj pretrazi.

## Provere

`testovi.py` pokriva: dominaciju, Pareto slojeve, meru razređenosti,
ukrštanje, mutaciju, početnu populaciju, ažuriranje, determinizam po
semenu, mali budžet, oba uslova završetka, nitnu verziju (bez duplikata,
budžet tačan, bez ostavljenih None ocena), paralelno rangiranje kroz
pool (isti slojevi kao sekvencijalno), topologije, izbor kandidata za
slanje i procesnu verziju (podela budžeta, spajanje arhivi).

## Okruženje

<!-- TODO: popuniti svoj hardver i verzije pre predaje -->
- Procesor: (popuniti — model i broj jezgara)
- Python: 3.x, biblioteke: numpy (standardna biblioteka u ostalom)
