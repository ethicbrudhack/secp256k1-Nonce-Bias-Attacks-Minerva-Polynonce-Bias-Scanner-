#!/usr/bin/sage
"""
bias_scan.sage -- SZYBKI detektor PRAWDZIWEGO biasu nonce k (BEZ klucza prywatnego).

Zasada: Gaussian Heuristic na kracie HNP.
    ratio = lambda1(realne) / GH_losowa      (<1 => krotki wektor => mozliwy bias)

ALE samo ratio<1 bywa FALSZYWE: przy duzym bias (male B) trywialny wektor
(0,...,0,W) i skalowanie zanizaja norme NIEZALEZNIE od danych. Dlatego liczymy
TEST KONTROLNY (null model): te same podpisy, ale z LOSOWO przetasowanymi z_i.
Tasowanie z niszczy prawdziwa relacje HNP  k = z*s^-1 + (r*s^-1)*d,  zachowujac
rozklady r,s,z. Sygnal jest PRAWDZIWY tylko gdy:
    ratio_real  <<  mediana(ratio_shuffled)      (z-score ponizej progu)

Dzieki temu odrzucamy artefakty kraty i raportujemy wylacznie realny bias.

Uzycie:
  sage bias_scan.sage --demo                      # walidacja na znanym biasie
  sage bias_scan.sage --file klucz_159072.json
  sage bias_scan.sage --dir /katalog --max-m 48 --trials 8 --z-thr 3.0
"""
from sage.all import *
import json, os, sys, time, argparse, math, random

n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
n = int(n)

# bias=256 usuniety: B=1 -> C*n ~ n^2 (512 bit), kosztowne i bez sensu
# (oznaczaloby pelna znajomosc k). Max realny bias = 224.
BIAS_LEVELS = [4, 8, 12, 16, 32, 48, 64, 80, 96, 128, 160, 192, 224]

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def modinv(a):
    return int(pow(int(a) % n, -1, n))

def gaussian_heuristic(dim, log_det):
    """Oczekiwana dlugosc najkrotszego wektora losowej kraty (wejscie: ln(det))."""
    log_gh = 0.5 * math.log(dim / (2.0 * math.pi * math.e)) + log_det / dim
    return math.exp(log_gh)

def lambda1_ratio(rs, ss, bias, msb):
    """Buduje krate HNP na (rs,ss) i zwraca ratio = min_norm(+-W) / GH.
    Zwraca None gdy brak wektora rozwiazania (ostatnia wsp. != +-W)."""
    m = len(rs)
    if m < 4:
        return None
    B = 1 << (256 - bias)
    C = max(1, n // B)
    W = B
    dim = m + 2

    if msb:
        t = rs; u = ss
    else:
        inv2b = modinv(pow(2, bias, n))
        t = [(x * inv2b) % n for x in rs]
        u = [(x * inv2b) % n for x in ss]

    M = matrix(ZZ, dim, dim)
    for i in range(m):
        M[i, i] = C * n
        M[m, i]     = (C * t[i]) % (C * n)
        M[m + 1, i] = (-C * u[i]) % (C * n)
    M[m, m] = 1
    M[m + 1, m + 1] = W
    Bred = M.LLL()

    # log(det) kraty = m*log(C*n) + log(W)   (macierz trojkatna w blokach)
    log_det = m * math.log(C * n) + math.log(W)
    gh = gaussian_heuristic(dim, log_det)

    min_norm = None
    for row in Bred.rows():
        last = int(row[m + 1])
        if abs(last) != W:
            continue                      # nie koduje rozwiazania HNP -> pomijamy
        norm = math.sqrt(sum(int(x) ** 2 for x in row))
        if norm > 0 and (min_norm is None or norm < min_norm):
            min_norm = norm
    if min_norm is None or gh <= 0:
        return None
    return min_norm / gh

def prep_all(sigs):
    """Zwraca (rs, ss) = (r*s^-1, z*s^-1) dla WSZYSTKICH podpisow."""
    rs, ss = [], []
    for r, s, z in sigs:
        if r == 0 or s == 0:
            continue
        sinv = modinv(s)
        rs.append((r * sinv) % n)
        ss.append((z * sinv) % n)
    return rs, ss

def test_level(rs, ss, zs_pool, bias, msb, trials, rng):
    """Zwraca (ratio_real, mean_null, std_null, zscore) dla danego (bias, typ).

    Null model: te same r,s, ale u_i przeliczone z PRZETASOWANYCH z (rozbicie
    relacji HNP). ratio_real << null => bias PRAWDZIWY.
    """
    ratio_real = lambda1_ratio(rs, ss, bias, msb)
    if ratio_real is None:
        return None
    null_vals = []
    m = len(rs)
    for _ in range(trials):
        zs_perm = zs_pool[:]
        rng.shuffle(zs_perm)
        ss_null = [(zs_perm[i]) % n for i in range(m)]  # z*s^-1 juz w puli
        rv = lambda1_ratio(rs, ss_null, bias, msb)
        if rv is not None:
            null_vals.append(rv)
    if not null_vals:
        return (ratio_real, None, None, None)
    mean = sum(null_vals) / len(null_vals)
    var = sum((v - mean) ** 2 for v in null_vals) / len(null_vals)
    std = math.sqrt(var)
    z = (mean - ratio_real) / std if std > 1e-12 else (99.0 if ratio_real < mean else 0.0)
    return (ratio_real, mean, std, z)

def _scan_window(rs, ss, args, rng, win_info, best_holder):
    """Testuje wszystkie (bias,typ) na JEDNYM oknie (rs,ss). Aktualizuje best_holder."""
    zs_pool = ss[:]  # pula z*s^-1 do tasowania (null)
    for bias in BIAS_LEVELS:
        need = max(4, int(256 / bias * 4 / 3))
        if len(rs) < need:
            continue
        for msb in (True, False):
            out = test_level(rs, ss, zs_pool, bias, msb, args.trials, rng)
            if out is None:
                continue
            ratio_real, mean, std, z = out
            genuine = (z is not None and z >= args.z_thr and ratio_real < args.ratio_thr)
            rec = {"bias": bias, "type": "MSB" if msb else "LSB",
                   "ratio": round(ratio_real, 4),
                   "null_mean": round(mean, 4) if mean is not None else None,
                   "zscore": round(z, 2) if z is not None else None,
                   "genuine": bool(genuine), **win_info}
            if genuine and (best_holder[0] is None or ratio_real < best_holder[0]["ratio"]):
                best_holder[0] = rec
            if args.verbose:
                flag = " <== PRAWDZIWY" if genuine else ""
                log(f"    {win_info.get('window_desc','')} bias={bias:3d} {rec['type']} "
                    f"ratio={rec['ratio']:.3f} null={rec['null_mean']} z={rec['zscore']}{flag}")

def scan_key(key_xy, sigs, args, rng):
    """Skanuje klucz. Zwraca najlepszy PRAWDZIWY sygnal albo None.

    Domyslnie: pierwsze args.max_m podpisow (szybko).
    --full-grind: PRZESUWANE okna po args.max_m przez CALY zbior podpisow
                  (lapie bias lokalny/okresowy w dowolnym miejscu).
    """
    rs_all, ss_all = prep_all(sigs)
    total = len(rs_all)
    if total < 4:
        return None
    m = min(args.max_m, total)
    best_holder = [None]

    if not args.full_grind:
        # jedno okno z poczatku
        _scan_window(rs_all[:m], ss_all[:m], args, rng,
                     {"window_start": 0, "window_desc": ""}, best_holder)
        return best_holder[0]

    # PELNE PRZEMIELANIE: przesuwane okna przez caly zbior
    step = max(1, int(m * args.step))
    nwin = (total - m) // step + 1
    if args.max_windows and nwin > args.max_windows:
        nwin = args.max_windows
    if args.verbose:
        log(f"  [grind] {total} podpisow -> {nwin} okien po m={m} (krok={step})")
    for w in range(nwin):
        a = w * step
        b = a + m
        if b > total:
            break
        _scan_window(rs_all[a:b], ss_all[a:b], args, rng,
                     {"window_start": a, "window_desc": f"okno#{w}@{a}"}, best_holder)
        if (w + 1) % 200 == 0:
            log(f"  [grind] okno {w+1}/{nwin} (przemielono {b}/{total} podpisow)")
    return best_holder[0]

def demo():
    """Generuje podpisy z ZNANYM biasem i sprawdza czy detektor je lapie."""
    try:
        from ecdsa import SigningKey, SECP256k1
    except Exception:
        log("DEMO wymaga biblioteki 'ecdsa'.  pip install ecdsa"); return
    import hashlib as hl
    log("DEMO: podpisy z ZNANYM biasem MSB, walidacja detektora.")
    rng = random.Random(1234)
    sk = SigningKey.generate(curve=SECP256k1)

    class A:  # atrapa args
        max_m = 48; trials = 8; z_thr = 3.0; ratio_thr = 0.8; verbose = True
        full_grind = False; step = 0.5; max_windows = 0
    args = A()

    for bias_bits in (0, 4, 8):  # 0 = kontrola (zdrowe nonce, ma NIC nie wykryc)
        sigs = []
        for _ in range(60):
            msg = rng.getrandbits(256).to_bytes(32, "big")
            z = int(hl.sha256(msg).hexdigest(), 16)
            if bias_bits == 0:
                k = rng.getrandbits(256) % (n - 1) + 1
            else:
                k = rng.getrandbits(256 - bias_bits) or 1
            sig = sk.sign(msg, k=k, hashfunc=hl.sha256)
            r = int.from_bytes(sig[:32], "big"); s = int.from_bytes(sig[32:], "big")
            sigs.append((r, s, z))
        label = "ZDROWE (kontrola)" if bias_bits == 0 else f"bias MSB={bias_bits}b"
        log(f"--- {label} ---")
        best = scan_key("demo", sigs, args, rng)
        log(f"  => wynik: {best}")

def main():
    ap = argparse.ArgumentParser(description="Szybki detektor prawdziwego biasu nonce (HNP GH + null model).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--file")
    src.add_argument("--dir")
    src.add_argument("--demo", action="store_true")
    ap.add_argument("--out", default="bias_report.json")
    ap.add_argument("--max-m", type=int, default=48, help="podpisow na krate. domyslnie 48")
    ap.add_argument("--trials", type=int, default=8, help="prob null (permutacje z). domyslnie 8")
    ap.add_argument("--z-thr", type=float, default=3.0, help="min z-score dla 'prawdziwy'. dom. 3.0")
    ap.add_argument("--ratio-thr", type=float, default=0.8, help="max ratio dla 'prawdziwy'. dom. 0.8")
    ap.add_argument("--full-grind", action="store_true",
                    help="PRZESUWANE okna przez CALY zbior podpisow (lapie bias lokalny)")
    ap.add_argument("--step", type=float, default=0.5,
                    help="krok okna jako ulamek m przy --full-grind (0<step<=1). dom. 0.5")
    ap.add_argument("--max-windows", type=int, default=0,
                    help="limit okien na klucz przy --full-grind (0=wszystkie)")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--resume", action="store_true",
                    help="wznow: pomijaj pliki juz przerobione (wg <out>.progress)")
    ap.add_argument("--progress-every", type=int, default=25,
                    help="loguj postep co N plikow. dom. 25")
    ap.add_argument("--min-size-mb", type=float, default=0.0,
                    help="przetwarzaj TYLKO pliki >= X MB (bieg na duze pliki). dom. 0")
    ap.add_argument("--max-size-mb", type=float, default=0.0,
                    help="przetwarzaj TYLKO pliki < X MB (bieg na male pliki). 0=bez limitu")
    ap.add_argument("--biggest-first", action="store_true",
                    help="sortuj pliki wg rozmiaru MALEJACO (najpierw najwieksze)")
    args = ap.parse_args()

    rng = random.Random(42)
    if args.demo:
        demo(); return

    import glob
    files = [args.file] if args.file else sorted(glob.glob(os.path.join(args.dir, "*.json")))
    if not files:
        log("Brak plikow."); sys.exit(1)

    # --- Filtr rozmiaru pliku (bieg na duze / male pliki osobno) ---
    if args.min_size_mb > 0 or args.max_size_mb > 0 or args.biggest_first:
        lo = args.min_size_mb * 1048576
        hi = args.max_size_mb * 1048576
        sized = []
        for f in files:
            try:
                sz = os.path.getsize(f)
            except OSError:
                continue
            if sz < lo:
                continue
            if hi > 0 and sz >= hi:
                continue
            sized.append((sz, f))
        if args.biggest_first:
            sized.sort(key=lambda x: -x[0])   # najpierw najwieksze
        files = [f for _, f in sized]
        log(f"Po filtrze rozmiaru: {len(files)} plikow "
            f"(min={args.min_size_mb}MB, max={args.max_size_mb or 'inf'}MB, "
            f"biggest_first={args.biggest_first})")
        if not files:
            log("Brak plikow po filtrze."); sys.exit(1)

    # --- Zapis na biezaco (atomowo) + wznawianie ---
    prog_fn = args.out + ".progress"
    hits_fn = args.out + ".hits.jsonl"   # kazdy hit = jedna linia (nigdy nie tracimy)

    def atomic_dump(obj, path):
        tmp = path + ".tmp"
        json.dump(obj, open(tmp, "w"), indent=2)
        os.replace(tmp, path)

    report = []
    done = set()
    if args.resume:
        if os.path.exists(prog_fn):
            done = set(l.strip() for l in open(prog_fn) if l.strip())
            log(f"WZNAWIANIE: {len(done)} plikow juz przerobionych, pomijam.")
        if os.path.exists(hits_fn):
            for line in open(hits_fn):
                line = line.strip()
                if line:
                    try: report.append(json.loads(line))
                    except Exception: pass
            log(f"WZNAWIANIE: wczytano {len(report)} wczesniejszych trafien.")

    prog_f = open(prog_fn, "a")   # dopisujemy nazwy przerobionych plikow
    hits_f = open(hits_fn, "a")   # dopisujemy trafienia (JSON Lines)

    t0 = time.time()
    n_files = len(files)
    for idx, fn in enumerate(files, 1):
        if args.resume and fn in done:
            continue
        try:
            data = json.load(open(fn))
        except Exception as e:
            log(f"BLAD {fn}: {e}")
            prog_f.write(fn + "\n"); prog_f.flush()
            continue
        for key_xy, sig_list in data.items():
            sigs = [(int(r, 16), int(s, 16), int(z, 16)) for r, s, z in sig_list]
            if len(sigs) < 4:
                continue
            if args.verbose:
                log(f"Klucz {key_xy[:24]}.. ({len(sigs)} sig)")
            best = scan_key(key_xy, sigs, args, rng)
            if best:
                rec = {"key": key_xy, "m_total": len(sigs), "file": os.path.basename(fn), **best}
                report.append(rec)
                # 1) natychmiast dopisz do .hits.jsonl (crash-safe)
                hits_f.write(json.dumps(rec) + "\n"); hits_f.flush()
                # 2) odswiez posortowany raport zbiorczy
                atomic_dump(sorted(report, key=lambda x: x["ratio"]), args.out)
                log(f"  *** PRAWDZIWY BIAS: {key_xy[:24]}.. bias~{best['bias']}b "
                    f"{best['type']} ratio={best['ratio']} z={best['zscore']} start={best.get('window_start')}")
        # plik zaliczony -> zapisz postep
        prog_f.write(fn + "\n"); prog_f.flush()
        if idx % args.progress_every == 0:
            el = time.time() - t0
            log(f"[postep] {idx}/{n_files} plikow, trafien: {len(report)}, {el:.0f}s "
                f"({el/idx:.1f}s/plik)")

    prog_f.close(); hits_f.close()
    atomic_dump(sorted(report, key=lambda x: x["ratio"]), args.out)
    log(f"KONIEC ({time.time()-t0:.0f}s). Kluczy z PRAWDZIWYM biasem: {len(report)} -> {args.out}")

if __name__ == "__main__":
    main()
