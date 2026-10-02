"""Entry 53d: recompute with 16 phases the regions where the 4-phase average was shown to be unreliable
(Entry 53c): all alpha_ac < 2, and alpha_ac >= 2 with S >= 3. Reuses e53c_boundary.run (cache_e43_ph16/)."""
import glob
import json
import sys
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e53c_boundary import C16, run  # noqa: E402

if __name__ == "__main__":
    C16.mkdir(exist_ok=True)
    rows = [json.loads(open(f).read()) for f in glob.glob(str(HERE / "cache_e43_ph4" / "*.json"))]
    todo = [(r["alpha"], r["aac"], r["S"]) for r in rows if r["S"] > 0 and (r["aac"] < 2 or r["S"] >= 3)]
    todo = [k for k in todo if not (C16 / f"a{k[0]}_q{k[1]}_S{k[2]}.json").exists()]
    print(f"{len(todo)} points to compute", flush=True)
    with Pool(11) as pool:
        for i, _ in enumerate(pool.imap_unordered(run, todo, chunksize=2)):
            if i % 50 == 0:
                print(i, flush=True)
    print("done", flush=True)
