import json, sys, time, resource, torch
sys.path.insert(0, "<repo>")
from model.a.train import TrainConfig, _load_all_windows, _window_loss, _emissions_for
from model.a.encoder import CandidateA
torch.set_num_threads(1)
cfg = TrainConfig.from_dict(json.load(open("<scratch>/run/smoke_config.json")))
t0=time.perf_counter(); ex, stats = _load_all_windows(cfg); t_load=time.perf_counter()-t0
print(f"load: {len(ex)} windows in {t_load:.1f}s; rss_MiB={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f}")
ex = sorted(ex, key=lambda e: e.n)
picks = [ex[len(ex)//10], ex[len(ex)//2], ex[9*len(ex)//10], ex[-1]]
model = CandidateA(); dev = torch.device("cpu"); dt = torch.float64
for e in picks:
    t0=time.perf_counter()
    with torch.no_grad(): em = _emissions_for(model, e, dev, dt)
    t_enc=time.perf_counter()-t0
    t0=time.perf_counter(); loss = _window_loss(model, e, dev, dt); t_fwd=time.perf_counter()-t0
    t0=time.perf_counter(); loss.backward(); t_bwd=time.perf_counter()-t0
    print(f"n={e.n:6d} enc_fwd={t_enc:.3f}s loss_fwd(enc+chain)={t_fwd:.3f}s bwd={t_bwd:.3f}s "
          f"total={t_fwd+t_bwd:.3f}s per_kb={(t_fwd+t_bwd)/(e.n/1000):.3f}s rss_MiB={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f}")
