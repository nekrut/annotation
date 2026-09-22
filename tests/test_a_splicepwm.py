"""Tests for the splice-site PWM of section 3.5 increment 2
(``model.a.splicepwm``): the estimator counts the junction context it says it
counts and normalises against the window background; the bias is the same
score computed position by position, honours the two edge conventions, is
additive with the dinucleotide bias and the hard mask, and survives a JSON
round trip. Skips without torch, like the other candidate-A tensor suites.
"""
import importlib.util
import json
import math
import random
import tempfile
import unittest
from pathlib import Path

try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ModuleNotFoundError:
    _HAS_TORCH = False

if _HAS_TORCH and importlib.util.find_spec("model.a.splicepwm") is not None:
    from model.a import pooled, splicepwm
    from model.a.encoder import DecoderParams
    HAS = True
else:
    HAS = False


class _Ex:
    """The two attributes :func:`model.a.splicepwm.fit` reads of a
    ``WindowExample``."""

    def __init__(self, window, intron_ranges):
        self.window = window
        self.intron_ranges = intron_ranges


def _uniform_pwm(donor_span=(-1, 2), acceptor_span=(-2, 1), **kw):
    """A PWM whose columns are set by hand, for exact-score assertions."""
    nd = donor_span[1] - donor_span[0]
    na = acceptor_span[1] - acceptor_span[0]
    return splicepwm.SplicePWM(
        donor=tuple(tuple(float(10 * j + i) for i in range(4)) for j in range(nd)),
        acceptor=tuple(tuple(float(-10 * j - i) for i in range(4)) for j in range(na)),
        background=(0.25, 0.25, 0.25, 0.25),
        donor_span=donor_span, acceptor_span=acceptor_span, **kw)


def _naive_scores(x, pwm):
    """Reference: the PWM score of every boundary, by definition."""
    D = pooled.CHANNEL_ORDER.index("donor") if hasattr(pooled, "CHANNEL_ORDER") else None
    out = {"donor": [], "acceptor": []}
    for name, mat, span in (("donor", pwm.donor, pwm.donor_span),
                            ("acceptor", pwm.acceptor, pwm.acceptor_span)):
        for t in range(len(x)):
            s = 0.0
            for j in range(span[1] - span[0]):
                p = t + span[0] + j
                if 0 <= p < len(x):
                    i = splicepwm.BASE_INDEX.get(x[p].upper())
                    if i is not None:
                        s += mat[j][i]
            out[name].append(s)
    return out


@unittest.skipUnless(HAS, "needs torch and model.a.splicepwm")
class Estimation(unittest.TestCase):
    def test_counts_exactly_the_declared_context(self):
        # One intron [4, 9): donor boundary 4, acceptor boundary 9.
        x = "AAAAGTAAGCCCCCC"
        ex = _Ex(x, [(4, 9)])
        pwm = splicepwm.fit([ex], donor_span=(-1, 2), acceptor_span=(-2, 1),
                            pseudocount=0.0, min_sites=1)
        self.assertEqual(pwm.sites, (1, 1))
        # donor columns: x[3]='A', x[4]='G', x[5]='T'
        for j, base in enumerate("AGT"):
            col = pwm.donor[j]
            i = splicepwm.BASE_INDEX[base]
            # frequency 1 for the observed base, 0 (-inf in log) elsewhere
            self.assertAlmostEqual(col[i], math.log(1.0 / pwm.background[i]))
            for k in range(4):
                if k != i:
                    self.assertEqual(col[k], float("-inf"))
        # acceptor columns: x[7]='A', x[8]='G', x[9]='C'
        for j, base in enumerate("AGC"):
            i = splicepwm.BASE_INDEX[base]
            self.assertAlmostEqual(pwm.acceptor[j][i], math.log(1.0 / pwm.background[i]))

    def test_background_is_the_window_composition(self):
        x = "AACCCGGGGTTTTTTT"  # 2 A, 3 C, 5 G, ... whatever it is
        ex = _Ex(x, [(4, 9)])
        pwm = splicepwm.fit([ex], min_sites=1)
        for base, i in splicepwm.BASE_INDEX.items():
            self.assertAlmostEqual(pwm.background[i], x.count(base) / len(x))

    def test_uninformative_column_scores_about_zero(self):
        # Junction contexts drawn from the same composition as the windows:
        # every column's log-odds must be near 0, so the PWM adds nothing.
        rng = random.Random(7)
        exs = []
        for _ in range(400):
            w = "".join(rng.choice("ACGT") for _ in range(200))
            exs.append(_Ex(w, [(60, 140)]))
        pwm = splicepwm.fit(exs, min_sites=1)
        for col in pwm.donor + pwm.acceptor:
            for v in col:
                self.assertLess(abs(v), 0.35)

    def test_non_acgt_and_out_of_window_columns_are_not_counted(self):
        ex = _Ex("NNGTAAGNNCC", [(2, 7)])
        pwm = splicepwm.fit([ex], donor_span=(-4, 2), acceptor_span=(-1, 4),
                            pseudocount=1.0, min_sites=1)
        # donor offsets -4 and -3 fall off the window, -2/-1 are 'N': the four
        # columns saw no base at all, so each is the pure pseudocount --
        # log(0.25 / background), i.e. no count information, not a flat score.
        flat = [math.log(0.25 / b) for b in pwm.background]
        for j in range(4):
            for i in range(4):
                self.assertAlmostEqual(pwm.donor[j][i], flat[i])
        # offset 0 saw 'G' and must differ from the uninformative column.
        self.assertNotAlmostEqual(pwm.donor[4][splicepwm.BASE_INDEX["G"]], flat[
            splicepwm.BASE_INDEX["G"]])

    def test_min_sites_refuses_a_thin_matrix(self):
        with self.assertRaises(ValueError):
            splicepwm.fit([_Ex("ACGTACGTAC", [(3, 8)])], min_sites=2)

    def test_empty_background_refused(self):
        with self.assertRaises(ValueError):
            splicepwm.fit([_Ex("NNNNNNNNNN", [(3, 8)])], min_sites=1)


@unittest.skipUnless(HAS, "needs torch and model.a.splicepwm")
class BalancedEstimation(unittest.TestCase):
    """:func:`model.a.splicepwm.fit_grouped`: one weight per source, not per
    junction."""

    def test_single_group_reproduces_the_pooled_matrix(self):
        ex = [_Ex("AAAAGTAAGCCCCCC", [(4, 9)]), _Ex("CCGCGTGAGTTTTAA", [(4, 10)])]
        pooled_pwm = splicepwm.fit(ex, donor_span=(-1, 2), acceptor_span=(-2, 1),
                                   min_sites=1)
        bal = splicepwm.fit_grouped([("one", ex)], donor_span=(-1, 2),
                                    acceptor_span=(-2, 1), min_sites=1)
        self.assertEqual(bal.weighting, "balanced")
        self.assertEqual(pooled_pwm.weighting, "pooled")
        self.assertEqual(bal.sites, pooled_pwm.sites)
        for got, want in zip(bal.background, pooled_pwm.background):
            self.assertAlmostEqual(got, want)
        for name in ("donor", "acceptor"):
            for cb, cp in zip(getattr(bal, name), getattr(pooled_pwm, name)):
                for a, b in zip(cb, cp):
                    self.assertAlmostEqual(a, b)

    def test_a_thin_source_counts_as_much_as_a_thick_one(self):
        # 100 junctions with a C at the donor column against 1 with a G:
        # pooled is decided by the majority, balanced splits the column.
        thick = [_Ex("ACGTCACGTACGTAT", [(4, 9)]) for _ in range(100)]
        thin = [_Ex("ACGTGACGTACGTAT", [(4, 9)])]
        span = (0, 1)
        pooled_pwm = splicepwm.fit(thick + thin, donor_span=span,
                                   acceptor_span=span, min_sites=1)
        bal = splicepwm.fit_grouped([("thick", thick), ("thin", thin)],
                                    donor_span=span, acceptor_span=span,
                                    min_sites=1)
        c, g = splicepwm.BASE_INDEX["C"], splicepwm.BASE_INDEX["G"]
        pooled_gap = pooled_pwm.donor[0][c] - pooled_pwm.donor[0][g]
        bal_gap = bal.donor[0][c] - bal.donor[0][g]
        self.assertGreater(pooled_gap, 3.0)
        # Each source contributes (n+1)/(n+4) to its own base, so the two
        # smoothed frequencies average to a column that still prefers C, but
        # the thin source has pulled the preference down by most of it.
        self.assertGreater(bal_gap, 0.0)
        self.assertLess(bal_gap, pooled_gap / 3.0)
        thick_f = (100 + 1.0) / (100 + 4.0)
        thin_f = (0 + 1.0) / (1 + 4.0)
        want = (thick_f + thin_f) / 2.0
        self.assertAlmostEqual(bal.donor[0][c],
                               math.log(want / bal.background[c]))

    def test_background_averages_the_sources_not_their_bases(self):
        # One long AT-rich source and one short GC-rich one: the balanced
        # background is the mean of the two compositions, halfway between.
        at = [_Ex("AAAATTTTAAAATTCG", [(4, 9)]) for _ in range(20)]   # A8 T6 C1 G1
        gc = [_Ex("CCCCGGGGCCCCGGAT", [(4, 9)])]                       # C8 G6 A1 T1
        bal = splicepwm.fit_grouped([("at", at), ("gc", gc)], min_sites=1)
        want = {"A": (8 + 1) / 32, "T": (6 + 1) / 32,
                "C": (1 + 8) / 32, "G": (1 + 6) / 32}
        for base, value in want.items():
            self.assertAlmostEqual(bal.background[splicepwm.BASE_INDEX[base]], value)
        pooled_pwm = splicepwm.fit(at + gc, min_sites=1)
        self.assertGreater(pooled_pwm.background[splicepwm.BASE_INDEX["A"]], 0.4)

    def test_site_counts_land_in_the_source_entries(self):
        thick = [_Ex("ACGTCACGTACGTAT", [(4, 9), (10, 13)]) for _ in range(3)]
        thin = [_Ex("ACGTGACGTACGTAT", [(4, 9)])]
        bal = splicepwm.fit_grouped(
            [("thick", thick), ("thin", thin)], min_sites=1,
            sources=[{"name": "thick", "dev_seqids": []},
                     {"name": "thin", "dev_seqids": []},
                     {"name": "absent", "dev_seqids": []}])
        by_name = {s["name"]: s for s in bal.sources}
        self.assertEqual(by_name["thick"]["train_sites"], 6)
        self.assertEqual(by_name["thin"]["train_sites"], 1)
        self.assertNotIn("train_sites", by_name["absent"])
        self.assertEqual(bal.sites, (7, 7))

    def test_a_source_with_no_junction_is_refused(self):
        ok = [_Ex("ACGTCACGTACGTAT", [(4, 9)])]
        with self.assertRaises(ValueError):
            splicepwm.fit_grouped([("ok", ok), ("empty", [])], min_sites=1)
        with self.assertRaises(ValueError):
            splicepwm.fit_grouped([], min_sites=1)
        # min_group_sites is the same guard, made explicit.
        with self.assertRaises(ValueError):
            splicepwm.fit_grouped([("ok", ok)], min_sites=1, min_group_sites=2)

    def test_weighting_survives_the_round_trip_and_defaults_to_pooled(self):
        bal = splicepwm.fit_grouped([("one", [_Ex("ACGTCACGTACGTAT", [(4, 9)])])],
                                    min_sites=1)
        d = bal.to_json()
        self.assertEqual(d["weighting"], "balanced")
        self.assertEqual(splicepwm.SplicePWM.from_json(d).weighting, "balanced")
        del d["weighting"]  # a matrix written before the field existed
        self.assertEqual(splicepwm.SplicePWM.from_json(d).weighting, "pooled")


class Bias(unittest.TestCase):
    def test_matches_the_position_by_position_score(self):
        rng = random.Random(3)
        x = "".join(rng.choice("ACGTN") for _ in range(120))
        pwm = _uniform_pwm(donor_span=(-3, 6), acceptor_span=(-20, 3))
        b = splicepwm.bias(x, pwm)
        want = _naive_scores(x, pwm)
        D = pooled.CHANNEL_ORDER.index("donor")
        A = pooled.CHANNEL_ORDER.index("acceptor")
        for t in range(len(x)):
            self.assertAlmostEqual(float(b[D, t]), want["donor"][t], places=9)
            self.assertAlmostEqual(float(b[A, t]), want["acceptor"][t], places=9)

    def test_only_the_splice_rows_are_touched(self):
        pwm = _uniform_pwm()
        b = splicepwm.bias("ACGTACGTACGT", pwm)
        D = pooled.CHANNEL_ORDER.index("donor")
        A = pooled.CHANNEL_ORDER.index("acceptor")
        for c in range(b.shape[0]):
            if c not in (D, A):
                self.assertTrue(bool(torch.all(b[c] == 0)), f"channel {c} not zero")

    def test_batch_pads_past_each_window(self):
        pwm = _uniform_pwm()
        xs = ["ACGTACGT", "TTTT", ""]
        b = splicepwm.bias_batch(xs, pwm, L=12)
        for i, x in enumerate(xs):
            self.assertTrue(bool(torch.all(b[i, :, len(x):] == 0)))
            if x:
                one = splicepwm.bias(x, pwm)
                self.assertTrue(bool(torch.allclose(b[i, :, :len(x)], one)))
        with self.assertRaises(ValueError):
            splicepwm.bias_batch(["A" * 13], pwm, L=12)

    def test_additive_with_the_dinucleotide_bias_and_the_mask(self):
        torch.manual_seed(0)
        dec = DecoderParams()
        with torch.no_grad():
            for p in dec.parameters():
                p.normal_(0.0, 1.0)
        x = "ACGTGTAAGCCCCCCCCCCAGGTACGT"
        pwm = _uniform_pwm(donor_span=(-3, 6), acceptor_span=(-20, 3))
        m = pooled.motif_bias(x, dec, canonical=True)
        p = splicepwm.bias(x, pwm)
        total = m + p
        # A masked site stays -inf however large the PWM score is: the mask is
        # the support, the PWM only ranks what the mask leaves.
        D = pooled.CHANNEL_ORDER.index("donor")
        masked = [t for t in range(len(x)) if bool(torch.isinf(m[D, t]))]
        self.assertTrue(masked)
        for t in masked:
            self.assertTrue(bool(torch.isinf(total[D, t])) and float(total[D, t]) < 0)
        unmasked = [t for t in range(len(x)) if not bool(torch.isinf(m[D, t]))]
        for t in unmasked:
            self.assertAlmostEqual(float(total[D, t]),
                                   float(m[D, t]) + float(p[D, t]), places=9)


@unittest.skipUnless(HAS, "needs torch and model.a.splicepwm")
class BalancedMissingColumns(unittest.TestCase):
    """stalin-0118 P2: a source that observed nothing in a column must leave
    the column's frequencies normalised, not shrink them."""

    def test_zero_pseudocount_column_stays_normalised(self):
        observed = _Ex("ACGTACGT", [(4, 6)])
        missing = _Ex("ACGTNACGT", [(4, 6)])   # the donor column is an N
        kw = dict(donor_span=(0, 1), acceptor_span=(0, 1), pseudocount=0.0,
                  min_sites=1)
        one = splicepwm.fit_grouped([("observed", [observed])], **kw)
        two = splicepwm.fit_grouped([("observed", [observed]),
                                     ("missing", [missing])], **kw)
        self.assertEqual(one.background, two.background)

        def mass(pwm):
            return sum(math.exp(v) * q for v, q in zip(pwm.donor[0], pwm.background))

        self.assertAlmostEqual(mass(one), 1.0)
        self.assertAlmostEqual(mass(two), 1.0)
        # and the surviving source decides the column on its own
        for a, b in zip(one.donor[0], two.donor[0]):
            self.assertAlmostEqual(a, b)

    def test_partially_missing_column_averages_over_contributors(self):
        # Two sources observe the column, a third does not: the average is
        # over the two, so adding an uninformative source cannot penalise
        # every base of the column.
        a = _Ex("ACGTACGT", [(4, 6)])
        b = _Ex("ACGTCCGT", [(4, 6)])
        blind = _Ex("ACGTNACGT", [(4, 6)])
        kw = dict(donor_span=(0, 1), acceptor_span=(0, 1), pseudocount=0.0,
                  min_sites=1)
        two = splicepwm.fit_grouped([("a", [a]), ("b", [b])], **kw)
        three = splicepwm.fit_grouped([("a", [a]), ("b", [b]),
                                       ("blind", [blind])], **kw)
        # backgrounds differ (the third source has its own composition), so
        # compare frequencies, which the log-odds are taken of
        for col_two, col_three in ((two.donor[0], three.donor[0]),):
            f2 = [math.exp(v) * q for v, q in zip(col_two, two.background)]
            f3 = [math.exp(v) * q for v, q in zip(col_three, three.background)]
            self.assertAlmostEqual(sum(f2), 1.0)
            self.assertAlmostEqual(sum(f3), 1.0)
            for x, y in zip(f2, f3):
                self.assertAlmostEqual(x, y)

    def test_a_column_no_source_observed_is_flat(self):
        blind = _Ex("ACGTNNGT", [(4, 6)])
        pwm = splicepwm.fit_grouped([("blind", [blind])], donor_span=(0, 1),
                                    acceptor_span=(0, 1), pseudocount=0.0,
                                    min_sites=1)
        self.assertEqual(pwm.donor[0], (0.0,) * 4)

    def test_default_pseudocount_is_unaffected(self):
        # Every source contributes to every column when smoothed, so the
        # committed matrices' estimator is the same before and after the fix.
        a = _Ex("ACGTACGT", [(4, 6)])
        blind = _Ex("ACGTNACGT", [(4, 6)])
        pwm = splicepwm.fit_grouped([("a", [a]), ("blind", [blind])],
                                    donor_span=(0, 1), acceptor_span=(0, 1),
                                    min_sites=1)
        f = [math.exp(v) * q for v, q in zip(pwm.donor[0], pwm.background)]
        self.assertAlmostEqual(sum(f), 1.0)


@unittest.skipUnless(HAS, "needs torch and model.a.splicepwm")
class Scale(unittest.TestCase):
    """The section 3.7 scalar on the matrix."""

    def test_scale_multiplies_every_score(self):
        x = "ACGTGTAAGTCCCCCCCCCCCCCCCCAGGTACGT"
        pwm = _uniform_pwm(donor_span=(-3, 6), acceptor_span=(-20, 3))
        one = splicepwm.bias(x, pwm)
        half = splicepwm.bias(x, pwm, scale=0.5)
        self.assertTrue(bool(torch.allclose(half, one * 0.5)))
        loud = splicepwm.bias(x, pwm, scale=2.5)
        self.assertTrue(bool(torch.allclose(loud, one * 2.5)))

    def test_scale_zero_is_exactly_no_pwm(self):
        x = "ACGTGTAAGTCCCCCCCCCCCCCCCCAGGTACGT"
        pwm = _uniform_pwm(donor_span=(-3, 6), acceptor_span=(-20, 3))
        z = splicepwm.bias(x, pwm, scale=0.0)
        self.assertTrue(bool(torch.all(z == 0)))

    def test_scale_zero_survives_an_infinite_entry(self):
        # pseudocount 0 puts -inf on an unseen base; 0 * -inf is nan, which
        # would poison the whole decode, so scale=0 must short-circuit.
        ex = [_Ex("AAAAGTAAGCCCCCCAG", [(4, 15)])]
        pwm = splicepwm.fit(ex, donor_span=(0, 2), acceptor_span=(-2, 0),
                            pseudocount=0.0, min_sites=1)
        self.assertTrue(any(math.isinf(v) for col in pwm.donor for v in col))
        z = splicepwm.bias("AAAAGTAAGCCCCCCAG", pwm, scale=0.0)
        self.assertTrue(bool(torch.all(z == 0)))

    def test_negative_scale_is_refused(self):
        pwm = _uniform_pwm()
        with self.assertRaises(ValueError):
            splicepwm.bias("ACGTACGT", pwm, scale=-1.0)
        with self.assertRaises(ValueError):
            pwm.tables(scale=float("nan"))


@unittest.skipUnless(HAS, "needs torch and model.a.splicepwm")
class Serialisation(unittest.TestCase):
    def test_round_trip(self):
        rng = random.Random(1)
        exs = [_Ex("".join(rng.choice("ACGT") for _ in range(300)), [(80, 200)])
               for _ in range(50)]
        pwm = splicepwm.fit(exs, sources=[{"name": "x", "dev_seqids": ["NC_1"]}],
                            min_sites=1)
        with tempfile.TemporaryDirectory() as d:
            path = str(Path(d) / "pwm.json")
            pwm.dump(path)
            back = splicepwm.SplicePWM.load(path)
        self.assertEqual(back, pwm)
        self.assertTrue(bool(torch.allclose(splicepwm.bias("ACGTACGTACGT" * 4, back),
                                            splicepwm.bias("ACGTACGTACGT" * 4, pwm))))

    def test_shape_is_validated(self):
        with self.assertRaises(ValueError):
            splicepwm.SplicePWM(donor=((0.0,) * 4,), acceptor=((0.0,) * 4,),
                                background=(0.25,) * 4,
                                donor_span=(-3, 6), acceptor_span=(0, 1))
        with self.assertRaises(ValueError):
            splicepwm.SplicePWM(donor=((0.0,) * 3,), acceptor=((0.0,) * 4,),
                                background=(0.25,) * 4,
                                donor_span=(0, 1), acceptor_span=(0, 1))
        with self.assertRaises(ValueError):
            splicepwm.SplicePWM(donor=((0.0,) * 4,), acceptor=((0.0,) * 4,),
                                background=(0.5,) * 2,
                                donor_span=(0, 1), acceptor_span=(0, 1))

    def test_base_order_is_checked(self):
        d = _uniform_pwm().to_json()
        d["bases"] = "ACGU"
        with self.assertRaises(ValueError):
            splicepwm.SplicePWM.from_json(d)


if __name__ == "__main__":
    unittest.main()
