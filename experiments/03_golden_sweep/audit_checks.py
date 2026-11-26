#!/usr/bin/env python3
"""Independent mathematical and result audit for the jewel sweep.

Does not import the search implementation or use planted answers to solve.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import random
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def violations(state, clauses):
    return sum(not any(bool(state & (1 << (abs(lit) - 1))) == (lit > 0)
                       for lit in clause) for clause in clauses)


def exact_bitset(n, clauses):
    """Exhaustively evaluate every assignment using Python integer bitsets."""
    count = 1 << n
    all_assignments = (1 << count) - 1
    masks = []
    for j in range(n):
        block = 1 << j
        repeat = all_assignments // ((1 << (2 * block)) - 1)
        masks.append(repeat * (((1 << block) - 1) << block))
    surviving = all_assignments
    for clause in clauses:
        clause_mask = 0
        for lit in clause:
            mask = masks[abs(lit) - 1]
            clause_mask |= mask if lit > 0 else all_assignments ^ mask
        surviving &= clause_mask
    witness = (surviving & -surviving).bit_length() - 1 if surviving else None
    return {"label": "SAT" if surviving else "UNSAT", "witness": witness,
            "solution_count": surviving.bit_count(), "assignments_checked": count}


def exact_dpll(n, clauses, max_nodes=200000):
    """Independent complete branching with explicit UNKNOWN on node limit."""
    encoded = []
    for clause in clauses:
        positive = negative = 0
        for lit in clause:
            if lit > 0:
                positive |= 1 << (lit - 1)
            else:
                negative |= 1 << (-lit - 1)
        encoded.append((positive, negative))
    nodes = 0
    universe = (1 << n) - 1

    def visit(assigned, true):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            return None, None
        while True:
            false = assigned & ~true
            residual = []
            positive_units = negative_units = 0
            for positive, negative in encoded:
                if positive & true or negative & false:
                    continue
                unassigned = (positive | negative) & ~assigned
                if not unassigned:
                    return False, None
                if unassigned & (unassigned - 1) == 0:
                    if positive & unassigned:
                        positive_units |= unassigned
                    else:
                        negative_units |= unassigned
                else:
                    residual.append((positive & unassigned, negative & unassigned))
            if positive_units & negative_units:
                return False, None
            units = positive_units | negative_units
            if not units:
                break
            assigned |= units
            true |= positive_units
        if not residual:
            return True, true & universe
        scores = defaultdict(lambda: [0.0, 0.0])
        for positive, negative in residual:
            weight = 2.0 ** -(positive | negative).bit_count()
            for polarity, bits in enumerate((negative, positive)):
                while bits:
                    bit = bits & -bits
                    scores[bit][polarity] += weight
                    bits ^= bit
        bit = max(scores, key=lambda b: sum(scores[b]))
        prefer_true = scores[bit][1] >= scores[bit][0]
        unknown = False
        for value in (prefer_true, not prefer_true):
            status, witness = visit(assigned | bit, true | bit if value else true & ~bit)
            if status is True:
                return status, witness
            unknown |= status is None
        return (None if unknown else False), None

    label, witness = visit(0, 0)
    return {"label": "SAT" if label is True else "UNSAT" if label is False else "UNKNOWN",
            "witness": witness, "nodes": nodes, "node_cap": max_nodes}


def transform(clauses, n, permutation, flips, state, weights, order):
    """new y_i = old x_permutation[i] XOR flips_i; reorder clauses too."""
    inverse = {old: new for new, old in enumerate(permutation)}
    changed = []
    for clause in clauses:
        new_clause = []
        for lit in clause:
            i = inverse[abs(lit) - 1]
            sign = (lit > 0) ^ bool(flips & (1 << i))
            new_clause.append((i + 1) * (1 if sign else -1))
        changed.append(new_clause)
    new_state = sum((((state >> old) & 1) ^ ((flips >> i) & 1)) << i
                    for i, old in enumerate(permutation))
    return [changed[i] for i in order], new_state, [weights[i] for i in order]


def energy(state, clauses, weights):
    return sum(w for clause, w in zip(clauses, weights) if violations(state, [clause]))


def invariant_checks():
    rng = random.Random(7162026)
    phi = (1 + math.sqrt(5)) / 2
    alias_error = reversed_error = 0.0
    for j in range(513):
        a, b, c = j / phi, j / (1 / phi), j / phi ** 2
        alias_error = max(alias_error, abs(math.sin(math.pi * (a - b))))
        reversed_error = max(reversed_error, abs(math.sin(math.pi * (a + c))))
    assert alias_error < 1e-10 and reversed_error < 1e-10
    positive_checks = transformations = dpll_comparisons = 0
    for trial in range(32):
        n = 4 + trial % 7
        clauses = [[(i + 1) * rng.choice((-1, 1)) for i in rng.sample(range(n), 3)]
                   for _ in range(round(4.2 * n))]
        weights = [0.2 + 2 * rng.random() for _ in clauses]
        permutation = rng.sample(range(n), n)
        flips = rng.randrange(1 << n)
        order = rng.sample(range(len(clauses)), len(clauses))
        for state in range(1 << min(n, 6)):
            e = energy(state, clauses, weights)
            assert (e == 0) == (violations(state, clauses) == 0)
            positive_checks += 1
            new_clauses, new_state, new_weights = transform(
                clauses, n, permutation, flips, state, weights, order)
            assert violations(state, clauses) == violations(new_state, new_clauses)
            assert math.isclose(e, energy(new_state, new_clauses, new_weights), abs_tol=1e-10)
            # Every original one-bit edge maps to a new one-bit edge.
            for variable in range(n):
                _, neighbor, _ = transform(clauses, n, permutation, flips,
                                          state ^ (1 << variable), weights, order)
                assert (neighbor ^ new_state).bit_count() == 1
            transformations += 1
        brute = exact_bitset(n, clauses)
        independent = exact_dpll(n, clauses)
        assert brute["label"] == independent["label"]
        if independent["label"] == "SAT":
            assert violations(independent["witness"], clauses) == 0
        dpll_comparisons += 1
    strict_greedy_obstructions = 0
    for trial in range(65):
        fixed_weights = [1.0] * 7 if trial == 0 else [0.1 + 10 * rng.random() for _ in range(7)]
        index = min(range(7), key=fixed_weights.__getitem__)
        forbidden = list(range(7))
        forbidden[0], forbidden[index] = forbidden[index], forbidden[0]
        clauses = [[(j + 1) * (-1 if state & (1 << j) else 1) for j in range(3)]
                   for state in forbidden]
        assert [s for s in range(8) if violations(s, clauses) == 0] == [7]
        start_energy = energy(0, clauses, fixed_weights)
        assert start_energy > 0
        assert all(energy(1 << j, clauses, fixed_weights) >= start_energy - 1e-12 for j in range(3))
        strict_greedy_obstructions += 1
    return {"positive_zero_set_cases": positive_checks,
            "variable_permutation_gauge_and_clause_reorder_cases": transformations,
            "dpll_compared_with_exhaustive_instances": dpll_comparisons,
            "golden_frequency_alias_max_error": alias_error,
            "golden_frequency_reversal_max_error": reversed_error,
            "fixed_weight_strict_greedy_obstructions_checked": strict_greedy_obstructions,
            "obstruction_theorem": "For any seven positive weights fixed independently of clause contents, put the clause forbidding 000 at a minimum-weight index and use all clauses forbidding assignments other than 111. Then 111 is uniquely satisfying while strict one-bit descent stalls at 000. This obstruction does not apply to adaptive weights, sideways moves, nonlocal moves, or content-dependent choices.",
            "golden_identity": "1/phi = phi-1; 1/phi^2 = 1-1/phi. Thus j/phi and j/(1/phi) agree modulo integers; j/phi^2 reverses that rotation. Tile phases must also be reflected/rephased before claiming full schedule equivalence."}


def main():
    started = time.perf_counter()
    report = {"invariants": invariant_checks(), "mismatches": [], "limitations": [
        "Finite successes and favorable scaling plots do not prove a worst-case polynomial bound.",
        "Heuristic timeout or trap is UNKNOWN, never an UNSAT certificate.",
        "Runs sharing a formula are correlated; configuration selection introduces multiple-comparison bias.",
        "A held-out comparison must keep the training-selected configuration fixed.",
        "Golden-ratio reciprocal/power schedules can be rotation aliases rather than independent constants.",
        "A tiling layout does not ensure SAT locality; clause-to-tile permutations test representation sensitivity."]}
    dataset_file, runs_file = ROOT / "dataset.json", ROOT / "raw_runs.csv"
    if dataset_file.exists() and runs_file.exists():
        data = json.loads(dataset_file.read_text())
        instances = {str(i["id"]): i for i in data["instances"]}
        configurations = {str(c["id"]): c for c in data["configurations"]}
        geometry_rows = 0
        schedule_trials = 0
        rng = random.Random(984216)
        lengths = sorted({instance["m"] for instance in instances.values()})
        for tiling_name, tiling in data["tilings"].items():
            assert len(tiling["angles"]) >= max(lengths)
            assert len(tiling["neighbors"]) == len(tiling["angles"])
            assert all(math.isfinite(a) for a in tiling["angles"])
            for row in tiling["neighbors"]:
                assert all(0 <= neighbor < len(tiling["angles"]) for neighbor in row)
            for m in lengths:
                weights = [0.1 + rng.random() * 10 for _ in range(m)]
                for j in range(m):
                    retained = [v for v in tiling["neighbors"][j] if v < m and v != j]
                    denominator = len(retained) + 1
                    coefficients = [0.8 + 0.2 / denominator] + [0.2 / denominator] * len(retained)
                    assert all(v >= 0 for v in coefficients)
                    assert math.isclose(sum(coefficients), 1, abs_tol=1e-12)
                    result = 0.8 * weights[j] + 0.2 * (weights[j] + sum(weights[k] for k in retained)) / denominator
                    assert min(weights) <= result <= max(weights)
                    geometry_rows += 1
        for config in configurations.values():
            if config["method"] not in ("static", "cos", "341", "diffusion"):
                continue
            tiling = data["tilings"][config["tiling"]]
            assert config["constant"] > 0
            for update in (0, 1, 37, 99):
                raw_boosts = []
                for j in range(max(lengths)):
                    angle = tiling["angles"][j] + 2 * math.pi * (j + 1 - update) / config["constant"]
                    if config["method"] in ("static", "cos"):
                        value = 1 + config["amplitude"] * math.cos(angle)
                    else:
                        value = 0.25 + 2 * (1 + math.cos(angle)) ** 2
                    assert value > 0 and math.isfinite(value)
                    raw_boosts.append(value)
                    schedule_trials += 1
                mean = sum(raw_boosts) / len(raw_boosts)
                assert math.isclose(sum(v / mean for v in raw_boosts) / len(raw_boosts), 1, abs_tol=1e-12)
        report["geometry_and_schedule_checks"] = {
            "truncated_diffusion_rows_checked": geometry_rows,
            "positive_schedule_values_checked": schedule_trials,
            "diffusion": "Each retained-neighbor+self averaging row is nonnegative and sums to one; .8I+.2P therefore preserves strictly positive weights.",
            "boost_control": "Mean raw boost is normalized to one over violated clauses; this removes average increment-size confounding at each trap."}
        report["limitations"].append("Configuration IDs seed different tie-breaking RNG streams; rotation-alias configurations can therefore have different sampled search paths without constituting distinct rotation mechanisms.")
        exact = {}
        for key, instance in instances.items():
            result = (exact_bitset(instance["n"], instance["clauses"])
                      if instance["n"] <= 16 else exact_dpll(instance["n"], instance["clauses"]))
            exact[key] = result
            expected = instance.get("sat_label", "UNKNOWN")
            if result["label"] != "UNKNOWN" and expected != "UNKNOWN" and result["label"] != expected:
                report["mismatches"].append({"instance": key, "expected": expected, "independent": result["label"]})
            if result["label"] == "SAT":
                assert violations(result["witness"], instance["clauses"]) == 0
            if instance.get("witness") is not None:
                assert violations(instance["witness"], instance["clauses"]) == 0
        report["independent_exact_labels"] = exact
        report["label_counts"] = dict(Counter(x["label"] for x in exact.values()))
        rows = list(csv.DictReader(runs_file.open()))
        checked_witnesses = 0
        scores = defaultdict(lambda: [0, 0, 0])
        for row in rows:
            instance = instances[str(row["instance"])]
            state = int(row["final_state"])
            assert 0 <= state < 1 << instance["n"], row
            actual = violations(state, instance["clauses"])
            success = str(row["success"]).lower() in ("1", "true")
            assert actual == int(row["final_violations"]), row
            assert success == (actual == 0), row
            checked_witnesses += success
            assert int(row["evaluations"]) <= int(data.get("budget", data.get("budget_assignment_evaluations"))), row
            if instance["split"] == "train" and exact[str(row["instance"])]["label"] == "SAT":
                score = scores[str(row["config"])]
                score[0] += success; score[1] += 1; score[2] += int(row["evaluations"])
        report["run_records_checked"] = len(rows)
        report["SAT_witnesses_checked"] = checked_witnesses
        experiments = [key for key in scores if configurations[key]["method"] not in
                       ("greedy", "restarts", "breakout", "walksat")]
        if experiments:
            best = max(experiments, key=lambda key: (scores[key][0] / scores[key][1],
                                                      -scores[key][2] / scores[key][1], key))
            failures = [row for row in rows if str(row["config"]) == best
                        and instances[str(row["instance"])]["split"] == "test"
                        and exact[str(row["instance"])]["label"] == "SAT"
                        and str(row["success"]).lower() not in ("1", "true")]
            report["independently_selected_training_best"] = configurations[best]
            report["selected_training_scores"] = scores[best]
            report["selected_heldout_SAT_failures"] = len(failures)
            if failures:
                row = failures[0]
                instance = instances[str(row["instance"])]; truth = exact[str(row["instance"])]
                report["concrete_heldout_counterexample"] = {
                    "configuration": configurations[best], "instance": instance, "failed_run": row,
                    "independent_satisfying_witness": truth["witness"],
                    "interpretation": "A concrete failure under the stated finite search cap, not a proof that the heuristic cannot eventually solve it."}
        stress_summary_file = ROOT / "permutation_stress_summary.json"
        stress_runs_file = ROOT / "permutation_stress_runs.csv.gz"
        if stress_summary_file.exists() and stress_runs_file.exists():
            stress = json.loads(stress_summary_file.read_text())
            maps = {(int(t["variant"]), str(t["instance"])): t for t in stress["transformations"]}
            transformed = {}
            rng = random.Random(628251)
            equiv_checks = 0
            for key, mapping in maps.items():
                original = instances[key[1]]
                permutation, gauge = mapping["variable_permutation"], mapping["gauge_flips"]
                order = mapping["clause_order"]
                assert sorted(permutation) == list(range(original["n"]))
                assert sorted(order) == list(range(original["m"]))
                assert set(gauge) <= {0, 1}
                clauses = [[(permutation[abs(lit) - 1] + 1) *
                            (1 if lit > 0 else -1) * (-1 if gauge[abs(lit) - 1] else 1)
                            for lit in original["clauses"][j]] for j in order]
                transformed[key] = clauses
                # Distinct convention from our generic invariant checker: here old j maps to new permutation[j].
                def forward(state):
                    return sum((bool(state & (1 << j)) ^ bool(gauge[j])) << permutation[j]
                               for j in range(original["n"]))
                for state in [exact[key[1]]["witness"]] + [rng.randrange(1 << original["n"]) for _ in range(32)]:
                    assert violations(state, original["clauses"]) == violations(forward(state), clauses)
                    equiv_checks += 1
            with gzip.open(stress_runs_file, "rt") as stream:
                stress_rows = list(csv.DictReader(stream))
            originals = {(row["instance"], row["config"], row["start_index"]): row for row in rows}
            aggregation = defaultdict(lambda: [0, 0])
            identities = 0
            for row in stress_rows:
                key = int(row["variant"]), row["instance"]
                state = int(row["final_state"])
                instance = instances[row["instance"]]
                assert 0 <= state < 1 << instance["n"]
                bad = violations(state, transformed[key])
                success = row["success"] in ("1", "true", "True")
                assert success == (bad == 0)
                assert bad == int(row["final_violations"])
                assert int(row["evaluations"]) <= int(data["budget_assignment_evaluations"])
                mapping = maps[key]
                old_state = sum((((state >> mapping["variable_permutation"][j]) & 1) ^ mapping["gauge_flips"][j]) << j
                                for j in range(instance["n"]))
                assert bad == violations(old_state, instance["clauses"])
                if int(row["variant"]) == 0:
                    old = originals[row["instance"], row["config"], row["start_index"]]
                    assert all(row[k] == v for k, v in old.items())
                    identities += 1
                else:
                    aggregation[row["config"]][0] += success
                    aggregation[row["config"]][1] += 1
            for item in stress["aggregate_permuted"]:
                assert aggregation[str(item["config"])] == [item["successes"], item["runs"]]
            report["permutation_stress_audit"] = {
                "transformations_verified": len(maps), "assignment_equivalences_checked": equiv_checks,
                "all_terminal_states_checked": len(stress_rows), "identity_replays_matching_exactly": identities,
                "new_permuted_runs": len(stress_rows) - identities,
                "counts_successes_runs": {configurations[k]["name"]: v for k, v in aggregation.items()},
                "interpretation": "Training winner remains frozen; starts, witnesses, variable labels, gauge bits, and clauses transform consistently. Reordered representations change bounded-heuristic outcomes, not SAT truth."}
    else:
        report["status"] = "Mathematical invariants audited; dataset/raw runs not yet present. Rerun after sweep."
    report["elapsed_seconds"] = time.perf_counter() - started
    (ROOT / "audit_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in
                     ("independent_exact_labels", "concrete_heldout_counterexample")}, indent=2))
    assert not report["mismatches"], report["mismatches"]


if __name__ == "__main__":
    main()
