#!/usr/bin/env python3
"""An octagonal cut-and-project patch and the four-bit projected SAT cube.

Only the lattice patch is window-filtered.  All 16 SAT assignments are retained.
The physical projection does not define the SAT adjacency: a SAT move flips one
bit, regardless of how close two displayed physical coordinates happen to be.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path


R = 1 / math.sqrt(2)
PHYSICAL_COLUMNS = ((1.0, 0.0), (R, R), (0.0, 1.0), (-R, R))
INTERNAL_COLUMNS = ((1.0, 0.0), (-R, R), (0.0, -1.0), (R, R))
WINDOW_SHIFT = (0.071, 0.113)
PATCH_RADIUS = 4.5
ENUMERATION_BOUND = 5
CANDIDATE_PADDING = 0.5
PHASON_X_RANGE = (-0.35, 0.35)


def project(vector, columns):
    return tuple(sum(value * column[d] for value, column in zip(vector, columns)) for d in range(2))


def cross(origin, a, b):
    return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])


def convex_hull(points):
    points = sorted(set(points))
    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 1e-12:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 1e-12:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def polygon_area(polygon):
    return abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(polygon, polygon[1:] + polygon[:1]))) / 2


BASE_WINDOW = convex_hull(project(corner, INTERNAL_COLUMNS) for corner in itertools.product((-0.5, 0.5), repeat=4))
WINDOW = [(p[0] + WINDOW_SHIFT[0], p[1] + WINDOW_SHIFT[1]) for p in BASE_WINDOW]


def inside_window(point, tolerance=1e-10):
    return all(cross(a, b, point) >= -tolerance for a, b in zip(WINDOW, WINDOW[1:] + WINDOW[:1]))


def reconstruct(physical, internal):
    return tuple(0.5 * sum(physical[d] * PHYSICAL_COLUMNS[j][d] + internal[d] * INTERNAL_COLUMNS[j][d] for d in range(2)) for j in range(4))


def generate_geometry():
    candidates = []
    candidate_bounds = [
        [min(p[d] for p in WINDOW) - CANDIDATE_PADDING, max(p[d] for p in WINDOW) + CANDIDATE_PADDING]
        for d in range(2)
    ]
    for lattice in itertools.product(range(-ENUMERATION_BOUND, ENUMERATION_BOUND + 1), repeat=4):
        physical = project(lattice, PHYSICAL_COLUMNS)
        if math.hypot(*physical) > PATCH_RADIUS + 1e-10:
            continue
        internal = project(lattice, INTERNAL_COLUMNS)
        if all(low <= internal[d] <= high for d, (low, high) in enumerate(candidate_bounds)):
            candidates.append({"lattice": list(lattice), "physical": list(physical), "internal": list(internal)})
    candidates.sort(key=lambda point: (round(math.hypot(*point["physical"]), 10), math.atan2(point["physical"][1], point["physical"][0]), tuple(point["lattice"])))

    def index_with_edges(records):
        lookup = {tuple(point["lattice"]): i for i, point in enumerate(records)}
        edges = []
        for i, point in enumerate(records):
            point["id"] = i
            for direction in range(4):
                neighbor = point["lattice"].copy()
                neighbor[direction] += 1
                if tuple(neighbor) in lookup:
                    edges.append([i, lookup[tuple(neighbor)]])
        return edges

    candidate_edges = index_with_edges(candidates)
    points = [dict(point) for point in candidates if inside_window(point["internal"])]
    edges = index_with_edges(points)

    boolean_vertices = []
    for assignment_id in range(16):
        bits = [(assignment_id >> j) & 1 for j in range(4)]
        physical = project(bits, PHYSICAL_COLUMNS)
        internal = project(bits, INTERNAL_COLUMNS)
        boolean_vertices.append({
            "id": assignment_id,
            "bits": bits,
            "label": "".join(map(str, bits)),
            "physical": list(physical),
            "internal": list(internal),
            "insideWindow": inside_window(internal),
        })
    boolean_edges = [[i, i ^ (1 << j)] for i in range(16) for j in range(4) if i < (i ^ (1 << j))]
    maximum_internal_norm = max(math.hypot(*p) for p in BASE_WINDOW) + math.hypot(*WINDOW_SHIFT)
    coordinate_bound = (PATCH_RADIUS + maximum_internal_norm) / 2
    maximum_candidate_internal_norm = math.hypot(*(max(abs(low), abs(high)) for low, high in candidate_bounds))
    candidate_coordinate_bound = (PATCH_RADIUS + maximum_candidate_internal_norm) / 2

    result = {
        "metadata": {
            "description": "Finite disk crop of a regular octagonal model set from a rank-four lattice; independently, the unfiltered four-bit cube under the same projection.",
            "physicalColumns": [list(p) for p in PHYSICAL_COLUMNS],
            "internalColumns": [list(p) for p in INTERNAL_COLUMNS],
            "windowDefinition": "P_internal([-1/2,1/2]^4) + windowShift",
            "windowShift": list(WINDOW_SHIFT),
            "patchRadius": PATCH_RADIUS,
            "enumerationBound": ENUMERATION_BOUND,
            "provenAbsoluteLatticeCoordinateBound": coordinate_bound,
            "provenAbsoluteCandidateLatticeCoordinateBound": candidate_coordinate_bound,
            "candidateInternalBounds": candidate_bounds,
            "candidatePadding": CANDIDATE_PADDING,
            "phasonXRange": list(PHASON_X_RANGE),
            "phasonConvention": "Translate the centered baseWindow by (windowShift[0] + phasonX, windowShift[1]); retain a candidate when its internal point lies inside.",
            "windowArea": polygon_area(BASE_WINDOW),
            "projectedLatticeCovolume": 4.0,
            "modelSetDensity": polygon_area(BASE_WINDOW) / 4,
            "bitConvention": "bits[j] = (id >> j) & 1; label lists x0,x1,x2,x3 left to right, least significant bit first",
            "satAdjacency": "One Boolean bit flip; booleanEdges records the 32 cube edges.",
            "warning": "Window membership is illustrative for SAT vertices. No SAT assignments are filtered out.",
        },
        "baseWindow": [list(p) for p in BASE_WINDOW],
        "window": [list(p) for p in WINDOW],
        "patch": {"points": points, "edges": edges},
        "candidates": {"points": candidates, "edges": candidate_edges},
        "booleanVertices": boolean_vertices,
        "booleanEdges": boolean_edges,
    }
    validate(result)
    return result


def validate(geometry):
    assert len(BASE_WINDOW) == 8
    assert math.isclose(polygon_area(BASE_WINDOW), 2 * (1 + math.sqrt(2)), abs_tol=1e-10)
    assert ENUMERATION_BOUND > geometry["metadata"]["provenAbsoluteLatticeCoordinateBound"]
    assert ENUMERATION_BOUND > geometry["metadata"]["provenAbsoluteCandidateLatticeCoordinateBound"]
    assert len(geometry["candidates"]["points"]) < 1000
    assert len(geometry["booleanVertices"]) == 16
    assert len(geometry["booleanEdges"]) == 32
    assert len({tuple(round(v, 10) for v in p["physical"]) for p in geometry["booleanVertices"]}) == 16
    for point in geometry["patch"]["points"]:
        assert inside_window(point["internal"])
        assert math.hypot(*point["physical"]) <= PATCH_RADIUS + 1e-10
        assert max(abs(a - b) for a, b in zip(reconstruct(point["physical"], point["internal"]), point["lattice"])) < 1e-10
    for point in geometry["candidates"]["points"]:
        assert math.hypot(*point["physical"]) <= PATCH_RADIUS + 1e-10
        assert all(low <= point["internal"][d] <= high for d, (low, high) in enumerate(geometry["metadata"]["candidateInternalBounds"]))
        assert max(abs(a - b) for a, b in zip(reconstruct(point["physical"], point["internal"]), point["lattice"])) < 1e-10
    for point in geometry["booleanVertices"]:
        assert max(abs(a - b) for a, b in zip(reconstruct(point["physical"], point["internal"]), point["bits"])) < 1e-10
    for pointset in (geometry["patch"], geometry["candidates"]):
        for i, j in pointset["edges"]:
            a, b = pointset["points"][i], pointset["points"][j]
            assert sum(abs(x - y) for x, y in zip(a["lattice"], b["lattice"])) == 1
            assert math.isclose(math.dist(a["physical"], b["physical"]), 1, abs_tol=1e-10)
    for shift in PHASON_X_RANGE:
        for corner in BASE_WINDOW:
            shifted = (corner[0] + WINDOW_SHIFT[0] + shift, corner[1] + WINDOW_SHIFT[1])
            assert all(low <= shifted[d] <= high for d, (low, high) in enumerate(geometry["metadata"]["candidateInternalBounds"]))
    for i, j in geometry["booleanEdges"]:
        assert (i ^ j).bit_count() == 1


if __name__ == "__main__":
    geometry = generate_geometry()
    destination = Path(__file__).with_name("jewel_geometry.json")
    destination.write_text(json.dumps(geometry, indent=2) + "\n")
    print(json.dumps({
        "file": str(destination),
        "patch_points": len(geometry["patch"]["points"]),
        "patch_edges": len(geometry["patch"]["edges"]),
        "candidate_points": len(geometry["candidates"]["points"]),
        "candidate_edges": len(geometry["candidates"]["edges"]),
        "boolean_vertices": len(geometry["booleanVertices"]),
        "boolean_edges": len(geometry["booleanEdges"]),
        "max_absolute_lattice_coordinate_bound": geometry["metadata"]["provenAbsoluteLatticeCoordinateBound"],
        "checks": "window, coordinate reconstruction, octagonal area, unit lattice edges, injective Boolean projection, enumeration completeness bound",
    }))
