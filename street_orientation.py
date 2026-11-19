# -*- coding: utf-8 -*-
"""Orientacion de conductos de calle segun el DEM, coherente a nivel de red.

Sin dependencias de QGIS (se puede probar con Python plano).

Metodo: "priority-flood" sobre el grafo de calles, sembrado desde los
OUTFALLS (y, en componentes sin outfall, desde su nodo de menor cota).
Cada nodo recibe una cota rellenada z' = max(z, z' del nodo desde el que se
alcanzo) y una distancia de flujo al outlet. Los ruidos locales del DEM
(pozos/lomos) y los tramos planos quedan con z' igual y se resuelven por
distancia al outlet, no por una comparacion aislada de dos cotas.

Orden total de nodos: (z' cuantizada, distancia, id). Cada tramo se orienta
del nodo mayor al menor => red sin ciclos, todo nodo drena a un outlet, y
los OUTFALLS quedan siempre aguas abajo de su unico conducto.
"""
import heapq
from collections import defaultdict

FLAT_DZ_M = 0.10     # diferencias de cota rellenada menores a esto cuentan como plano


def orient_edges(edges, node_z, outfalls, flat_dz=FLAT_DZ_M):
    """edges: [(a, b, largo_m)]; node_z: {nodo: cota}; outfalls: iterable de nodos.

    Devuelve (flip, stats): flip[i] = True si el tramo i debe invertirse
    (de a->b pasa a b->a) para que quede de mayor a menor cota global."""
    adj = defaultdict(list)
    for a, b, length in edges:
        if a == b:
            continue
        adj[a].append((b, max(float(length), 1e-6)))
        adj[b].append((a, max(float(length), 1e-6)))

    outfalls = {n for n in outfalls if n in adj}
    zf, dist, seen = {}, {}, set()
    heap = []

    def push_seed(n):
        zf[n] = node_z[n]
        dist[n] = 0.0
        heapq.heappush(heap, (zf[n], 0.0, n))

    for n in sorted(outfalls):
        push_seed(n)

    n_comp_seeds = 0
    nodes_sorted = sorted(adj)

    def flood():
        while heap:
            z, d, n = heapq.heappop(heap)
            if n in seen:
                continue
            seen.add(n)
            for m, length in adj[n]:
                if m in zf:
                    continue
                zf[m] = max(node_z[m], z)
                dist[m] = d + length
                heapq.heappush(heap, (zf[m], dist[m], m))

    flood()
    # componentes sin outfall: se siembra el nodo de menor cota (grado 1 si hay)
    for start in nodes_sorted:
        if start in zf:
            continue
        comp, stack = [], [start]
        mark = {start}
        while stack:
            u = stack.pop()
            comp.append(u)
            for v, _l in adj[u]:
                if v not in mark and v not in zf:
                    mark.add(v)
                    stack.append(v)
        ends = [u for u in comp if len(adj[u]) == 1] or comp
        push_seed(min(ends, key=lambda u: (node_z[u], u)))
        n_comp_seeds += 1
        flood()

    step = max(float(flat_dz), 1e-9)

    def key(n):
        return (int(zf[n] // step), dist[n], n)

    flip, n_flip, n_flat, n_adverse = [], 0, 0, 0
    for a, b, _l in edges:
        if a == b:
            flip.append(False)
            continue
        f = key(a) < key(b)          # a es el menor => hoy fluye hacia arriba
        flip.append(f)
        n_flip += int(f)
        if abs(zf[a] - zf[b]) < step:
            n_flat += 1
        lo, hi = (a, b) if f else (b, a)   # hi->lo es el sentido final (from=hi, to=lo)
        if node_z[hi] < node_z[lo] - step:
            n_adverse += 1
    stats = {"flipped": n_flip, "flat": n_flat, "adverse_left": n_adverse,
             "comp_seeds": n_comp_seeds, "edges": len(edges)}
    return flip, stats
