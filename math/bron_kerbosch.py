def get_maximal_cliques(neighbors: dict):
    """
    Get all maximal cliques in an undirected graph
    Args:
        neighbors: dict mapping nodes to sets of neighbors

    Returns:
        Generator of frozensets for all maximal cliques
    """
    return _BK(set(neighbors), neighbors)


def _BK(P, neighbors, R=frozenset(), X=frozenset()):
    if not P and not X:
        yield R
    else:
        u = max(P | X, key=lambda p: len(P & neighbors[p]))
        for v in P - neighbors[u]:
            yield from _BK(P & neighbors[v], neighbors, R=R | {v}, X=X & neighbors[v])
            P = P - {v}
            X = X | {v}
