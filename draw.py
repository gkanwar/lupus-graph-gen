import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

class Edge:
    def __init__(self, i, j, particle):
        self.i = i
        self.j = j
        self.particle = particle

class Node:
    edges: list[Edge]
    tag: str | None
    x: np.ndarray | None
    fixed: bool
    def __init__(self):
        self.edges = []
        self.tag = None
        self.x = np.array([0.0, 0.0])
        self.fixed = False

def make_nodes(n):
    return [Node() for _ in range(n)]

class Particle:
    charged: bool
    style: str # 'wavy', 'dashed', 'solid'
    weight: float # heavier lines affect layout more
    def __init__(self, charged, style, weight):
        self.charged = charged
        self.style = style
        self.weight = weight

class Graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges

def connect(n1, n2, particle):
    edge = Edge(n1, n2, particle)
    n1.edges.append(edge)
    n2.edges.append(edge)
    return edge

def make_clean_fig():
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    # ax.set_axis_off()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_aspect(1.0)
    ax.grid()
    for tic in ax.xaxis.get_major_ticks():
        tic.tick1line.set_visible(False)
        tic.tick2line.set_visible(False)
    for tic in ax.yaxis.get_major_ticks():
        tic.tick1line.set_visible(False)
        tic.tick2line.set_visible(False)
    return fig, ax

def make_arc_positions(n, radius, *, min_theta=-np.pi/2, max_theta=np.pi/2):
    """Make n evenly spaced points on an arc"""
    th = np.linspace(min_theta, max_theta, endpoint=True, num=n+2)
    pts = np.stack([np.cos(th)*radius, np.sin(th)*radius], axis=-1)
    return pts[1:-1]

def draw_node(pt, *, ax):
    ax.plot(*pt, marker='o', color='k')

def draw_edge(x1, x2, particle, *, ax):
    style = dict(marker='', color='k')
    rc_context = dict()
    if particle.style == 'solid':
        style['linestyle'] = '-'
    elif particle.style == 'dashed':
        style['linestlye'] = '--'
    elif particle.style == 'wavy':
        style['linestyle'] = '-'
        rc_context['path.sketch'] = (5, 15, 1)
    else:
        raise ValueError(style['linestyle'])
    with mpl.rc_context(rc_context):
        ax.plot(*np.transpose([x1, x2]), **style)

def apply_spring_layout(graph):
    """Distribute nodes according to spring layout, preserving fixed nodes"""
    inds = [i for i,node in enumerate(graph.nodes) if not node.fixed]
    # minimize the energy sum_{xy} weight*(x-y)^2/2
    # for y fixed, this reduces to weight*x^2/2 - weight*x*y
    # this can be cast as a simple Gaussian minimization problem:
    # X* = argmin_X X.W.X^T - X.v^T - v.X^T
    # in component notation
    # {xi}* = argmin_{xi} xi Wij xj - 2 xi vi
    # differentiating and solving,
    # 0 = Wij xj + xj Wji - 2 vi
    # vi = Wij xj
    # ===> X* = W^-1 v
    W = np.zeros((len(inds), len(inds)), dtype=np.float64)
    v = np.zeros((len(inds),2), dtype=np.float64)
    for i,ni in enumerate(inds):
        node = graph.nodes[ni]
        for edge in node.edges:
            node2 = edge.i if edge.j is node else edge.j
            if node2.fixed:
                # print(f'adding {node2.x=} to {i=}')
                v[i] += edge.particle.weight*node2.x
                W[i,i] += edge.particle.weight
            else:
                # print(f'found unfixed edge')
                j = inds.index(graph.nodes.index(node2))
                W[i,i] += edge.particle.weight
                W[i,j] -= edge.particle.weight
                W[j,i] -= edge.particle.weight
                W[j,j] += edge.particle.weight
    assert np.allclose(np.transpose(W), W)
    # print(f'{W=}')
    # print(f'{v=}')
    x = np.linalg.pinv(W) @ v
    # CHECK:
    # E = np.einsum('ia,ij,ja', x, W, x) - 2*np.sum(x*v)
    # for _ in range(10):
    #     dx = np.random.normal(size=x.shape)*0.01
    #     xp = x + dx
    #     Ep = np.einsum('ia,ij,ja', xp, W, xp) - 2*np.sum(xp*v)
    #     print(f'{E} vs {Ep}')
    #     assert E <= Ep
    # print(f'final assingment {x=}')
    for pt,ni in zip(x, inds):
        assert not graph.nodes[ni].fixed
        graph.nodes[ni].x = pt

def layout(graph):
    """Automatic layout into PDF figure"""
    fig, ax = make_clean_fig()
    in_nodes = list(filter(lambda n: n.tag == 'in', graph.nodes))
    out_nodes = list(filter(lambda n: n.tag == 'out', graph.nodes))
    other_nodes = list(filter(lambda n: n.tag is None, graph.nodes))
    assert len(in_nodes) + len(out_nodes) + len(other_nodes) == len(graph.nodes)
    # stage 1: position in + out nodes on a circle
    RAD = 5.0
    in_pos = make_arc_positions(len(in_nodes), RAD)
    in_pos[:,0] *= -1
    out_pos = make_arc_positions(len(out_nodes), RAD)
    for pt,node in zip(in_pos, in_nodes):
        node.x = pt
        node.fixed = True
    for pt,node in zip(out_pos, out_nodes):
        node.x = pt
        node.fixed = True
    # stage 2: spring layout
    apply_spring_layout(graph)
    # stage 3: bend propagators?
    # stage 4: render propagators and nodes
    for edge in graph.edges:
        draw_edge(edge.i.x, edge.j.x, edge.particle, ax=ax)
    for node in graph.nodes:
        draw_node(node.x, ax=ax)
    return fig

def charged_particle():
    return Particle(True, 'solid', 5.0)
def photon():
    return Particle(False, 'wavy', 1.0)

def make_qed_t_channel():
    e = charged_particle()
    g = photon()
    nodes = make_nodes(6)
    in1, in2, out1, out2, v1, v2 = nodes
    in1.tag = in2.tag = 'in'
    out1.tag = out2.tag = 'out'
    edges = []
    edges.append(connect(in1, v1, e))
    edges.append(connect(v1, out1, e))
    edges.append(connect(v1, v2, g))
    edges.append(connect(in2, v2, e))
    edges.append(connect(v2, out2, e))
    return Graph(nodes, edges)

def main():
    fig = layout(make_qed_t_channel())
    plt.show()

if __name__ == '__main__':
    main()
