"""
Agnieszka Szkutek, 208619
"""

import networkx as nx
import numpy as np
import random as rnd
import matplotlib.pyplot as plt
# import pylab
from timeit import default_timer as timer


# def save_graph(G, pathname, i):
#     fig = pylab.figure()
#     pos = nx.spring_layout(G)
#     size = 200
#
#     green = []
#     red = []
#     for node in G.nodes():
#         if G.node[node]['vote']:
#             green.append(node)
#         else:
#             red.append(node)
#
#     nx.draw_networkx_nodes(G, pos, nodelist=green, node_color='g', node_size=size, alpha=0.8)
#     nx.draw_networkx_nodes(G, pos, nodelist=red, node_color='r', node_size=size, alpha=0.8)
#
#     nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
#
#     s = pathname + "_" + str(i) + ".png"
#     pylab.savefig(s)
#     pylab.close(fig)


def new_graph(g, N):
    if g == 'complete':
        return nx.complete_graph(N)
    elif g == 'ba':
        return nx.barabasi_albert_graph(N, 4)
    elif g == 'ws1':
        return nx.watts_strogatz_graph(N, 4, 0.01)
    elif g == 'ws2':
        return nx.watts_strogatz_graph(N, 4, 0.2)


def MCsimulations(g, N, P, q=3, MC=100, MC2=1):
    magnP = []
    timeP = []
    # nodes = np.arange(0, N, 1)

    for p in P:
        final_magn = 0
        timeline0 = np.zeros(MC)
        for _ in range(MC2):
            timeline = one_timeline_simulation(g, N, p, q, MC)

            final_magn += timeline[-1]
            timeline0 = np.add(timeline0, timeline)

        magnP.append(final_magn / MC2)
        timeP.append(np.divide(timeline0, MC2))
    return magnP, timeP


def one_timeline_simulation(g, N, p, q, MC):
    G = new_graph(g, N)
    votes = np.ones(N, int)

    timeline = np.zeros(MC)
    timeline[0] = 1.
    for i in range(MC - 1):
        for _ in range(N - 1):
            qvoter_NNgroup(G, N, q, p, votes)
        timeline[i + 1] = qvoter_NNgroup(G, N, q, p, votes)

    return timeline


def qvoter_NNgroup(G, N, q, p, votes):  # zad1
    voter = rnd.randint(0, N - 1)

    if rnd.random() < p:  # independent
        if rnd.random() < 0.5:
            votes[voter] = - votes[voter]
    else:  # conformist
        subset = np.zeros(q, int)
        neighbours = G.neighbors(voter)
        for k in range(q):
            subset[k] = rnd.choice(neighbours)

        #####
        # sv_sum = 0
        # for node in subset:
        #     sv_sum += votes[node]
        #
        # if abs(sv_sum) == q:  # q-panel is unanimous
        #     votes[voter] = votes[subset[0]]
        first_vote = votes[subset[0]]
        i = 0
        while i < q and first_vote == votes[subset[i]]:
            i += 1
        if i == q:
            votes[voter] = first_vote

    return sum(votes) / N  # magnetization(N, votes)


def plot_results(x, y):
    plt.figure()
    plt.plot(x, y)
    plt.show()


def zad2(rep=1):
    N = 100
    MC = 1000
    Q = [3, 4]
    P = list(np.arange(0.0, 0.51, 0.02))
    for q in Q:
        m1, t1 = MCsimulations('complete', N, P, q, MC=MC, MC2=rep)
        m2, t2 = MCsimulations('bs', N, P, q, MC=MC, MC2=rep)
        m3, t3 = MCsimulations('ws1', N, P, q, MC=MC, MC2=rep)
        m4, t4 = MCsimulations('ws2', N, P, q, MC=MC, MC2=rep)


def zad3():
    zad2(100)


def zad4(N, MC, rep):
    Q = [3, 4]
    P_short = [0.0, 0.2, 0.5, 0.7]

    g = 'ws1'

    for q in Q:
        single_run_WS('t', g, N, P_short, q, MC, MC2=1, show=False)
        single_run_WS('t', g, N, P_short, q, MC, MC2=rep, show=False)
        # print(timer() - start)

        # plt.show()


def final_magnetization(N, MC, rep, q):
    # TODO ustawic wiecej p w P
    P_long = np.arange(0.1, 0.5, 0.02)

    # PART 1
    G = ['complete', 'ba', 'ws1', 'ws2']
    plt.figure()
    for g in G:
        m, t = MCsimulations(g, N, P_long, q, MC=MC, MC2=rep)
        plt.plot(P_long, m, '-')
    plt.legend(['complete graph', 'BA(100,4)', 'WS(100,4,0.01)', 'WS(100,4,0.2)'])
    plt.xlabel("independence probability")
    plt.ylabel("final magnetization")
    plt.title("Final magnetization for different topologies with q=" + str(q))
    plt.savefig("final_magn_" + str(q))
    plt.close()


def zad5a(N, MC, rep):
    print("part 1:")
    final_magnetization(N, MC, rep, 3)
    print(timer() - start)
    print("part 2:")
    final_magnetization(N, MC, rep, 4)
    print(timer() - start)


def zad5b(N, MC, rep):
    # TODO ustawic wiecej p w P
    P_long = np.arange(0.1, 0.5, 0.02)
    Q = [3, 4]

    plt.figure()
    leg = []
    for q in Q:
        m, t = MCsimulations('ws1', N, P_long, q, MC=MC, MC2=rep)
        plt.plot(P_long, m, '*-')
        leg.append('q = ' + str(q))

        print(timer() - start)

    plt.legend(leg)
    plt.xlabel("independence probability")
    plt.ylabel("final magnetization")
    plt.title("Final magnetization for WS(100,4,0.01)")

    plt.savefig("final_magn_ws001")
    plt.close()

    print(timer() - start)


def single_run_WS(plot_type, graph, N, P, q, MC, MC2=1, show=True):
    m, t = MCsimulations(graph, N, P, q, MC=MC, MC2=MC2)

    plt.figure()
    if plot_type == 'm':
        plt.plot(P, m, '*')
        plt.xlabel("independence probability")
    else:
        for res in t:
            plt.plot(res)
        plt.legend(P, loc='upper right')
        plt.xlabel("time")

    plt.ylabel("magnetization")
    if MC2 == 1:
        plt.title("magnetization of WS(100,4,0.01), single run, q=" + str(q))
    else:
        plt.title("magnetization of WS(100,4,0.01), average, q=" + str(q))

    if show is True:
        plt.show()
    else:
        plt.savefig("zad4_q=" + str(q) + "_rep=" + str(MC2) + ".png")


if __name__ == "__main__":
    # N = 100
    # MC = 1000
    # rep = 100

    N = 100
    MC = 10
    rep = 10

    start = timer()

    zad4(N, MC, rep)
    print(timer() - start)
    print("koniec zad 4")

    zad5a(N, MC, rep)
    print("koniec zad 5a")

    zad5b(N, MC, rep)
    print("koniec zad 5b")
