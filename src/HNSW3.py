from time import time
import hnswlib
import numpy as np
import queue

class HNSW3:
    def __init__(self, mat, M, ef, k, base_step):
        self.M = int(M)
        self.ef = int(ef)
        self.k = int(k)
        self.base_step = int(base_step)
        self.fit(mat)
        self.g = []

    def get_base(self, node, mx_lev):
        set_vis = set()
        base = queue.Queue()
        ret = [node]
        base.put([node, 0])
        set_vis.add(node)
        while base.empty() is False:
            cur, step = base.get()
            for nxt in self.g[cur]:
                if nxt in set_vis:
                    continue
                if step < mx_lev:
                    set_vis.add(nxt)
                    base.put([nxt, step + 1])
                    ret.append(nxt)
        return ret

    def shrink_point(self, es, mx_nodes, mx_lev=3, rate=1.5):
        self.g = [[] for i in range(mx_nodes)]
        t_label = [-1 for i in range(mx_nodes)]
        for edge in es:
            s, t = edge
            self.g[s].append(t)
            self.g[t].append(s)
        vis = [False for i in range(mx_nodes)]
        for node in range(mx_nodes):
            if vis[node]:
                continue
            base = self.get_base(node, mx_lev)
            sum_edges = 0
            central_nodes = set()
            for b in base:
                c_nodes = [x for x in self.g[b] if x in base]
                sum_edges += len(c_nodes)
                for x in c_nodes:
                    if vis[x] is False and len(c_nodes) > 1:
                        central_nodes.add(x)
            if sum_edges / len(base) >= rate:
                for x in central_nodes:
                    vis[x] = True
                    t_label[x] = node
        return t_label

    def fit(self, mat):
        self.labels = {}
        if len(mat) <= self.k:
            self.labels = dict(zip(mat[:, 0:1].flatten(), np.zeros(len(mat), dtype=int)))
            print("!!!!!!!len<=k")
            return

        # hnsw
        t = time()
        num_elements = len(mat)
        dim = mat.shape[1] - 1
        data = mat[:, 1:]
        data_labels = mat[:, 0:1].flatten()
        p = hnswlib.Index(space='l2', dim=dim)
        p.init_index(max_elements=num_elements, ef_construction=self.ef, M=self.M)
        p.add_items(data, data_labels)
        p.set_ef(self.ef)
        labels, distances = p.knn_query(data, k=self.k)
        print("hnsw_time:", time() - t)

        # convert to igraph directly
        t = time()
        es = []
        for i in range(num_elements):
            s = labels[i][0]
            for j in range(1, self.k):
                e = labels[i][j]
                es.append([s, e])

        # g = ig.Graph(n=num_elements, edges=es, directed=True)
        # print("igraph_time:", time() - t)

        # leiden
        # t = time()
        # if self.le == 0:
        #     part = leidenalg.ModularityVertexPartition(g)
        # elif self.le == 1:
        #     part = leidenalg.CPMVertexPartition(g, resolution_parameter=self.resolution)
        # elif self.ef == 2:
        #     part = leidenalg.RBConfigurationVertexPartition(g, resolution_parameter=self.resolution)
        # optimiser = leidenalg.Optimiser()
        # diff = optimiser.optimise_partition(part)
        # print("leiden_time:", time() - t)

        # label
        # t = time()
        # label = np.zeros(num_elements, dtype=int)
        # n_clusters = len(part)
        # new_part = np.array(part, dtype=list)
        #
        # for i in range(n_clusters):
        #     # label[new_part[i]] = i
        #     if len(new_part[i]) >= self.minPts:
        #         label[list(new_part[i])] = i
        #     else:
        #         label[list(new_part[i])] = -1
        label = self.shrink_point(es, num_elements, mx_lev=self.base_step)
        self.labels = dict(zip(data_labels, label))

        for i in range(num_elements):
            if self.labels[i] == -1:
                del self.labels[i]
        print("subnetting_cluster_time:", time() - t)
