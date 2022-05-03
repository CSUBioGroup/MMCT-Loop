from time import time
import hnswlib
import numpy as np
import leidenalg
import igraph as ig

class HNSW3:
    def __init__(self, mat, M, ef, k, le, resolution, minPts):
        self.M = int(M)
        self.ef = int(ef)
        self.k = int(k)
        self.le = int(le)
        self.resolution = float(resolution)
        self.minPts = int(minPts)
        self.fit(mat)

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

        g = ig.Graph(n=num_elements, edges=es, directed=True)
        print("igraph_time:", time() - t)

        # leiden
        t = time()
        if self.le == 0:
            part = leidenalg.ModularityVertexPartition(g)
        elif self.le == 1:
            part = leidenalg.CPMVertexPartition(g, resolution_parameter=self.resolution)
        elif self.ef == 2:
            part = leidenalg.RBConfigurationVertexPartition(g, resolution_parameter=self.resolution)
        optimiser = leidenalg.Optimiser()
        diff = optimiser.optimise_partition(part)
        print("leiden_time:", time() - t)

        # label
        t = time()
        label = np.zeros(num_elements, dtype=int)
        n_clusters = len(part)
        new_part = np.array(part, dtype=list)

        for i in range(n_clusters):
            # label[new_part[i]] = i
            if len(new_part[i]) >= self.minPts:
                label[list(new_part[i])] = i
            else:
                label[list(new_part[i])] = -1
        self.labels = dict(zip(data_labels, label))
        for i in range(num_elements):
            if self.labels[i] == -1:
                del self.labels[i]
        print("labels_time:", time() - t)
