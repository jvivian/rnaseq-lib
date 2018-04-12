#!/usr/bin/python -tt

# -*- coding: utf-8 -*-

"""
An implementation of the TriMap method for dimensionality reduction.

Usage:
    Y = trimap(X, num_dims, kin, kout, krand, eta, Yinit)

Parameters:
    num_dim   ------   output dimension
    kin       ------   number of knn points
    kout      ------   number of outliers (num triplets per point = kin * kout)
    krand     ------   number of random triplets per point
    eta       ------   initial learning rate
    Yinit     ------   initial solution

requirements:
    numpy
    sklearn
    annoy
    numba

Reference:
    Ehsan Amid and Manfred K. Warmuth, A more globally accurate dimensionality
    reduction method using triplets, arXiv preprint arXiv:1803.00854, 2018.
"""

import numba
from annoy import AnnoyIndex
import sys
from sklearn.neighbors import NearestNeighbors as knn
from sklearn.decomposition import TruncatedSVD
import numpy as np
import time


@numba.njit()
def euclid_dist(x1, x2):
    result = 0.0
    for i in range(x1.shape[0]):
        result += (x1[i] - x2[i]) ** 2
    return np.sqrt(result)


@numba.njit()
def rejection_sample(n_samples, max_int, rejects):
    result = np.empty(n_samples, dtype=np.int64)
    for i in range(n_samples):
        reject_sample = True
        while reject_sample:
            j = np.random.randint(max_int)
            for k in range(i):
                if j == result[k]:
                    break
            for k in range(rejects.shape[0]):
                if j == rejects[k]:
                    break
            else:
                reject_sample = False
        result[i] = j
    return result


def knn_annoy(X, num_neighbors=50):
    n = X.shape[0]
    tree = AnnoyIndex(X.shape[1])
    for i in xrange(n):
        tree.add_item(i, X[i, :])
    tree.build(10)
    nbrs = np.empty((n, num_neighbors), dtype=np.int64)
    distances = np.empty((n, num_neighbors), dtype=np.float64)
    for i in xrange(n):
        nbrs[i, :] = tree.get_nns_by_item(i, num_neighbors)
        for j in xrange(num_neighbors):
            distances[i, j] = tree.get_distance(i, nbrs[i, j])
    return (nbrs, distances, tree)


@numba.njit('i8[:,:](f8[:,:],i8[:,:],f8[:,:], i8,i8)', parallel=True, nogil=True)
def sample_knn_triplets(P, nbrs, distances, kin, kout):
    n, num_neighbors = nbrs.shape
    triplets = np.empty((n * kin * kout, 3), dtype=np.int64)
    for i in xrange(n):
        sort_indices = np.argsort(-P[i, :])
        for j in range(kin):
            sim = nbrs[i, sort_indices[j + 1]]
            samples = rejection_sample(kout, n, sort_indices[j + 1:])
            for k in range(kout):
                index = i * kin * kout + j * kout + k
                out = samples[k]
                triplets[index, 0] = i
                triplets[index, 1] = sim
                triplets[index, 2] = out
    return triplets


@numba.njit('f8[:,:](f8[:,:],i8,f8[:])', parallel=True, nogil=True)
def sample_random_triplets(X, krand, sig):
    n = X.shape[0]
    rand_triplets = np.empty((n * krand, 4), dtype=np.float64)
    for i in xrange(n):
        for j in range(krand):
            sim = np.random.choice(n)
            while sim == i:
                sim = np.random.choice(n)
            out = np.random.choice(n)
            while out == i or out == sim:
                out = np.random.choice(n)
            p_sim = np.exp(-euclid_dist(X[i, :], X[sim, :]) ** 2 / (sig[i] * sig[sim]))
            if p_sim < 1e-20:
                p_sim = 1e-20
            p_out = np.exp(-euclid_dist(X[i, :], X[out, :]) ** 2 / (sig[i] * sig[out]))
            if p_out < 1e-20:
                p_out = 1e-20
            if p_sim < p_out:
                sim, out = out, sim
                p_sim, p_out = p_out, p_sim
            rand_triplets[i * krand + j, 0] = i
            rand_triplets[i * krand + j, 1] = sim
            rand_triplets[i * krand + j, 2] = out
            rand_triplets[i * krand + j, 3] = p_sim / p_out
    return rand_triplets


@numba.njit('f8[:,:](f8[:,:],f8[:],i8[:,:])', parallel=True, nogil=True)
def find_p(distances, sig, nbrs):
    n, num_neighbors = distances.shape
    P = np.zeros((n, num_neighbors), dtype=np.float64)
    for i in range(n):
        for j in range(num_neighbors):
            P[i, j] = np.exp(-distances[i, j] ** 2 / sig[i] / sig[nbrs[i, j]])
    return P


@numba.njit('f8[:](i8[:,:],f8[:,:],i8[:,:],f8[:],f8[:])', parallel=True, nogil=True)
def find_weights(triplets, P, nbrs, distances, sig):
    num_triplets = triplets.shape[0]
    weights = np.empty(num_triplets, dtype=np.float64)
    for t in xrange(num_triplets):
        i = triplets[t, 0]
        sim = 0
        while (nbrs[i, sim] != triplets[t, 1]):
            sim += 1
        p_sim = P[i, sim]
        p_out = np.exp(-distances[t] ** 2 / (sig[i] * sig[triplets[t, 2]]))
        if p_out < 1e-20:
            p_out = 1e-20
        weights[t] = p_sim / p_out
    return weights


def generate_triplets(X, kin, kout, krand):
    n, dim = X.shape
    num_neighbors = max(kin, 150)
    exact = n <= 1e4 or dim <= 50
    if exact:  # do exact knn search
        if dim > 50:
            X = TruncatedSVD(n_components=50, random_state=0).fit_transform(X)
        knn_tree = knn(n_neighbors=num_neighbors, algorithm='auto').fit(X)
        distances, nbrs = knn_tree.kneighbors(X)
        distances = np.empty((n, num_neighbors), dtype=np.float64)
        for i in xrange(n):
            for j in xrange(num_neighbors):
                distances[i, j] = euclid_dist(X[i, :], X[nbrs[i, j], :])
    else:  # use annoy
        tree = AnnoyIndex(dim)
        for i in xrange(n):
            tree.add_item(i, X[i, :])
        tree.build(10)
        nbrs = np.empty((n, num_neighbors), dtype=np.int64)
        distances = np.empty((n, num_neighbors), dtype=np.float64)
        for i in xrange(n):
            nbrs[i, :] = tree.get_nns_by_item(i, num_neighbors)
            for j in xrange(num_neighbors):
                distances[i, j] = tree.get_distance(i, nbrs[i, j])
    print("found nearest neighbors")
    sig = np.maximum(np.mean(distances[:, 10:20], axis=1), 1e-20)  # scale parameter
    P = find_p(distances, sig, nbrs)
    triplets = sample_knn_triplets(P, nbrs, distances, kin, kout)
    num_triplets = triplets.shape[0]
    outlier_dist = np.empty(num_triplets, dtype=np.float64)
    if exact:
        for t in xrange(num_triplets):
            outlier_dist[t] = np.sqrt(np.sum((X[triplets[t, 0], :] - X[triplets[t, 2], :]) ** 2))
    else:
        for t in xrange(num_triplets):
            outlier_dist[t] = tree.get_distance(triplets[t, 0], triplets[t, 2])
    weights = find_weights(triplets, P, nbrs, outlier_dist, sig)
    if krand > 0:
        rand_triplets = sample_random_triplets(X, krand, sig)
        rand_weights = rand_triplets[:, -1]
        rand_triplets = rand_triplets[:, :-1].astype(np.int64)
        triplets = np.vstack((triplets, rand_triplets))
        weights = np.hstack((weights, rand_weights))
    weights /= np.max(weights)
    weights += 0.0001
    return (triplets, weights)


@numba.njit('void(f8[:,:],f8[:,:],f8)', parallel=True, nogil=True)
def update_embedding(Y, grad, lr):
    n, dim = Y.shape
    for i in xrange(n):
        for d in xrange(dim):
            Y[i, d] -= lr * grad[i, d]


@numba.njit('f8[:,:](f8[:,:],i8,i8,i8[:,:],f8[:])', parallel=True, nogil=True)
def trimap_grad(Y, kin, kout, triplets, weights):
    n, dim = Y.shape
    num_triplets = triplets.shape[0]
    grad = np.zeros((n, dim), dtype=np.float64)
    y_ij = np.empty(dim, dtype=np.float64)
    y_ik = np.empty(dim, dtype=np.float64)
    num_viol = 0.0
    loss = 0.0
    num_knn_triplets = n * kin * kout
    for t in xrange(num_triplets):
        i = triplets[t, 0]
        j = triplets[t, 1]
        k = triplets[t, 2]
        if (t % kout) == 0 or (t >= num_knn_triplets):  # update y_ij, d_ij
            d_ij = 1.0
            d_ik = 1.0
            for d in xrange(dim):
                y_ij[d] = Y[i, d] - Y[j, d]
                y_ik[d] = Y[i, d] - Y[k, d]
                d_ij += y_ij[d] ** 2
                d_ik += y_ik[d] ** 2
        else:
            d_ik = 1.0
            for d in xrange(dim):
                y_ik[d] = Y[i, d] - Y[k, d]
                d_ik += y_ik[d] ** 2
        if (d_ij > d_ik):
            num_viol += 1.0
        loss += weights[t] * 1.0 / (1.0 + d_ik / d_ij)
        w = weights[t] / (d_ij + d_ik) ** 2
        for d in xrange(dim):
            gs = y_ij[d] * d_ik * w
            go = y_ik[d] * d_ij * w
            grad[i, d] += gs - go
            grad[j, d] -= gs
            grad[k, d] += go
    last = np.zeros((1, dim), dtype=np.float64)
    last[0] = loss
    last[1] = num_viol
    return np.vstack((grad, last))


def trimap(X, num_dims=2, kin=50, kout=5, krand=5, eta=10000.0, Yinit=[]):
    t = time.time()
    n, dim = X.shape
    print "running TriMap on %d points with dimension %d" % (n, dim)
    print("pre-processing")
    X -= np.min(X)
    X /= np.max(X)
    X -= np.mean(X, axis=0)
    triplets, weights = generate_triplets(X, kin, kout, krand)
    print("sampled triplets")

    if np.size(Yinit) > 0:
        Y = Yinit
    else:
        Y = np.random.normal(size=[n, num_dims]) * 0.0001
    C = np.inf
    tol = 1e-7
    num_iter = 1500
    num_triplets = float(triplets.shape[0])

    print("running TriMap")
    for itr in range(num_iter):
        old_C = C
        grad = trimap_grad(Y, kin, kout, triplets, weights)
        C = grad[-1, 0]
        num_viol = grad[-1, 1]

        # update Y
        update_embedding(Y, grad, eta * n / num_triplets)
        # Y -= (eta * n/num_triplets) * grad[:-1,:]

        # update the learning rate
        if old_C > C + tol:
            eta = eta * 1.05
        else:
            eta = eta * 0.5

        if (itr + 1) % 100 == 0:
            print 'Iteration: %4d, Loss: %3.3f, Violated triplets: %0.4f' % (
            itr + 1, C, num_viol / num_triplets * 100.0)
    print "Elapsed time %s" % (time.time() - t)
    return Y


def load_known_size(filename):
    first = True
    with open(filename) as f:
        for irow, line in enumerate(f):
            if first:
                nrow, ncol = [int(a) for a in line.split(',')]
                x = np.empty((nrow, ncol), dtype=np.double)
                first = False
            else:
                x[irow - 1, :] = [float(a) for a in line.split(',')]
    return x


def main():
    filename = sys.argv[1]
    X = load_known_size(filename)
    X = X.astype(np.float64)
    Y = trimap(X, 2, 50, 5)
    np.savetxt('result.txt', Y)


if __name__ == '__main__':
    main()
