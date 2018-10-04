from __future__ import division
import conformal_blocks.cbbundle as cbd
import cbclient.cbclient as cbc
import random

def experiment():
    rank = 5
    level = 4
    num_points = 6
    tries = 10

    client = cbc.CBClient()
    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True, exact=False)
    A_l = liealg.get_weights(level)
    print("Weight", "Rank", "Divisor")
    for i in range(tries):
        weights = [random.choice(A_l) for i in range(num_points)]
        rk = client.request_rank(liealg.get_type(), liealg.rank, weights, level)
        print(weights, rk)

if __name__ == '__main__':
    experiment()