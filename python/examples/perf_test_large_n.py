from __future__ import division
import conformal_blocks.cbbundle as cbd
import cbclient.cbclient as cbc

def experiment():
    rank = 5
    level = 4
    num_points = 100

    client = cbc.CBClient()
    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True, exact=False)
    print("Weight", "Rank", "Divisor")
    for wt in liealg.get_weights(level):
        cbb = cbd.SymmetricConformalBlocksBundle(client, liealg, wt, num_points, level)
        if cbb.get_rank() == 0:
            #print(wt, 0)
            continue
        #divisor = cbb.get_symmetrized_divisor()
        print(wt, cbb.get_rank())

if __name__ == '__main__':
    experiment()