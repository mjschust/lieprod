from __future__ import print_function

import grpc
import cblocks_pb2
import cblocks_pb2_grpc
from conformal_blocks import lie, cbbundle

def run():
    with grpc.insecure_channel('localhost:50051') as channel:
        stub = cblocks_pb2_grpc.CBlocksStub(channel)
        rank = 6
        level = 4
        num_points = 10

        liealg = cbbundle.TypeALieAlgebra(rank, store_fusion=True, exact=False)
        print("Weight", "Rank")
        rkPromises = []
        for wt in liealg.get_weights(level):
            cbb = cbbundle.SymmetricConformalBlocksBundle(liealg, wt, num_points, level)
            rkPromises = rkPromises + [sendRankRequest(stub, cbb)]
            #if rk == 0:
            #    continue
            #print(wt, rk)

        for rkPromise in rkPromises:
            rk = rkPromise.result().result
            if rk == 0:
                continue
            print(rk)

def sendRankRequest(stub, cbb):
    r = cbb.liealg.rank
    l = cbb.level
    n = len(cbb.weights)
    wt = cblocks_pb2.Weight(coords=cbb.weights[0])
    alg = cblocks_pb2.LieAlgebra(type=0, rank=r)
    req = cblocks_pb2.SymConformalBlocksRequest(algebra=alg, weight=wt, num_points=n, level=l)
    promise = stub.ComputeRank.future(req)
    return promise


if __name__ == '__main__':
    run()