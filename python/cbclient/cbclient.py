from __future__ import print_function

import grpc

import cblocks_pb2
import cblocks_pb2_grpc


def run():
    with grpc.insecure_channel('localhost:50051') as channel:
        stub = cblocks_pb2_grpc.CBlocksStub(channel)
        r = 2
        l = 2
        n = 5
        wt = cblocks_pb2.Weight(coords=[1,1])
        alg = cblocks_pb2.LieAlgebra(type=0, rank=r)
        req = cblocks_pb2.SymConformalBlocksRequest(algebra=alg, weight=wt, num_points=n, level=l)
        response = stub.ComputeRank(req)
    print("Greeter client received: " + str(response.result))


if __name__ == '__main__':
    run()