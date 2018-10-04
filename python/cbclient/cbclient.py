from __future__ import print_function

import grpc
import cblocks_pb2
import cblocks_pb2_grpc

class CBClient(object):
    def __init__(self, host='localhost', port='50051'):
        self.channel = grpc.insecure_channel(host + ":" + port)
        self.stub = cblocks_pb2_grpc.CBlocksStub(self.channel)

    def request_rank(self, type, rank, weights, level):
        wts = [cblocks_pb2.Weight(coords=weight) for weight in weights]
        alg = cblocks_pb2.LieAlgebra(type=type.value, rank=rank)
        req = cblocks_pb2.ConformalBlocksRequest(
            algebra=alg, 
            weights=wts,
            level=level)
        response = self.stub.ComputeRank(req)
        if not response.big_result:
            return response.result
        else:
            return long(response.big_result, 16)

    def request_sym_rank(self, type, rank, weight, num_points, level):
        wt = cblocks_pb2.Weight(coords=weight)
        alg = cblocks_pb2.LieAlgebra(type=type.value, rank=rank)
        req = cblocks_pb2.SymConformalBlocksRequest(
            algebra=alg, 
            weight=wt, 
            num_points=num_points, 
            level=level)
        response = self.stub.SymComputeRank(req)
        if not response.big_result:
            return response.result
        else:
            return long(response.big_result, 16)
        

    def __del__(self):
        self.channel.close()
        print("Channel closed.")
