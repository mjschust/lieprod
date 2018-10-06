from __future__ import print_function, division

import grpc, fractions
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
        reply = self.stub.ComputeRank(req)
        return self._process_int_reply(reply)

    def request_sym_rank(self, type, rank, weight, num_points, level):
        wt = cblocks_pb2.Weight(coords=weight)
        alg = cblocks_pb2.LieAlgebra(type=type.value, rank=rank)
        req = cblocks_pb2.SymConformalBlocksRequest(
            algebra=alg, 
            weight=wt, 
            num_points=num_points, 
            level=level)
        reply = self.stub.SymComputeRank(req)
        return self._process_int_reply(reply)

    def request_sym_divisor(self, type, rank, weight, num_points, level):
        wt = cblocks_pb2.Weight(coords=weight)
        alg = cblocks_pb2.LieAlgebra(type=type.value, rank=rank)
        req = cblocks_pb2.SymConformalBlocksRequest(
            algebra=alg, 
            weight=wt, 
            num_points=num_points, 
            level=level)
        reply = self.stub.SymComputeDivisor(req)

        return self._process_vector_reply(reply)

    def _process_vector_reply(self, reply):
        return [self._process_rat_reply(x) for x in reply.coords]

    def _process_rat_reply(self, reply):
        numerator = self._process_int_reply(reply.numerator)
        denominator = self._process_int_reply(reply.denominator)
        return fractions.Fraction(numerator, denominator)

    def _process_int_reply(self, reply):
        if not reply.big_result:
            return reply.result
        else:
            return long(reply.big_result, 16)
        

    def __del__(self):
        self.channel.close()
        print("Channel closed.")
