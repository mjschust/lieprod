package main

import (
	"log"
	"math/big"
	"net"

	"github.com/mjschust/cblocks/bundle"
	pb "github.com/mjschust/cblocks/cbservice"
	"github.com/mjschust/cblocks/lie"
	"golang.org/x/net/context"
	"google.golang.org/grpc"
	"google.golang.org/grpc/reflection"
)

const (
	port = ":50051"
)

type server struct{}

func (s *server) ComputeRank(ctx context.Context, cbr *pb.ConformalBlocksRequest) (*pb.IntReply, error) {
	cbb := constructCBBundle(cbr)
	rslt := cbb.Rank()

	return constructIntReply(rslt), nil
}

func (s *server) SymComputeRank(ctx context.Context, cbr *pb.SymConformalBlocksRequest) (*pb.IntReply, error) {
	cbb := processSymCBRequest(cbr)
	rslt := cbb.Rank()

	return constructIntReply(rslt), nil
}

func (s *server) SymComputeDivisor(ctx context.Context, cbr *pb.SymConformalBlocksRequest) (*pb.VectorReply, error) {
	cbb := processSymCBRequest(cbr)
	rslt := cbb.SymmetrizedDivisor()

	return constructVectorReply(rslt), nil
}

func constructCBBundle(cbr *pb.ConformalBlocksRequest) bundle.CBBundle {
	rank := int(cbr.Algebra.Rank)
	level := int(cbr.Level)
	wts := make([]lie.Weight, len(cbr.Weights))
	for i := range cbr.Weights {
		wt := make([]int, rank)
		for j := range wt {
			wt[j] = int(cbr.Weights[i].Coords[j])
		}
		wts[i] = wt
	}

	alg := lie.NewAlgebra(lie.NewTypeARootSystem(rank))
	return bundle.NewCBBundle(alg, wts, level)
}

func processSymCBRequest(cbr *pb.SymConformalBlocksRequest) bundle.SymCBBundle {
	rank := int(cbr.Algebra.Rank)
	level := int(cbr.Level)
	n := int(cbr.NumPoints)
	wt := make([]int, len(cbr.Weight.Coords))
	for i := range wt {
		wt[i] = int(cbr.Weight.Coords[i])
	}

	alg := lie.NewAlgebra(lie.NewTypeARootSystem(rank))
	return bundle.NewSymmetricCBBundle(alg, wt, level, n)
}

func constructVectorReply(rslt []*big.Rat) *pb.VectorReply {
	coords := make([]*pb.RatReply, len(rslt))
	for i := range rslt {
		coords[i] = constructRatReply(rslt[i])
	}

	return &pb.VectorReply{Coords: coords}
}

func constructRatReply(rslt *big.Rat) *pb.RatReply {
	numerator := constructIntReply(rslt.Num())
	denominator := constructIntReply(rslt.Denom())
	return &pb.RatReply{Numerator: numerator, Denominator: denominator}
}

func constructIntReply(rslt *big.Int) *pb.IntReply {
	if rslt.IsInt64() {
		return &pb.IntReply{Result: rslt.Int64()}
	}

	return &pb.IntReply{BigResult: rslt.Text(16)}
}

func main() {
	lis, err := net.Listen("tcp", port)
	if err != nil {
		log.Fatalf("failed to listen: %v", err)
	}
	s := grpc.NewServer()
	pb.RegisterCBlocksServer(s, &server{})
	// Register reflection service on gRPC server.
	reflection.Register(s)
	if err := s.Serve(lis); err != nil {
		log.Fatalf("failed to serve: %v", err)
	}
}
