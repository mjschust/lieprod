package main

import (
	"log"
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

func (s *server) Sum(ctx context.Context, in *pb.Weight) (*pb.IntReply, error) {
	coords := in.GetCoords()
	var sum int64
	for i := range coords {
		sum += coords[i]
	}
	return &pb.IntReply{Result: sum}, nil
}

func (s *server) ComputeRank(ctx context.Context, cbr *pb.SymConformalBlocksRequest) (*pb.IntReply, error) {
	rank := int(cbr.Algebra.Rank)
	level := int(cbr.Level)
	n := int(cbr.NumPoints)
	wt := make([]int, len(cbr.Weight.Coords))
	for i := range wt {
		wt[i] = int(cbr.Weight.Coords[i])
	}

	alg := lie.NewAlgebra(lie.NewTypeARootSystem(rank))
	cbb := bundle.NewSymmetricCBBundle(alg, wt, level, n)
	rslt := cbb.Rank()

	if rslt.IsInt64() {
		return &pb.IntReply{Result: rslt.Int64()}, nil
	}

	return &pb.IntReply{BigResult: rslt.Text(16)}, nil
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
