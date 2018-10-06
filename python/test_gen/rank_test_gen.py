from __future__ import division
import conformal_blocks.cbbundle as cbd
import cProfile, time, random, subprocess
import cbclient.cbclient as cbc

def experiment():
    """
    Generates Macaulay 2 test cases, runs them using Swinarski's program, then outputs new unit tests if successful
    :return: Null
    """
    rank = 3
    level = 3
    num_points = 3
    tries = 10

    client = cbc.CBClient()
    liealg = cbd.TypeALieAlgebra(rank)
    A_l = liealg.get_weights(level)
    m2file = open("TestRank.m2", "w")
    m2file.write("loadPackage(\"ConformalBlocks\");\n")
    m2file.write("sl_" + str(rank+1) + " = simpleLieAlgebra(\"A\", " + str(rank) + ");\n")
    test_cases = []
    for i in range(tries):
        weights = [random.choice(A_l) for i in range(num_points)]
        test_cases.append(weights)
        cbb = cbd.ConformalBlocksBundle(client, liealg, weights, level)
        wt_str = "{"
        for wt in weights:
            if len(wt) == 1:
                wt_str += "{" + str(wt)[1] + "}, "
            else:
                wt_str += "{" + str(wt)[1:-1] + "}, "
        wt_str = wt_str[:-2] + "}"

        m2file.write("V = conformalBlockVectorBundle(sl_" + str(rank+1) + ", " + str(level)  + ", " + wt_str + ", 0);\n")
        m2file.write("if " + str(cbb.get_rank()) + " != conformalBlockRank(V) then error(\"Bundle " + "(sl_" + str(rank+1) + ", " + str(level)  + ", " + wt_str + ") incorrect rank\");\n")

    m2file.write("print(\"OK\");\n")
    m2file.close()

    test_out = subprocess.check_output(["M2", "--script", "TestRank.m2"])
    if test_out == "OK\n":
        for case in test_cases:
            cbb = cbd.ConformalBlocksBundle(client, liealg, case, level)
            wts_str = "[]lie.Weight{"
            for wt in case:
                if len(wt) == 1:
                    wts_str += "lie.Weight{" + str(wt)[1] + "}, "
                else:
                    wts_str += "lie.Weight{" + str(wt)[1:-1] + "}, "
            wts_str = wts_str[:-2] + "}, "
            
            print("{lie.NewTypeARootSystem(" + str(rank) + "), " + wts_str + str(level) + ", big.NewInt(" + str(cbb.get_rank()) + ")},")
        print("OK")
    else:
        print(test_out)


if __name__ == '__main__':
    experiment()
