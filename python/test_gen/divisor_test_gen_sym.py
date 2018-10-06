from __future__ import division
import conformal_blocks.cbbundle as cbd
import cProfile, time, random, subprocess, math
import cbclient.cbclient as cbc

def experiment():
    """
    Generates Macaulay 2 test cases, runs them using Swinarski's program, then outputs new unit tests if successful
    :return: Null
    """
    rank = 3
    level = 3
    num_points = 6

    client = cbc.CBClient()
    liealg = cbd.TypeALieAlgebra(rank)
    A_l = liealg.get_weights(level)
    m2file = open("TestRank.m2", "w")
    m2file.write("loadPackage(\"ConformalBlocks\");\n")
    m2file.write("sl_" + str(rank+1) + " = simpleLieAlgebra(\"A\", " + str(rank) + ");\n")
    test_cases = []
    for wt in A_l:
        cbb = cbd.SymmetricConformalBlocksBundle(client, liealg, wt, num_points, level)
        if cbb.get_rank() == 0:
            continue

        test_cases.append(wt)
        wt_str = "{"
        for i in range(num_points):
            if len(wt) == 1:
                wt_str += "{" + str(wt)[1] + "}, "
            else:
                wt_str += "{" + str(wt)[1:-1] + "}, "
        wt_str = wt_str[:-2] + "}"

        div = cbb.get_symmetrized_divisor()
        div_str = "{"
        for coord in div:
            div_str += str(coord*math.factorial(num_points)) + ", "
        div_str = div_str[:-2] + "}"

        m2file.write("V = conformalBlockVectorBundle(sl_" + str(rank+1) + ", " + str(level)  + ", " + wt_str + ", 0);\n")
        m2file.write("if " + div_str + " != coefficientList symmetrizedConformalBlockDivisor(V) then error(\"Bundle " + "(sl_" + str(rank+1) + ", " + str(level)  + ", " + wt_str + ") incorrect rank\");\n")

    m2file.write("print(\"OK\");\n")
    m2file.close()

    test_out = subprocess.check_output(["M2", "--script", "TestRank.m2"])
    if test_out == "OK\n":
        for wt in test_cases:
            cbb = cbd.SymmetricConformalBlocksBundle(client, liealg, wt, num_points, level)
            wt_str = ""
            if len(wt) == 1:
                wt_str += "{" + str(wt)[1] + "}, "
            else:
                wt_str += "{" + str(wt)[1:-1] + "}, "
            
            div = cbb.get_symmetrized_divisor()
            div_str = "[]*big.Rat{"
            for coord in div:
                div_str += "big.NewRat(" + str(coord.numerator) + ", " + str(coord.denominator) + ")" + ", "
            div_str = div_str[:-2] + "}"
            
            print("{lie.NewTypeARootSystem(" + str(rank) + "), lie.Weight" + wt_str + str(level) + ", " + str(num_points) + ", " + div_str + "},")
        print("OK")
    else:
        print(test_out)


if __name__ == '__main__':
    experiment()
