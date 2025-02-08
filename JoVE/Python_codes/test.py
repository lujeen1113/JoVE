from dwave_qbsolv import QBSolv
def Solver(H,J):
    response = QBSolv().sample_ising(H,J)
    sample=list(response.samples())
    return sample


