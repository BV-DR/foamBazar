import PyFoam
from ideFoam.inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import DictProxy
from os.path import join
from copy import deepcopy
"""
  Convenience class to simply write "fvSheme"
"""


GAMG_prec_1 = DictProxy()
GAMG_prec_1["preconditioner"] = "GAMG"
GAMG_prec_1["tolerance"] = 1e-7
GAMG_prec_1["relTol"] = 0
GAMG_prec_1["smoother"] = "DICGaussSeidel"
GAMG_prec_1["nPreSweeps"] = 0
GAMG_prec_1["nPostSweeps"] = 2
GAMG_prec_1["nFinestSweeps"] = 2
GAMG_prec_1["cacheAgglomeration"] = "false"
GAMG_prec_1["nCellsInCoarsestLevel"] = 10
GAMG_prec_1["agglomerator"] = "faceAreaPair"
GAMG_prec_1["mergeLevels"] = 1


GAMG_prec_2 = DictProxy()
GAMG_prec_2["preconditioner"] = "GAMG"
GAMG_prec_2["tolerance"] = 1e-7
GAMG_prec_2["relTol"] = 0
GAMG_prec_2["nVcycles"] = 2
GAMG_prec_2["smoother"] = "DICGaussSeidel"
GAMG_prec_2["nPreSweeps"] = 2
GAMG_prec_2["nPostSweeps"] = 2
GAMG_prec_2["nFinestSweeps"] = 2
GAMG_prec_2["cacheAgglomeration"] = "true"
GAMG_prec_2["nCellsInCoarsestLevel"] = 10
GAMG_prec_2["agglomerator"] = "faceAreaPair"
GAMG_prec_2["mergeLevels"] = 1

GAMG_1 = DictProxy()
GAMG_1["solver"] = "GAMG"
GAMG_1["tolerance"] = 1e-8
GAMG_1["relTol"] = 0
GAMG_1["smoother"] = "DIC"
GAMG_1["nPreSweeps"] = 0
GAMG_1["nPostSweeps"] = 2
GAMG_1["nFinestSweeps"] = 2
GAMG_1["cacheAgglomeration"] = "true"
GAMG_1["nCellsInCoarsestLevel"] = 10
GAMG_1["agglomerator"] = "faceAreaPair"
GAMG_1["mergeLevels"] = 1
GAMG_1["minIter"] = 2


GAMG_2 = DictProxy()
GAMG_2["solver"] = "GAMG"
GAMG_2["tolerance"] = 1e-8
GAMG_2["relTol"] = 0
GAMG_2["smoother"] = "GaussSeidel"
GAMG_2["cacheAgglomeration"] = "true"
GAMG_2["nCellsInCoarsestLevel"] = 10
GAMG_2["agglomerator"] = "faceAreaPair"
GAMG_2["mergeLevels"] = 1


SMOOTHSOLVER_1 = DictProxy()
SMOOTHSOLVER_1["solver"] = "smoothSolver"
SMOOTHSOLVER_1["smoother"] = "GaussSeidel"
SMOOTHSOLVER_1["tolerance"] = 1e-8
SMOOTHSOLVER_1["relTol"] = 0
SMOOTHSOLVER_1["nSweeps"] = 2
SMOOTHSOLVER_1["minIter"] = 2

PCG_1 = DictProxy()
PCG_1["solver"] = "PCG"
PCG_1["preconditioner"] = GAMG_prec_2
PCG_1["tolerance"] = 1e-8
PCG_1["relTol"] = 0
PCG_1["maxIter"] = 1000
PCG_1["minIter"] = 2


linearSolverDict = {
                     "PCG_1" : PCG_1,
                     "GAMG_1" : GAMG_1,
                     "SMOOTHSOLVER_1" : SMOOTHSOLVER_1,
                   }


class FvSolution(ReadWriteFile) :
    """
        FvSchemes dictionnary
    """

    @classmethod
    def Build(cls , case, fsiTol = 1e-8, useEulerCells=False, nOuterCorrectors=5, nInnerCorrectors = 4, application = "foamStar",
                    pressureSolver = "PCG_1",
                    velocitySolver = "SMOOTHSOLVER_1",
                    airDamping = 2.5
                    ) :

        res = cls( name = join(case, getFilePath("fvSolution") ), read = False )

        solvers = DictProxy()

        #Alpha
        alpha = DictProxy()
        alpha["nAlphaCorr"] = 3
        alpha["nAlphaSubCycles"] = 1
        alpha["cAlpha"] = 0.3
        alpha["MULESCorr"] = "yes"
        alpha["nLimiterIter"] = 5
        alpha["alphaApplyPrevCorr"] = "no"
        alpha["solver"] = "smoothSolver"
        alpha["smoother"] = "symGaussSeidel"
        alpha["tolerance"] = 1e-8
        alpha["relTol"] = 0
        alpha["minIter"] = 2
        solvers['"alpha.water.*"'] = alpha

        #Pressure
        solvers['p_rgh'] = deepcopy(linearSolverDict[pressureSolver])
        solvers['"(p_rghFinal|pcorr|pcorrFinal)"'] = deepcopy(linearSolverDict[pressureSolver])

        #Velocity
        solvers['"(U|k|epsilon)"'] = deepcopy(linearSolverDict[velocitySolver])
        solvers['"(U|k|epsilon)Final"'] = deepcopy(linearSolverDict[velocitySolver])

        #CellDisplacement
        solvers['"cellDisplacement.*"'] = deepcopy(GAMG_2)

        #Iteration
        pimp = DictProxy()
        if useEulerCells:
           pimp["EulerCells"] = "EulerCells"
        pimp["momentumPredictor"] = "yes"
        pimp["nOuterCorrectors"] = nOuterCorrectors
        pimp["fsiTol"] = fsiTol
        pimp["fsiMaxIter"] = 23
        pimp["nCorrectors"] = nInnerCorrectors
        pimp["nNonOrthogonalCorrectors"] = 0
        pimp["correctPhi"] = "no"
        pimp["moveMeshOuterCorrectors"] = "yes"
        pimp["airDamping"] = airDamping

        relax = { "equations" : { "U": 1, "UFinal" : 1, "p_rgh" : 1, "p_rghFinal" : 1 } }
        res["solvers"] = solvers
        res["PIMPLE"] = pimp
        res["relaxationFactors"] = relax
        return res


if __name__ == "__main__" :
   print(FvSolution.Build("test" , application = "foamStar"))
