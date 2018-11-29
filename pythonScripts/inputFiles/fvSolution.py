import PyFoam
from inputFiles import ReadWriteFile
from PyFoam.Basics.DataStructures import DictProxy
from os.path import join

"""
  Convenience class to simply write "fvSheme"
"""

class FvSolution(ReadWriteFile) :
    """
        FvSchemes dictionnary
    """
    
    @classmethod
    def Build(cls , case, fsiTol = 1e-8, useEuler=False, nOuterCorrectors=5, version = "foamStar") :

        res = cls ( name = join(case, "system" , "fvSolution" ), read = False )

        #-------- ddtSchemes
        solvers = DictProxy()
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
        
        pcorr = DictProxy()
        pcorr["solver"] = "PCG"
        prec = DictProxy()
        prec["preconditioner"] = "GAMG"
        prec["tolerance"] = 1e-7
        prec["relTol"] = 0
        prec["smoother"] = "DICGaussSeidel"
        prec["nPreSweeps"] = 0
        prec["nPostSweeps"] = 2
        prec["nFinestSweeps"] = 2
        prec["cacheAgglomeration"] = "false"
        prec["nCellsInCoarsestLevel"] = 10
        prec["agglomerator"] = "faceAreaPair"
        prec["mergeLevels"] = 1
        pcorr["preconditioner"] = prec
        pcorr["tolerance"] = 1e-8
        pcorr["relTol"] = 0
        pcorr["maxIter"] = 1000
        pcorr["minIter"] = 2
        solvers['"pcorr.*"'] = pcorr
        
        prgh = DictProxy()
        prgh["solver"] = "GAMG"
        prgh["tolerance"] = 1e-8
        prgh["relTol"] = 0
        prgh["smoother"] = "DIC"
        prgh["nPreSweeps"] = 0
        prgh["nPostSweeps"] = 2
        prgh["nFinestSweeps"] = 2
        prgh["cacheAgglomeration"] = "true"
        prgh["nCellsInCoarsestLevel"] = 10
        prgh["agglomerator"] = "faceAreaPair"
        prgh["mergeLevels"] = 1
        prgh["minIter"] = 2
        solvers['p_rgh'] = prgh

        prghf = DictProxy()
        prghf["solver"] = "PCG"
        prec = DictProxy()
        prec["preconditioner"] = "GAMG"
        prec["tolerance"] = 1e-7
        prec["relTol"] = 0
        prec["nVcycles"] = 2
        prec["smoother"] = "DICGaussSeidel"
        prec["nPreSweeps"] = 2
        prec["nPostSweeps"] = 2
        prec["nFinestSweeps"] = 2
        prec["cacheAgglomeration"] = "true"
        prec["nCellsInCoarsestLevel"] = 10
        prec["agglomerator"] = "faceAreaPair"
        prec["mergeLevels"] = 1
        prghf["preconditioner"] = prec
        prghf["tolerance"] = 1e-8
        prghf["relTol"] = 0
        prghf["maxIter"] = 1000
        prghf["minIter"] = 2
        solvers['p_rghFinal'] = prghf
        
        uke = DictProxy()
        uke["solver"] = "smoothSolver"
        uke["smoother"] = "GaussSeidel"
        uke["tolerance"] = 1e-8
        uke["relTol"] = 0
        uke["nSweeps"] = 2
        uke["minIter"] = 2
        solvers['"(U|k|epsilon)"'] = uke
        solvers['"(U|k|epsilon)Final"'] = uke
        
        cdisp = DictProxy()
        cdisp["solver"] = "GAMG"
        cdisp["tolerance"] = 1e-8
        cdisp["relTol"] = 0
        cdisp["smoother"] = "GaussSeidel"
        cdisp["cacheAgglomeration"] = "true"
        cdisp["nCellsInCoarsestLevel"] = 10
        cdisp["agglomerator"] = "faceAreaPair"
        cdisp["mergeLevels"] = 1
        solvers['"cellDisplacement.*"'] = cdisp

        pimp = DictProxy()
        if not useEuler: pimp["EulerCells"] = "EulerCells"
        pimp["momentumPredictor"] = "yes"
        pimp["nOuterCorrectors"] = nOuterCorrectors
        pimp["fsiTol"] = fsiTol
        pimp["fsiMaxIter"] = 23
        pimp["nCorrectors"] = 4
        pimp["nNonOrthogonalCorrectors"] = 0
        pimp["correctPhi"] = "no"
        pimp["moveMeshOuterCorrectors"] = "yes"
        
        relax = { "equations" : { "U": 1, "UFinal" : 1, "p_rgh" : 1, "p_rghFinal" : 1 } }
                
        res["solvers"] = solvers
        res["PIMPLE"] = pimp
        res["relaxationFactors"] = relax
        return res
        

if __name__ == "__main__" :
   print(FvSolution.Build("test" , version = "foamStar"))
