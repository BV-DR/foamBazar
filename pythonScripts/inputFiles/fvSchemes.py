from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import DictProxy
from os.path import join
from inputFiles.compatOF import alpha, p_rgh

"""
  Convenience class to simply write "fvSheme"
"""

class FvSchemes(WriteParameterFile) :
    """
        FvSchemes dictionnary
    """
    def __init__(self , case, version = "foamStar", prsJump = False,  orthogonalCorrection = False, blendCN = 0.9, simType = "steady", limitedGrad=False):

        WriteParameterFile.__init__(self,  name = join(case, "system" , "fvSchemes" )  )

        #-------- ddtSchemes
        ddt = DictProxy()
        if simType.lower()=="steady":
            ddt["default"]          = "steadyState"
            ddt["ddt(U)"]           = "Euler"
            ddt["ddt(alpha.water)"] = "Euler"
        elif simType.lower()=="euler":
            ddt["default"] = "Euler"
            ddt["ddt(U)"]  = "Euler"
        else:
            ddt["default"] = "CrankNicolson {}".format(blendCN)
            ddt["ddt(U)"]  = "Euler"
        self["ddtSchemes"] = ddt

        #-------- gradSchemes
        grad = DictProxy()
        if limitedGrad:
            grad["default"] = "cellLimited leastSquares 1"
            grad["limitedGrad"] = "cellLimited Gauss linear 1"
        else:
            grad["default"] = "Gauss linear"
            if prsJump : grad["snGradCorr(pd)"] = "interfaceGauss linear"
        self["gradSchemes"] = grad

        #-------- divSchemes
        if simType=="steady": bounded = "bounded"
        else: bounded = ""
        div = DictProxy()
        div["div(rhoPhi,U)"] = "{} Gauss linearUpwind GradU".format(bounded)

        if version == "foamStar" :
            div["div(phi,alpha)"]   = "Gauss vanLeer"            #vanLeer01DC not in openCFD version
            div["div(phirb,alpha)"] = "Gauss interfaceCompression" #From Vuko's opinion "Gauss linear" should not be used
            if limitedGrad:
                div["div(phi,k)"] = "Gauss linearUpwind limitedGrad"
                div["div(phi,omega)"] = "Gauss linearUpwind limitedGrad"
            div["div((muEff*dev(T(grad(U)))))"] = "Gauss linear"
            div["div(((rho*nuEff)*dev2(T(grad(U)))))"] = "Gauss linear"
        else :
            div["div(phi,alpha)"]   = "Gauss vanLeer01DC"
            div["div(phirb,alpha)"] = "Gauss vofCompression" #From Vuko's opinion "Gauss linear" should not be used
            if version == "swenseFoam" :
                div["div(phi,UDiff)"]        = "Gauss linearUpwind Gauss linear"
                div["div(phi,levelSetDiff)"] = "Gauss vanLeerDC"
                div["div(phi,UInc)"]         = "Gauss linear"
                div["div(phi,levelSetInc)"]  = "Gauss linear"
        self["divSchemes"] = div

        #-------- laplacianSchemes
        lapl = DictProxy()
        if orthogonalCorrection:
            if type(orthogonalCorrection) == "float" :
                lapl["default"] = "Gauss linear limited {}".format(orthogonalCorrection)
                if prsJump : 
                    lapl["laplacian(rAU,pd)"] = "interfaceGauss linear {}".format(orthogonalCorrection)
                    lapl["laplacian(rAU,p)"] = "interfaceGauss linear {}".format(orthogonalCorrection)
            elif orthogonalCorrection.lower() == "implicit" : 
                lapl["default"] = "Gauss linear corrected"
            else : 
                print ("orthogonalCorrection" , orthogonalCorrection , "not recognized")
        else :
            lapl["default"] = "Gauss linear uncorrected"
            if prsJump : 
                lapl["laplacian(rAU,pd)"] = "interfaceGauss linear interfaceUncorrected"
                lapl["laplacian(rAU,p)"] = "interfaceGauss linear interfaceUncorrected"
        self["laplacianSchemes"] = lapl

        #-------- interpolationSchemes
        interp = DictProxy()
        interp["default"] = "linear"
        # interp["interpolate(y)"] = "linear"
        self["interpolationSchemes"] = interp

        #-------- snGradSchemes
        sng = DictProxy()
        if orthogonalCorrection :
            if type(orthogonalCorrection) == "float" :
                sng["default"] = "limited {}".format(orthogonalCorrection)
                if prsJump : sng["snGrad(pd)"] = "interfaceLimited 0.5"
            elif orthogonalCorrection.lower() == "implicit" : 
                sng["default"] = "corrected"
            else :
                sng["default"] = "uncorrected"
                if prsJump : sng["snGrad(pd)"] = "interfaceUncorrected"
        self["snGradSchemes"] = sng

        #-------- fluxRequired
        flux = DictProxy()
        flux["default"]      = "no"
        flux[p_rgh[version]] = ""
        flux["pcorr"]        = ""
        flux[alpha[version]] = ""
        self["fluxRequired"] = flux

if __name__ == "__main__" :
   print(FvSchemes("test" , version = "foamStar" , orthogonalCorrection = "implicit"))
