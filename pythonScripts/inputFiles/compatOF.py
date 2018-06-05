alpha = { "foamStar" : "alpha.water" ,
          "foamExtend" : "alpha1" , 
          "swenseFoam" : "alpha1"
        }


p_rgh = { "foamStar" : "p_rgh" ,
          "foamExtend" : "pd" ,
          "swenseFoam" : "pd"
        }


application = {
               "snappyHexMesh" : "snappyHexMesh",
               "foamStar"   : "foamStar" ,
               "foamExtend" : "navalFoam" ,
               "swenseFoam" : "swenseFoam"
              }

waveAlpha = {
               "foamStar"   : "waveAlpha" ,
               "foamExtend" : "waveAlpha" ,
               "swenseFoam" : "waveAlphaSwense" ,
            }

waveVelocity = {
   "foamStar"   : "waveVelocity" ,
   "foamExtend" : "waveVelocity" ,
   "swenseFoam" : "waveVelocitySwense" ,
}
              
surfaceElevation = { "foamStar" : "waveProbes" ,
                     "foamExtend" : "surfaceElevation" ,
                     "swenseFoam" : "surfaceElevation" ,
                    }
                    
pointDisp = {
              "foamStar"   : "pointDisplacement" ,
            }
                    
initWaveField = { "foamStar" : "CRSinitWaveField" ,
                  "foamExtend" : "initWaveField" ,
                  "swenseFoam" : "initWaveField" ,
                }
                
waveTypeDict = { "stokes5th" : {"foamStar" : "stokes5th" ,
                                "foamExtend" : "stokesFifth" ,
                                "swenseFoam" : "stokesFifth" ,
                               },
                 "noWaves"  : {"foamStar" : "noWaves" ,
                               "foamExtend" : "noWaves" ,
                               "swenseFoam" : "noWaves" ,
                               },
                 "streamFunction"  : {"foamStar" : "streamFunction" ,
                                "foamExtend" : "streamFunction" ,
                                "swenseFoam" : "streamFunction" ,
                               }
               }
               
water = {
            "foamStar" : "water|phase1",
            "foamExtend" : "phase1",
            "swenseFoam" : "phase1",
         }

air = {
            "foamStar" : "air|phase2",
            "foamExtend" : "phase2",
            "swenseFoam" : "phase2",
         }
         
         
foamStarPatch = {
                    "outlet" : "domainX0",
                    "inlet" : "domainX1",
                    "side" : "domainY1",
                    "side1" : "domainY0",
                    "side2" : "domainY1",
                    "bottom" : "domainZ0" ,
                    "top" : "domainZ1",
                    "structure" : "ship",
                    "auto0" : "domainX0",
                    "auto1" : "domainX1",
                    "auto2" : "domainY0",
                    "auto3" : "domainY1",
                    "auto4" : "domainZ0" ,
                    "auto5" : "domainZ1"
                }
                
defaultPatch = {
                    "inlet" : "inlet" ,
                    "outlet" : "outlet" ,
                    "side" : "side" ,
                    "symmetryPlane" : "symmetryPlane" ,
                    "bottom" : "bottom" ,
                    "top" : "top" ,
                    "structure" : "structure" ,
                    "frontBack" : "frontBack"   #In case of  2D computation
               }

namePatch = { 
                "default" : defaultPatch,
                "foamStar" : foamStarPatch
            }
 
