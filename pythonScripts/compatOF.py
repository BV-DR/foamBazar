alpha = { "foamStar" : "alpha.water" ,
          "foamExtend" : "alpha1" , 
          "swenseFoam" : "alpha1"
        }


p_rgh = { "foamStar" : "p_rgh" ,
          "foamExtend" : "pd" ,
          "swenseFoam" : "pd"
        }


application = {
               "foamStar" : "CRSinterDyMFoam" ,
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
            "foamStar" : "water",
            "foamExtend" : "phase1",
            "swenseFoam" : "phase1",
         }

air = {
            "foamStar" : "phase2|air",
            "foamExtend" : "phase2",
            "swenseFoam" : "phase2",
         }

