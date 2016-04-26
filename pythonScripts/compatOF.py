alpha = { "foamStar" : "alpha.water" ,
          "foamExtend" : "alpha1"
        }


p_rgh = { "foamStar" : "p_rgh" ,
          "foamExtend" : "pd"
        }


application = { "foamStar" : "CRSinterDyMFoam" ,
                "foamExtend" : "navalFoam" ,
              }
              
surfaceElevation = { "foamStar" : "waveProbes" ,
                      "foamExtend" : "surfaceElevation" ,
                    }
                    
initWaveField = { "foamStar" : "CRSinitWaveField" ,
                  "foamExtend" : "initWaveField" ,
                }
                
waveTypeDict = { "stokes5th" : {"foamStar" : "stokes5th" ,
                                "foamExtend" : "stokesFifth" ,
                               },
                 "noWaves"  : {"foamStar" : "noWaves" ,
                                "foamExtend" : "noWaves" ,
                               },
                 "streamFunction"  : {"foamStar" : "streamFunction" ,
                                "foamExtend" : "streamFunction" ,
                               }
               }
               
water = {
            "foamStar" : "water",
            "foamExtend" : "phase1",
         }

air = {
            "foamStar" : "phase2|air",
            "foamExtend" : "phase2",
         }

