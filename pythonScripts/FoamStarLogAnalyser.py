#!/usr/bin/python
import os
from os.path import join
import re
from PyFoam.LogAnalysis.FoamLogAnalyzer import FoamLogAnalyzer
from PyFoam.LogAnalysis.ExecutionTimeLineAnalyzer import GeneralExecutionLineAnalyzer
from PyFoam.LogAnalysis.ContinuityLineAnalyzer import GeneralContinuityLineAnalyzer
from PyFoam.LogAnalysis.LinearSolverLineAnalyzer import GeneralLinearSolverLineAnalyzer,GeneralLinearSolverIterationsLineAnalyzer
from PyFoam.LogAnalysis.ExecutionTimeLineAnalyzer import GeneralExecutionLineAnalyzer
from PyFoam.LogAnalysis.SimpleLineAnalyzer import GeneralSimpleLineAnalyzer
from PyFoam.LogAnalysis.DeltaTLineAnalyzer import GeneralDeltaTLineAnalyzer


class GeneralExecutionLineAnalyzer_foamStar(GeneralExecutionLineAnalyzer):
   """
   Correct GeneralExecutionLineAnalyzer for FoamStar format
   """

   def __init__(self,*args, **kwargs) :
       GeneralExecutionLineAnalyzer.__init__(self,*args, **kwargs)
       timeExpression = "^ExecutionTime = (.+) s  ClockTime = (.+) s  CurrExecTime = (.+) s \((.+)\)$"
       self.exp=re.compile(timeExpression)
       
class FoamStarLogAnalyser(FoamLogAnalyzer):
    """
    The analyzer for the FoamStar log file
    """
    def __init__(self, progress=False, doTimelines=False, doFiles=True, singleFile=True, startTime=None, endTime=None):
        """
        @param progress: Print time progress on console?
        @param doTimelines: generate timelines?
        @param doFiles: generate files?
        """
        FoamLogAnalyzer.__init__(self,progress=progress)

        self.addAnalyzer("Continuity",  GeneralContinuityLineAnalyzer(doTimelines=doTimelines, doFiles=doFiles, singleFile=singleFile,startTime=startTime,endTime=endTime))
        self.addAnalyzer("Linear",      GeneralLinearSolverLineAnalyzer(doTimelines=doTimelines,doFiles=doFiles,singleFile=singleFile,startTime=startTime,endTime=endTime))
        self.addAnalyzer("Execution",   GeneralExecutionLineAnalyzer_foamStar(doTimelines=doTimelines,doFiles=doFiles,singleFile=singleFile,startTime=startTime,endTime=endTime))
        self.addAnalyzer("DeltaT",      GeneralDeltaTLineAnalyzer(doTimelines=doTimelines,doFiles=doFiles,singleFile=singleFile,startTime=startTime,endTime=endTime))
        courantExpression="^Courant Number mean: (.+) max: (\S+).*$"
        self.addAnalyzer("Courant",GeneralSimpleLineAnalyzer("courant",courantExpression,titles=["mean","max"],doTimelines=doTimelines,doFiles=doFiles,singleFile=singleFile,startTime=startTime,endTime=endTime))
        alphaExpression= "^Phase-1 volume fraction = (.+)  Min\(alpha.water\) = (.+)  Max\(alpha.water\) = (.+)"
        self.addAnalyzer("Alpha",GeneralSimpleLineAnalyzer("alpha",alphaExpression,titles=["mean" , "min" , "max"],doTimelines=doTimelines,doFiles=doFiles,singleFile=singleFile,startTime=startTime,endTime=endTime))
        #self.addAnalyzer("Iterations",GeneralLinearSolverIterationsLineAnalyzer(doTimelines=doTimelines,doFiles=doFiles,singleFile=singleFile,startTime=startTime,endTime=endTime))


class LogAnalysis( object) :
   """
      Convenience class for quick plot of the log file
   """
   def __init__(self, logfile , outputPath = "none", analyse = False):

      self.analyser = FoamStarLogAnalyser()
      self.logfile = os.path.abspath(logfile)
      if outputPath == "none" :
         outputPath = join( os.path.dirname(logfile),os.path.basename(logfile)+"_analyzed")
      self.outputPath = outputPath
      if not os.path.exists(outputPath) : os.mkdir(outputPath)
      self.analyser.setDirectory(outputPath)

      if analyse  : self.analyse()

   def analyse(self) :
      print "Parsing " , self.logfile
      fh=open(self.logfile,'r')
      self.analyser.analyze(fh)


   def plot(self,   listPlot, iteration = False ):
      import numpy as np
      from matplotlib import pyplot as plt
      if not os.path.exists(self.outputPath) : self.analyse()
      fig, ax = plt.subplots()
      ax.grid(True)
      for f, c in listPlot :
         fname = join( self.outputPath , f )
         if os.path.exists(fname) :
            try  :
               data = np.loadtxt( fname  )
            except :
               data = np.loadtxt( fname, skiprows = 2 )
	    with open( fname, "r") as lf : labels = lf.readline().strip().split()[1:]
         else :
            print f , "does not exists"
	    print "Available quantities are : " , os.listdir(self.outputPath)
            return
         if iteration :  #Plot versus iteration
            ax.plot( data[:,c] , label = labels[c] )
	    ax.set_xlabel("Iteration")
         else :          #Plot versus time
            ax.plot( data[:,0] , data[:,c] , label = labels[c]  )
	    ax.set_xlabel("Time (s)")
	 ax.set_title( "{} - {}".format(os.path.basename(f) , labels[c] ) )
      plt.show()


if __name__ == "__main__" :

   import argparse
   parser = argparse.ArgumentParser(description='FoamStar log parser')
   parser.add_argument( '-plot' ,     help='Quantity to plot' , type = str,  default = "none")
   parser.add_argument( '-col' ,      help='Column to plot' , type = int,  default = 1)
   parser.add_argument( '-update',    help='update the log information' , action="store_true")
   parser.add_argument( '-ouputDir',  help='Directory to store parsed information' , type = str,  default = "none")
   parser.add_argument( '-iteration', help='Plot versus iteration (instead of time)' , action="store_true")
   parser.add_argument( "logfile" )
   args = parser.parse_args()
   a = LogAnalysis( args.logfile , args.ouputDir )
   if args.update : a.analyse()
   if args.plot != "none" : a.plot(  [(args.plot , args.col)] , iteration = args.iteration)
