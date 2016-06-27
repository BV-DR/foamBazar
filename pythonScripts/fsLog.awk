#! /usr/bin/awk -f

#########################################################################
# Filename: fsLog.awk                                                   #
# Date:     2016-June-22                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:  This awk-script prepares data in log-files for further      #
#           analyses. Data for each log file are kept in sub-folders    #
#           ./fsLog_<log-file>/   Data from stdin are kept in ./fsLog/  #
#           WARNING: files in existing sub-folders will be overwritten. #
#########################################################################

BEGIN {
    nSteps=0
    Time=0
    piping=0
}

END {
    newLineInFiles()
}

# create folder for each log-file
FNR==1 {
    # default folder for data from stdin
    logdir="./fsLog/"
    if (FILENAME!="-")
    {
        print "fsLog: Processing log-file: " FILENAME
        name=FILENAME
        sub(".*/", "", name)
        logdir="./fsLog_" name "/"
    }
    else
    {
        # flush at every new line to when piping from stdin
        piping=1
    }
    system("mkdir -p \"" logdir "\"");
}

function resetCounters()
{
    for (varName in counters) { counters[varName]=-1 }
}

function newLineInFiles()
{
    for (varName in files)
    {
        if (hasNewLines[varName])
        {
            print " " >> files[varName]
        }
        hasNewLines[varName]=0
    }
    if (piping)
    {
        for (varName in files) { fflush(files[varName]) }
    }
}

function checkFile(name, header)
{
    if (length(files[name]) == 0)
    {
        files[name]=(logdir name)
        counters[name]=-1
        print header > files[name]
    }
    if (counters[name]!=nSteps)
    {
        printf Time >> files[name]
        counters[name]=nSteps
    }
    hasNewLines[name]=1
}

# Extract value after columnSel
function extract(line,columnSel,outVar,a,b)
{
    a=index(line, columnSel)
    b=length(columnSel)
    split(substr(line, a+b),outVar)
    gsub("[,:()]","",outVar[1])
}

# Iteration separator
/^Time = / {
    nSteps++
    resetCounters()
    #
    extract($0, "Time = ", val)
    Time=val[1]
    name="Time"
    checkFile(name, "#Time: <value>")
    #
    newLineInFiles()
}

#Interface Courant Number mean: <value> max: <value>
/^Interface Courant Number / {
    name = "Courant_mean_interface"
    checkFile(name, "#Interface Co-Number (mean)")
    extract($0, "mean: ", val)
    printf "\t" val[1] >> files[name]
    name = "Courant_max_interface"
    checkFile(name, "#Interface Co-Number (max.)")
    extract($0, "max: ", val)
    printf "\t" val[1] >> files[name]
}

#Courant Number mean: <value> max: <value>
/^Courant Number / {
    name = "Courant_mean"
    checkFile(name, "#Co-Number (mean)")
    extract($0, "mean: ", val)
    printf "\t" val[1] >> files[name]
    name = "Courant_max"
    checkFile(name, "#Co-Number (max.)")
    extract($0, "max: ", val)
    printf "\t" val[1] >> files[name]
}

#PIMPLE: iteration <value>
/^PIMPLE: iteration / {
    name = "nIter_PIMPLE"
    checkFile(name, "#PIMPLE iter.")
    extract($0, "iteration ", val)
    printf "\t" val[1] >> files[name]
}

#fsi: <iter> residual: <value> (target: <value>)
/^fsi: / {
    name = "Res_fsi"
    checkFile(name, "#Res. fsi")
    extract($0, "residual: ", val)
    printf "\t" val[1] >> files[name]
}

#Phase-1 volume fraction = <value> Min(alpha.water) = <value>  Max(alpha.water) = <value>
/^Phase-1 volume fraction = / {
    name = "Phase_volume"
    checkFile(name, "#Phase vol.fraction")
    extract($0, "volume fraction = ", val)
    printf "\t" val[1] >> files[name]
    #
    name = "Phase_min"
    checkFile(name, "#alpha (min)")
    extract($0, "Min(alpha.water) = ", val)
    printf "\t" val[1] >> files[name]
    #
    name = "Phase_max"
    checkFile(name, "#alpha (max-1.0)")
    extract($0, "Max(alpha.water) = ", val)
    printf "\t" (val[1]-1.0) >> files[name]
}

# Extraction of any solved for variable
/:[ \t]*Solving for / {
    extract($0, "Solving for ", varNameVal)
    # no underscore in var. name due to restriction(s) in fsPlot.py
    gsub("_","",varNameVal[1])
    header="#Res. " varNameVal[1] " "
    #
    name="Res_init_" varNameVal[1]
    checkFile(name, header "(init.)")
    extract($0, "Initial residual = ", val)
    printf "\t" val[1] >> files[name]
    #
    name="Res_final_" varNameVal[1]
    checkFile(name, header "(final)")
    extract($0, "Final residual = ", val)
    printf "\t" val[1] >> files[name]
    #
    name="nIter_" varNameVal[1]
    checkFile(name, "#nIter. " varNameVal[1])
    extract($0, "No Iterations ", val)
    printf "\t" val[1] >> files[name]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

#time step continuity errors : sum local = <value>, global = <value>, cumulative = <value>
/^time step continuity errors :/ {
    name="contErr_local"
    checkFile(name, "#contErr (local)")
    extract($0, "sum local = ", val)
    printf "\t" val[1] >> files[name]
    #
    name="contErr_global"
    checkFile(name, "#contErr (global)")
    extract($0, "global = ", val)
    printf "\t" val[1] >> files[name]
    #
    name="contErr_cumu"
    checkFile(name, "#contErr (cumul.)")
    extract($0, "cumulative = ", val)
    printf "\t" val[1] >> files[name]
}

#ExecutionTime = <value> s  ClockTime = <value> s  CurrExecTime = <value> s (<value>)
/^ExecutionTime = / {
    name="timing_exec"
    checkFile(name, "#Exec. time (cumul.)")
    extract($0, "ExecutionTime = ", val)
    printf "\t" val[1] >> files[name]
    #
    name="timing_clock"
    checkFile(name, "#Clck. time")
    extract($0, "ClockTime = ", val)
    printf "\t" val[1] >> files[name]
}

#
/^ExecutionTime = .* CurrExecTime = / {
    name="timing_curr"
    checkFile(name, "#Exec. time (curr. step)")
    extract($0, "CurrExecTime = ", val)
    printf "\t" val[1] >> files[name]
}

#Execution time for mesh.update() = <value> s
/^Execution time for mesh/ {
    name="timing_meshUpdate"
    checkFile(name, "#Exec. time (mesh update)")
    extract($0, "mesh.update() = ", val)
    printf "\t" val[1] >> files[name]
}

#fluidForce: relax (f,m) : <value> <value> (<value value value>) (<value value value>)
/^fluidForce: relax \(f,m\) : / {
    header="#force "
    extract($0, "relax (f,m) : ", val)
    for (i in val) { gsub("[()]","", val[i]) }
    name="fluidForce_relax"
    checkFile(name, header "relaxCoeff.")
    printf "\t" val[1] >> files[name]
    #
    name="fluidForce_x"
    checkFile(name, header "fx")
    printf "\t" val[3] >> files[name]
    #
    name="fluidForce_y"
    checkFile(name, header "fy")
    printf "\t" val[4] >> files[name]
    #
    name="fluidForce_z"
    checkFile(name, header "fz")
    printf "\t" val[5] >> files[name]
    #
    header="#moment "
    #
    name="fluidMoment_relax"
    checkFile(name, header "relaxCoeff")
    printf "\t" val[2] >> files[name]
    #
    name="fluidMoment_x"
    checkFile(name, header "mx")
    printf "\t" val[6] >> files[name]
    #
    name="fluidMoment_y"
    checkFile(name, header "my")
    printf "\t" val[7] >> files[name]
    #
    name="fluidMoment_z"
    checkFile(name, header "mz")
    printf "\t" val[8] >> files[name]
}
