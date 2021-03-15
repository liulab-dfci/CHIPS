#!/usr/bin/env python
"""This script parses a set of _conservation.R and fills that into a
template (conserv_plotConserv.R.txt) to generate an R script 
(conserv_plotConserv.R) that, when run, generates the plot (as a png)
"""
import os
import sys
from optparse import OptionParser
from string import Template

def parseR(rscript):
    """Given a {run}_conservation.R, parse out the X,Y data--returns a tuple
    """
    f = open(rscript)
    #SKIP lines 1-3
    tmp = f.readline().strip()
    tmp = f.readline().strip()
    tmp = f.readline().strip()

    #LINE 4: READ in x -- ignored for now
    tmp = f.readline().strip()
    x = eval(tmp.replace("x<-c",""))

    #READ in y -- ignored for now
    tmp = f.readline().strip()
    y = eval(tmp.replace("y0<-c",""))

    #print(len(x), len(y))
    return (x,y)
    
def dataToString(var, data):
    """Given a tuple of data, and a name to save it as
    returns var <- c(data)
    """
    #convert data to strings
    d = [str(d) for d in data]
    return "%s <- c(%s)" % (var, ",".join(d))

def generateOneRun(*args):
    """Given a list (of 3) RUNNAMES extracted from the R model, generate
    the R-code/string to make a plot for the (3) runs
    
    NOTE: I know that this breaks everything I believe about separation of 
    conten/code, but oh well...
    """
    #UGLY HACK--need to handle case of 1, case of 2, case of 3
    #this should be in afile
    tmp1 = Template(
    """
    ymax<-max($r1)
    ymin<-min($r1)
    yquart<-(ymax-ymin)/4
    plot(x,$r1,type='l',xlab="Distance from Center (bp)", ylab="Avg. Phastcons",col=icolor[1], ylim=c(ymin-yquart, ymax+yquart), lwd=2.0)
    legend('topright', c("$r1"), lty=1,col=icolor, cex=1.0)
    """)

    tmp2 = Template(
    """
    ymax<-max($r1,$r2)
    ymin<-min($r1,$r2)
    yquart<-(ymax-ymin)/4
    plot(x,$r1,type='l',xlab="Distance from Center (bp)", ylab="Avg. Phastcons",col=icolor[1], ylim=c(ymin-yquart, ymax+yquart), lwd=2.0)
    lines(x,$r2,type='l',col=icolor[2], lwd=2.0)
    legend('topright', c("$r1","$r2"), lty=1,col=icolor, cex=1.0)
    """)

    tmp3 = Template(
    """
    ymax<-max($r1,$r2,$r3)
    ymin<-min($r1,$r2,$r3)
    yquart<-(ymax-ymin)/4
    plot(x,$r1,type='l',xlab="Distance from Center (bp)", ylab="Avg. Phastcons",col=icolor[1], ylim=c(ymin-yquart, ymax+yquart), lwd=2.0)
    lines(x,$r2,type='l',col=icolor[2], lwd=2.0)
    lines(x,$r3,type='l',col=icolor[3], lwd=2.0)
    legend('topright', c("$r1","$r2","$r3"), lty=1,col=icolor, cex=1.0)
    """)
    if len(args) >= 3:
        return tmp3.substitute(r1=args[0],r2=args[1],r3=args[2])
    elif len(args) >= 2:
        return tmp2.substitute(r1=args[0], r2=args[1])
    else:
        return tmp1.substitute(r1=args[0])

def main():
    usage = "USAGE: %prog -r [R file 1] -r [R file 2] ...-r [R file N] -t [path to conserv_plotConserv.R] -o [output.png]; OUTPUTS- conserv_plotConserv.R"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--rfiles", action="append", help="list of files")
    optparser.add_option("-t", "--template", help="conserv_plotConserv.R file path")
    optparser.add_option("-p", "--pngout", help="path to output pngs: conservationPlot01.png, conservationPlot02.png, etc. NOTE: no trailing slashes please!")
    optparser.add_option("-o", "--output", help="conserv_plotConserv.R output", default="./conserv_plotConserv.R")
    (options, args) = optparser.parse_args(sys.argv)

    if (not options.rfiles or not options.template or not options.pngout):
        optparser.print_help()
        sys.exit(-1)

    #GET run names
    runs = [r.split("/")[-1].replace("_conserv.R", "") for r in options.rfiles]
    #PARSE the data from the runs
    data = { runName: parseR(rscript) for (runName, rscript) in zip (runs,options.rfiles) }
    #compose the data string: RUN <- c(...YDATA...)
    data_str = "\n".join([dataToString(r, data[r][1]) for r in runs])
    
    #runNames = c(Run1, Run2,...)
    runNames = "c(%s)" % ",".join(runs)

    #for x labels, try to take the first elm's x-s
    x = [str(d) for d in data[runs[0]][0]]
    xvals = "c(%s)" % ",".join(x)

    #handle the plotting
    _nPerPlot=3
    #SPLIT into runs into groups of _nPerPlot (i.e. 3)
    groups = [runs[i:i+_nPerPlot] for i in range(0,len(runs),_nPerPlot)]
    plots = "\n".join([generateOneRun(*g) for g in groups])

    #PARSE the data
    #read in template
    f = open(options.template)
    tmp = Template(f.read())
    out = tmp.substitute(data=data_str, runNames=runNames, xvals=xvals, 
                         pngout=options.pngout, plots=plots)
    f.close()

    #WRITE output conserv_plotConserv.R
    f = open(options.output, "w")
    f.write(out)
    f.close()

if __name__=='__main__':
    main()
