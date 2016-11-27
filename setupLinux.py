#!/usr/bin/env python
#Setting up phyloGenerator's requires folder
#Will Pearse

print "Linux configuration script for phyloGenerator\n"
print "Pass, as additional command line arguments, the path where you've downloaded all your files"
print "\tand *which folder contains BEAST*"
print "\te.g., './setupLinux.py /home/will/phyloGenerator /home/will/phyloGenerator/BEAST\ v1.7.4'"
print "Make sure all programs are executable - if unsure, make them so"
print "\te.g., 'chmod +x NAMEOFPROGRAM'"
print "The resulting 'requires' folder must contain only output from this script"
print "Do not leave source code from the programs phyloGenerator uses in the same folder"
print "\tas phyloGenerator.py, or (for safety) in your output 'working directory'"
print "\tThis will cause obscure-looking errors from phyloGenerator!"

#Modules
import os, sys, subprocess

#Functions
def checkProgram(programs, programName, paths, requiresPath, errorMessage="", args=[]):
    for program in programs:
        for path in paths:
            try:
                call = [path + program]
                call.extend(args)
                subprocess.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocess.Popen(['ln', '-s', path + program, requiresPath + programName], stdout=subprocess.PIPE)
                return True
            except:
                pass
    else:
        print "!!!", programName, "could not be detected or setup"
        print "!!!", errorMessage, "\n"
        return False


def checkModule(module, errorMessage=""):
    try:
        __import__(module)
        return True
    except:
        print "...", module, "is not setup"
        print "!!!", errorMessage, "\n"
        return False

#Handling paths
paths = os.environ['PATH'].split(":")
for i,path in enumerate(paths):
    if path[-1] != "/":
        paths[i] += "/"

if len(sys.argv) > 1:
    for path in sys.argv:
        if path[-1] != "/":
            path += "/"
        paths.append(path)

#Setup the requires folder
subprocess.Popen(['rm', '-rf', 'requires'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.Popen(['mkdir', 'requires'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)


requiresPath = os.getcwd() + "/requires/"

#External programs
print "Checking and configuring external programs\n"
tracker = True
tracker = checkProgram(['mafft', 'MAFFT'], 'mafft', paths, requiresPath, "Run 'sudo apt-get install mafft'", ['--version']) and tracker
tracker = checkProgram(['muscle', 'MUSCLE'], 'muscle',  paths, requiresPath, "Run 'sudo apt-get install muscle'") and tracker and tracker
tracker = checkProgram(['prank', 'PRANK', 'prank.1'], 'prank',  paths, requiresPath, "Download from here: 'http://code.google.com/p/prank-msa/downloads/list'") and tracker
tracker = checkProgram(['clustalo', 'clustal-o', 'clustalomega', 'clustal-omega', 'clustalo-1.1.0-linux-64', 'clustalo-1.0.3-Ubuntu-i686'], 'clustalo', paths, requiresPath, "Download from here: 'http://www.clustal.org/omega/#Download'") and tracker
tracker = checkProgram(['metal'], 'metal', paths, requiresPath, "Download from here: 'http://kumiho.smith.man.ac.uk/blog/whelanlab/?page_id=396'") and tracker
tracker = checkProgram(['trimal'], 'trimal', paths, requiresPath, "Download from here: 'http://trimal.cgenomics.org/downloads' and compile inside the folder 'source' with 'make'") and tracker
tracker = checkProgram(['phylomatic'], 'phylomatic', paths, requiresPath, "Download from here: 'https://github.com/phylocom/phylocom/downloads' and compile inside the folder 'src' with 'make'") and tracker
tracker = checkProgram(['pathd8', 'PATHd8', 'PATHD8'], 'PATHd8', paths, requiresPath, "Download from here: 'http://www2.math.su.se/PATHd8/compile.html' and compile with 'cc PATHd8.c -O3 -lm -o PATHd8'") and tracker
tracker = checkProgram(['RAxML', 'raxml', 'raxmlHPC-SSE3', 'raxmlHPC'], 'raxml',  paths, requiresPath, "Download from here: 'https://github.com/stamatak/standard-RAxML' and compile the sequential version ('make -f Makefile.gcc')") and tracker
#BEAST requires some additional trickery...
beastPaths = [x+'bin/' for x in paths]
beastPaths.extend(paths)
tracker = checkProgram(['BEAST', 'beast'], 'beast',  beastPaths, requiresPath, "Install Java ('sudo apt-get install default-jre')\n!!! Download from here: 'http://beast.bio.ed.ac.uk/Main_Page'", ["--help"]) and tracker
tracker = checkProgram(['treeannotator', 'TREEANNOTATOR'], 'treeannotator',  beastPaths, requiresPath, "Download from here: 'http://beast.bio.ed.ac.uk/Main_Page'") and tracker

#Python module
print "Checking Python libraries\n"
tracker = checkModule("numpy", "Run 'sudo apt-get install python-numpy'") and tracker
tracker = checkModule("scipy", "Run 'sudo apt-get install python-scipy'") and tracker
tracker = checkModule("Bio", "Run 'sudo apt-get install python-biopython'") and tracker
tracker = checkModule("dendropy", "Run 'sudo pip install -U dendropy'") and tracker

#Clean-up
if tracker:
    print "\n\nCONGRATULATIONS!"
    print "phyloGenerator is setup. You should now be able to run it by typing './phyloGenerator.py'"
else:
    print "\n\nCOMMISERATIONS!"
    print "phyloGenerator is not setup. Above should be a list of instructions to follow; do this, then re-run this script"

