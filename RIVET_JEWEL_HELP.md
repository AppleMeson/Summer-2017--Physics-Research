
Cd ..Rivet tutorial
 go inside JEWEL
using jewel to genrate evnts:

open emacs: save bllank file as <filename>.dat
then
ctrl x 2
 open up emacs jewel-2.1.0.f to undeerstand what im changing

in the log file, change the folloowing parameters:

example log file (NOT COMPLETE, THIS IS A VERY BASIC EXAMPLE)

LOGFILE pp2p76TeV.log 
HEPMCFILE pp2p76TeV.hepmc
NEVENT 10000 
NJOB 97676 # random seed 
SQRTS 2760 # energy of collision
PROCESS 'PPJJ' # process type here, it is proton proton jet jet
# Changing the generation phase space
PTMIN 30
PTMAX 500
ETAMAX 3.0

ptmin ptmax and etamax- change these in the log file (.dat file) and 


THEN:

go back to terminal, inside JEWEL,  

./jewel-2.1.0-vac <filename>.dat

TO VIEW GENERATED EVENTS AND RESULTS, IN TERMINAL DO less <filename>.log ( you can see filename by doing ls -alrth in JEWEL

FOR RIVET

after doing everything above, go inside ANALYSIS in JEWEL (in terminal)
type ls

rivet-mkanalysis <filename2>
ls
open filename2.cc using emacs

there, you need to declare the appropriate variables and functions and classes and histograms

rivet-buildplugin RivetMY_FIRST_ANALYSIS.so MY_FIRST_ANALYSIS.cc---→ to compile code and the plot file
ACCESS THE PLOT FILE IN ANALYSIS AND MAKE THE PLOTS OF THE HISTOGRAMS THERE

--
rivet -a MY_FIRST_ANALYSIS ../pp2p76TeV.hepmc -H myanalysis_output_pp2p76tev.yoda --→ runs the code and makes the output YODA file ( which you can name).

rivet-mkhtml --mc-errs --outputdir=plots_mynewanalysis_jetstudies_pp2p76 --font=helvetica myanaly sis_output_pp2p76tev.yoda:'test run' ←----- this creates the plots and uses the options we gave to create it.

You can acces the .dat file for whichever plot you want, and you can change stuff, colours, plot style, LOOK FOR SYNTAX FPR THIS IN make plots-hepforge- online 

everytime you change a plot file, you dont need to redo the build plugin
you can just 
make-plots plotchanged.dat
refresh















cd ~/eos/user/a/adverma/
