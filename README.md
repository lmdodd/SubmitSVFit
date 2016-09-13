# SubmitSVFit
```
cmsrel CMSSW_8_0_14 #for 2016 analysis
cd CMSSW_8_0_14/src/
cmsenv
git cms-addpkg RecoMET/METPUSubtraction
git cms-addpkg DataFormats/METReco
git cms-merge-topic -u rfriese:mvamet80
mkdir RecoMET/METPUSubtraction/data
cd RecoMET/METPUSubtraction/data
wget https://github.com/rfriese/cmssw/raw/MVAMET2_beta_0.6/RecoMET/METPUSubtraction/data/weightfile.root
cd $CMSSW_BASE/src
git clone --recursive https://github.com/truggles/SubmitSVFit.git
cd SubmitSVFit
source recipe.sh
scram b -j 8
```

There are a number of current versions of the stand alone svFitter included here.
You can find them in ROOT/bin/SVFitStandAlone... with the names of their executables
defined ROOT/bin/BuildFile.xml.

Most current on being worked on: 
ROOT/bin/SVFitStandAloneFSA.cc

To run in interactive mode for example:
```
SVFitStandAloneFSA inputFile=coolInputFile.root newOutputFile=1 newFile=tmpOut.root doES=1
```

 - inputFile = obvious
 - newOutputFile = 0/1
   - 0 = update input file with svFit vars
   - 1 = output new file with original TTree and new svFit vars
 - newFile = name of output file, default is newFile.root if none specified
 - doES = apply energy scale adjustments providing nominal, shift UP and shift DOWN values of svFit
   - 0 = default, no shift
   - 1 = apply shifts
 - recoilType = type of recoile correction based on generator, MadGraph or amc@nlo
   - 0 = no recoil corrections (for all non-DYJets/WJets/Higgs samples)
   - 1 = amc@nlo recoil corrections
   - 2 = MadGraph recoil corrections
 - isWJets = this shifts the number of jets used in recoil corrections, it is critical for
WJets samples because we clear our jets to preven overlapping with out leptons, but
with WJets one of the leptons is a jet
   - 0 = non-WJets samples
   - 1 = WJets sample
 - metType = MVA MET vs. PF MET
   - 1 = Mva Met
   - -1 = PF Met

To submit jobs to condor:
```
cd test
python svFitSubmitter.py -dr -sd /hdfs/store/user/truggles/svFitTestSept03 -es=1 -r=2 -iswj=0 -mt=-1 --jobName svFitSept03Test
```

 - -dr = dryRun and outputs a command for you to run
 - -sd = select directory, the input directory, this will produce a list of all files in that directory to run over<BR>
       you must have your files in a /hdfs directory.
 - -es = apply energy scale, see above
 - --jobName = applys a job name to append to your file names and output directory structure
 - -r = recoileType
 - -iswj = isWJets
 - -mt = metType


To get your files from elsewhere to /hdfs do something like this:
```
gsido mkdir /hdfs/store/user/truggles/mySubmitDir
gsido rsync -ahP /nfs_scratch/truggles/httSept04skimMerged/*.root /hdfs/store/user/truggles/httSept04skimMerge/
```

It is VERY helpful to make sure that you have ~1000 events per file when running this on Condor.  Anything much larger will take forever,
especially if you run with Energy Shifts. If using FSA ntuples:
 - skim your ntuples without merging any files
 - do a controlled merge that hadds ~1000 events / output file

To do this controlled merge edit the file tools/controlledMerge.py and specify your
 - original directory
 - samples
 - TTreePath
 - output directory
 - and edit the event count / file if you would like to adjust it away from 1,000
Depending on your file naming convention, you may have to edit line 16<BR>

 
Then:
```
python tools/controlledMerge.py
```


