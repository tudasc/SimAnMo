rem SET SHARED_FLAGS=--ci 0.1 --pc --gl --nll --outpath "../outputs" --nt 8 --pcd --epmae 2.00 --epmic 0.0075
SET SHARED_FLAGS=--gr --ci 0.1 --pc  --logy --gl --ll --outpath "../outputs" --nt 8 --pcd --epmae 2.00 --epmic 0.0075
SET SHARED_ANN_FLAGS=--ann_steps_wo_mod 100000 --ann_steps 80 --ann_cooling_rate 0.999 --ann_target_temp 1e-13

start "" /wait ../bin/SimAnMoDriver.exe --cct "nnrrsscostcalculator" --st "exponentialsolution" --inputfile "../inputs/hashsieve.txt" --texfile "hashsieve" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

start "" /wait ../bin/SimAnMoDriver.exe --cct "nnrrsscostcalculator" --st "exponentialsolution"  --inputfile "../inputs/subsieve.txt" --texfile "subsieve" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

start "" /wait ../bin/SimAnMoDriver.exe --cct "nnrrsscostcalculator" --st "exponentialsolution"  --inputfile "../inputs/g6ksieve.txt" --texfile "g6ksieve" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

start "" /wait ../bin/SimAnMoDriver.exe --cct "nnrrsscostcalculator" --st "exponentialsolution"  --inputfile "../inputs/pbkz.txt" --texfile "pbkz" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

start "" /wait ../x64/Release/SimAnMoDriver.exe --cct "nnrrsscostcalculator" --st "exponentialsolution"  --inputfile "../inputs/enumtheory.txt" --texfile "enumtheory" --logy %SHARED_FLAGS% %SHARED_ANN_FLAGS%

start "" /wait ../x64/Release/SimAnMoDriver.exe --cct "nnrrsscostcalculator" --st "exponentialsolution"  --inputfile "../inputs/synt9.txt" --texfile "synt9" --logy  %SHARED_FLAGS% %SHARED_ANN_FLAGS% --epmae 3.19 --epmic 0.0001

PAUSE