rem SET SHARED_FLAGS=--ci 0.1 --pc --gl --nll --outpath "../outputs" --nt 8 --pcd --epmae 2.00 --epmic 0.0075
SET SHARED_FLAGS= --ci 0.1 --logy --gl --ll --outpath "../outputs/n" --nt 8  --epmae 2.50 --epmic 0.0075
SET SHARED_ANN_FLAGS=--ann_steps_wo_mod 100000 --ann_steps 50 --ann_cooling_rate 0.999 --ann_target_temp 1e-11
rem  --pc --pcd

SET COSTM=rss
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "hashsieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/hashsieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "subsieve%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/subsieve.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "g6ksieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/g6ksieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "pbkz%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/pbkz.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "enumtheory%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/enumtheory.txt" --logy %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "synt9%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/synt9.txt" --logy  %SHARED_FLAGS% %SHARED_ANN_FLAGS% --epmae 3.19 --epmic 0.0001
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "breakRSA%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/breakRSA.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

SET COSTM=rarsd
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "hashsieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/hashsieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "subsieve%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/subsieve.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "g6ksieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/g6ksieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "pbkz%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/pbkz.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "enumtheory%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/enumtheory.txt" --logy %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "synt9%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/synt9.txt" --logy  %SHARED_FLAGS% %SHARED_ANN_FLAGS% --epmae 3.19 --epmic 0.0001
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "breakRSA%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/breakRSA.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

SET COSTM=r2
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "hashsieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/hashsieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "subsieve%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/subsieve.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "g6ksieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/g6ksieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "pbkz%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/pbkz.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "enumtheory%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/enumtheory.txt" --logy %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "synt9%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/synt9.txt" --logy  %SHARED_FLAGS% %SHARED_ANN_FLAGS% --epmae 3.19 --epmic 0.0001
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "breakRSA%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/breakRSA.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

SET COSTM=rmse
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "hashsieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/hashsieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "subsieve%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/subsieve.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "g6ksieve%COSTM%" --st "exponentialsolution" --inputfile "../inputs/g6ksieve.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "pbkz%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/pbkz.txt"  %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "enumtheory%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/enumtheory.txt" --logy %SHARED_FLAGS% %SHARED_ANN_FLAGS%
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "synt9%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/synt9.txt" --logy  %SHARED_FLAGS% %SHARED_ANN_FLAGS% --epmae 3.19 --epmic 0.0001
start "" /wait ../bin/SimAnMoDriver.exe --cct "%COSTM%cost" --texfile "breakRSA%COSTM%" --st "exponentialsolution"  --inputfile "../inputs/breakRSA.txt" %SHARED_FLAGS% %SHARED_ANN_FLAGS%

PAUSE