#!/bin/bash

SHARED_FLAGS = --ci 0.1 --pc --gl --ll --outpath "../outputs" --nt 8 --pcd --epmae 2.00 --epmic 0.0075
SHARED_ANN_FLAGS = --ann_steps_wo_mod 8000000 --ann_steps 120 --ann_cooling_rate 0.999 --ann_target_temp 1e-13

../bin/SiMaNo --inputfile "../inputs/hashsieve.txt" --texfile "hashsieve" $SHARED_FLAGS $SHARED_ANN_FLAGS

../bin/SiMaNo --inputfile "../inputs/subsieve.txt" --texfile "subsieve" $SHARED_FLAGS $SHARED_ANN_FLAGS

../bin/SiMaNo --inputfile "../inputs/g6ksieve.txt" --texfile "g6ksieve" $SHARED_FLAGS $SHARED_ANN_FLAGS

../bin/SiMaNo --inputfile "../inputs/pbkz.txt" --texfile "pbkz" $SHARED_FLAGS $SHARED_ANN_FLAGS

../bin/SiMaNo --inputfile "../inputs/enumtheory.txt" --texfile "enumtheory" $SHARED_FLAGS $SHARED_ANN_FLAGS --logy

../bin/SiMaNo --inputfile "../inputs/synt9.txt" --texfile "synt9" $SHARED_FLAGS $SHARED_ANN_FLAGS --epmae 3.19 --epmic 0.0001 --logy