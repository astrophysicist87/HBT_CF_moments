#!/bin/bash

echo 'Started processing at' `date`

time nice -n 0 ./run_HBT_CF_moments &> tmp_wo_collapse.out

echo 'Finished processing and cleaning up at' `date`
