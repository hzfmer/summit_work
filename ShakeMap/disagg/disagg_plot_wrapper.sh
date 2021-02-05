#!/bin/bash

#./run_opensha.sh -Xmx12g org.opensha.sha.cybershake.plot.DisaggregationPlotter $@
./run_opensha.sh -Dcybershake.db.host=updated_study_15_12.sqlite -Xmx12g org.opensha.sha.cybershake.plot.DisaggregationPlotter $@
