
#
# Script to run DumpPostFitHistos on a single file on the batch
#
# This script can be tested with the following command:
#   $ source run_TruthPlot.sh input/ws.root batch/output/out.root job_test
# and submitted to lxbatch with submit_lxbatch.py.
#
# davide.gerbaudo@gmail.com
# Jul 2015

INPUT_FILE=$1
OUTPUT_FILE=$2
JOBID=$3
# # for /f "tokens=3,* delims= " %%a in ("%*") do set ALL_BUT_FIRST_THREE=%%b
# for /f "usebackq tokens=3*" %%i in (`echo %*`) DO @ set ALL_BUT_FIRST_THREE=%%j
# OTHER_OPTIONS=%ALL_BUT_FIRST_THREE%
OTHER_OPTION=$4

# set -e # exit on error
# set -u # exit on undefined variable

echo "Using these options:"
echo "INPUT_FILE       $INPUT_FILE     "
echo "OUTPUT_FILE      $OUTPUT_FILE    "
echo "JOBID            $JOBID          "
echo "OTHER_OPTION     $OTHER_OPTION   "


BASE_DIR='/afs/cern.ch/user/g/gerbaudo/work/public/hlfv/hlfv_fit_tests/print_postfit_yields'


IN_LOCAL_DIR=`dirname ${INPUT_FILE}`
IN_REMOTE_DIR="${BASE_DIR}/${IN_LOCAL_DIR}"

OUT_LOCAL_DIR=`dirname ${OUTPUT_FILE}`
OUT_REMOTE_DIR="${BASE_DIR}/${OUT_LOCAL_DIR}"

WORK_DIR=$(pwd)
echo "Working in ${WORK_DIR} with the following options"
# echo "LSF Job ID : ${LSF_JOBID}" # undefined?
date

echo "Using the input file ${INPUT_FILE}"
echo "Setting up root"
source ${BASE_DIR}/script/setup.sh
echo "Using root `root-config --version` from `which root`"
cd ${WORK_DIR}

mkdir -p ${IN_LOCAL_DIR}
cp -p ${BASE_DIR}/${INPUT_FILE} ${IN_LOCAL_DIR}/
mkdir -p ${OUT_LOCAL_DIR}

${BASE_DIR}/DumpPostFitHistos -i ${INPUT_FILE} -o ${OUTPUT_FILE} --transport-cov ${OTHER_OPTION} 2>&1 | tee ${OUTPUT_FILE/root/log}


echo "Done `date`"
echo "Local files:"
ls -ltrh
echo "Local output:"
ls -ltrh ${OUT_LOCAL_DIR}

mkdir -p ${OUT_REMOTE_DIR}
cp -p ${OUT_LOCAL_DIR}/* ${OUT_REMOTE_DIR}/
# no need to cleanup, lsf will do it
echo "Done (`date`)"
