prog=muvttmmclj
$FEASST_INSTALL_DIR_/tools/compile.sh $prog
cp $FEASST_INSTALL_DIR_/src/main $prog
export OMP_NUM_THREADS=4
mkdir -p tmp #directory for checkpoint files
./$prog
