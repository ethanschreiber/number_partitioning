# ------------------------------------------
# Make sure this is run from BinPacking dir
# ------------------------------------------
PWD=`pwd`
case $PWD in
  *BinPacking ) 
  ;;
  *) echo "Must run from the BinPacking Directory: $0"
  exit
  ;;    
esac
