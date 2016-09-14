if [ $# -ne 2 ]; then
		echo ""
		echo "   Usage $0 [command file] [num threads]"
		echo "      [command file] - A list of commands one per line."
    echo "      [num threads ] - The max number of threads to spawn at a time."
		echo ""
		exit;
fi


inputFile=$1
numThreads=$2

cat $inputFile | xargs -I CMD --max-procs=$numThreads bash -c CMD
