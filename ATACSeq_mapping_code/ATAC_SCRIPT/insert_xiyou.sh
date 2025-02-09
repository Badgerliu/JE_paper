#!/usr/bin/env bash

for foo in *.bam
	do
		echo $foo
		bar=$(echo ${foo} | sed 's/.bam//') # define the describer with sed function
		echo $bar  
		#peak calling	
		picard CollectInsertSizeMetrics I=$foo O=${bar}_insert_size_metrics.txt H=${bar}_insert_size_histogram.pdf M=0.5
	done
