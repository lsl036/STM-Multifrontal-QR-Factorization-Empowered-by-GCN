#!/bin/bash

# 测试论文数据计算时间
IFS_old=$IFS
IFS=$'\n'
for line in $(<test.txt)
do
	#echo "$line"
	# arr=($line) 
	# echo "${arr[0]}:${arr[1]}"
	var1=`echo $line|awk -F ' ' '{print $1}'`
	# echo $var1
	var2=`echo $line|awk -F ' ' '{print $2}'`
	# echo $var2
	# var3=`echo $line|awk -F ' ' '{print $3}'`
	./qrtest $var1 $var2 
done;
IFS=$IFS_old