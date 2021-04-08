#!/bin/bash
#cholesky
# 性能测试部分
# IFS_old=$IFS
# IFS=$'\n'
# for line in $(<sym_defmat.txt)
# do
# 	#echo "$line"
# 	# arr=($line) 
# 	# echo "${arr[0]}:${arr[1]}"
# 	var1=`echo $line|awk -F ' ' '{print $1}'`
# 	# echo $var1
# 	var2=`echo $line|awk -F ' ' '{print $2}'`
# 	# echo $var2
# 	./choltest_nesdis $var1 $var2
# done;
# IFS=$IFS_old

# Cholesky 写数据的部分
# IFS_old=$IFS
# IFS=$'\n'
# for line in $(<chol_reorder_select.txt)
# do
# 	#echo "$line"
# 	# arr=($line) 
# 	# echo "${arr[0]}:${arr[1]}"
# 	var1=`echo $line|awk -F ' ' '{print $1}'`
# 	# echo $var1
# 	var2=`echo $line|awk -F ' ' '{print $2}'`
# 	# echo $var2
# 	var3=`echo $line|awk -F ' ' '{print $3}'`
# 	./chol_saveorder $var1 $var2 $var3
# done;
# IFS=$IFS_old


#QR
# 性能测试部分
# IFS_old=$IFS
# IFS=$'\n'
# for line in $(<lunwen_data0301.txt)
# do
# 	#echo "$line"
# 	# arr=($line) 
# 	# echo "${arr[0]}:${arr[1]}"
# 	var1=`echo $line|awk -F ' ' '{print $1}'`
# 	# echo $var1
# 	var2=`echo $line|awk -F ' ' '{print $2}'`
# 	# echo $var2
# 	./qrtest_nesdis $var1 $var2
# done;
# IFS=$IFS_old

# 写数据的部分
# IFS_old=$IFS
# IFS=$'\n'
# for line in $(<lunwen_data.txt)
# do
# 	#echo "$line"
# 	# arr=($line) 
# 	# echo "${arr[0]}:${arr[1]}"
# 	var1=`echo $line|awk -F ' ' '{print $1}'`
# 	# echo $var1
# 	var2=`echo $line|awk -F ' ' '{print $2}'`
# 	# echo $var2
# 	# var3=`echo $line|awk -F ' ' '{print $3}'`
# 	./qrtest_write $var1 $var2 
# done;
# IFS=$IFS_old

# 测试论文数据计算时间
IFS_old=$IFS
IFS=$'\n'
for line in $(<reorder.txt)
do
	#echo "$line"
	# arr=($line) 
	# echo "${arr[0]}:${arr[1]}"
	var1=`echo $line|awk -F ' ' '{print $1}'`
	# echo $var1
	var2=`echo $line|awk -F ' ' '{print $2}'`
	# echo $var2
	# var3=`echo $line|awk -F ' ' '{print $3}'`
	./qrtest_all_time $var1 $var2 
done;
IFS=$IFS_old