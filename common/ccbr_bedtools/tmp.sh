get_ratio() {
min=$(($1>$2?$2:$1))
max=$(($1>$2?$1:$2))
r=`echo $max $min|awk '{printf("%.2f",$1/$2)}'`
echo "$r"
}
result=`get_ratio $1 $2`
echo $result

