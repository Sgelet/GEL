l=${3:-10}
k=$((l*$1/$2))
echo -n "["
for ((i=0 ; i<k; i++)) do echo -n "#"; done
for ((j=i ; j<l ; j++)) do echo -n " "; done
echo -n "] "
v=$((100*$1/$2))
echo -n "$v %" $'\r'