for item in $(find . -name *.gch)
do
    if [ -f $item ];
    then
        rm $item
        echo $item is removed
    fi
done