for f in $(find . -name '*.shader' );
do
    sed -i 's/    /\t/g' $f
    echo $f
done

