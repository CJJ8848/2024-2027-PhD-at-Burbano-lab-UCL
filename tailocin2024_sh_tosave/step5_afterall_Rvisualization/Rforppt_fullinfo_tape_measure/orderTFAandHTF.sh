#change TFA to TFA
TFA=allTFAnames.txt
TFAout=TFAorder.txt
input=treetopologyorder.txt
>"$TFAout"
while read -r name; do
    hp=$(less $TFA | grep $name | cut -d '|' -f2) 
    echo 'found' $name $hp
    echo  $name $hp >> $TFAout
done < $input
