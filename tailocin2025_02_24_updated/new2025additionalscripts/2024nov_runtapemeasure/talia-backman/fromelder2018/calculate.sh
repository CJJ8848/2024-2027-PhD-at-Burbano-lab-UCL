less metadata_1524_OTU51only.txt | grep 'OTU1' -w | awk -F, '{print $1}' | sed 's/plate/p/g' | sort > 1355_OTU1_OTU5names.txt
cat 1355OTU5.txt 1355_OTU1_OTU5names.txt| sort | uniq -c | grep '1 '