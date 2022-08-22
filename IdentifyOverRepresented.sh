# Run this from output directory of step 2 (metagenomic analysis)
# Requires list of virus positive samples (see github README)
# to run: bash IdentifyOverRepresented.sh

# Check input file:
inputTrue=`ls -l virusPositives.txt 2>/dev/null | wc -l`
	if [ $inputTrue = 0 ]
	then
		echo "No list of virus positive samples found"; exit 1
	fi

# Begin process
mkdir PotentialHosts
cd PotentialHosts
cp ../virusPositives.txt .
posCount=`cat virusPositives.txt | wc -l`
echo "You found $posCount virus positive samples" >> summary.txt; echo >> summary.txt
while read line; do ln -s ../$line*uniq . ; done < virusPositives.txt

# All taxa found across virus positive samples
cat *uniq | sort > all_taxon_occurences.txt

# Top hits
echo "Of these, the top 25 most frequently found taxa are (sampleCount, taxon):"  >> summary.txt
cat all_taxon_occurences.txt | uniq -c | sort -k1,1nr | head -n 25 | sed 's/^     //' | tr ' ' '\t' >> summary.txt

# Clear up
rm *uniq all_taxon_occurences.txt
